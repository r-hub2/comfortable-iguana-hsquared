# Robustness guards for the v0.1 formula parser: each previously leaked a
# cryptic base-R error or silently mis-handled the user's syntax. The fixes live
# in R/model-spec.R; these tests pin the named messages and correct labels.

hs_pr_pedigree <- function() {
  data.frame(
    id = c("a", "b", "c", "d"),
    sire = c(NA, NA, "a", "a"),
    dam = c(NA, NA, "b", "b"),
    stringsAsFactors = FALSE
  )
}

hs_pr_data <- function() {
  data.frame(
    y = c(1, 2, 3, 4),
    sex = c("M", "F", "M", "F"),
    age = c(4, 5, 6, 7),
    id = c("a", "b", "c", "d"),
    stringsAsFactors = FALSE
  )
}

# ---- #1 single-level / zero-row factor fixed effect ----

test_that("single-level factor fixed effect is rejected by name", {
  ped <- hs_pr_pedigree()
  dat <- hs_pr_data()
  dat$sex <- c("M", "M", "M", "M")

  expect_error(
    hsquared:::hs_build_model_spec(
      y ~ sex + animal(1 | id, pedigree = ped),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "fixed-effect term `sex` has only one observed level",
    fixed = TRUE
  )
})

test_that("zero-row model frame is rejected cleanly", {
  ped <- hs_pr_pedigree()
  dat <- hs_pr_data()[0, , drop = FALSE]

  expect_error(
    hsquared:::hs_build_model_spec(
      y ~ sex + animal(1 | id, pedigree = ped),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "fixed-effect model frame has zero rows",
    fixed = TRUE
  )
})

# ---- #2 offset() ----

test_that("offset() in the fixed RHS is rejected, not silently dropped", {
  ped <- hs_pr_pedigree()
  dat <- hs_pr_data()

  expect_error(
    hsquared:::hs_build_model_spec(
      y ~ offset(age) + sex + animal(1 | id, pedigree = ped),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "Offset terms (`offset(age)`) are planned, not implemented",
    fixed = TRUE
  )
})

# ---- #3 univariate response name ----

test_that("transformed univariate response keeps its deparsed label", {
  ped <- hs_pr_pedigree()
  dat <- hs_pr_data()

  spec <- hsquared:::hs_build_model_spec(
    log(y) ~ sex + animal(1 | id, pedigree = ped),
    data = dat,
    family = stats::gaussian(),
    REML = TRUE
  )

  expect_equal(spec$response$name, "log(y)")
  expect_false(spec$response$multivariate)
})

test_that("bare univariate response label is unchanged", {
  ped <- hs_pr_pedigree()
  dat <- hs_pr_data()

  spec <- hsquared:::hs_build_model_spec(
    y ~ sex + animal(1 | id, pedigree = ped),
    data = dat,
    family = stats::gaussian(),
    REML = TRUE
  )

  expect_equal(spec$response$name, "y")
})

# ---- #4 cbind() derived/transformed columns ----

test_that("derived cbind() trait columns are rejected by name", {
  ped <- hs_pr_pedigree()
  dm <- data.frame(
    y1 = c(1, 2, 3, 4),
    y2 = c(2, 3, 4, 5),
    id = c("a", "b", "c", "d"),
    stringsAsFactors = FALSE
  )

  expect_error(
    hsquared:::hs_build_model_spec(
      cbind(y1, y1 + y2) ~ animal(1 | id, pedigree = ped),
      data = dm,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "require bare trait columns inside `cbind()`",
    fixed = TRUE
  )

  expect_error(
    hsquared:::hs_build_model_spec(
      cbind(log(y1), y2) ~ animal(1 | id, pedigree = ped),
      data = dm,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "`log(y1)`",
    fixed = TRUE
  )
})

test_that("bare-column cbind() multivariate response still works", {
  ped <- hs_pr_pedigree()
  dm <- data.frame(
    y1 = c(1, 2, 3, 4),
    y2 = c(2, 3, 4, 5),
    id = c("a", "b", "c", "d"),
    stringsAsFactors = FALSE
  )

  spec <- hsquared:::hs_build_model_spec(
    cbind(y1, y2) ~ animal(1 | id, pedigree = ped),
    data = dm,
    family = stats::gaussian(),
    REML = TRUE
  )

  expect_true(spec$response$multivariate)
  expect_equal(spec$response$trait_names, c("y1", "y2"))
  expect_equal(dim(spec$response$values), c(4L, 2L))
})

# ---- #5 dot (.) ----

test_that("bare dot (.) RHS shorthand is rejected by name", {
  ped <- hs_pr_pedigree()
  dat <- hs_pr_data()

  expect_error(
    hsquared:::hs_build_model_spec(
      y ~ . + animal(1 | id, pedigree = ped),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "`.` all-columns shorthand is not supported",
    fixed = TRUE
  )
})

# ---- #6 iterative topological pedigree sort ----

test_that("deep linear pedigree (N = 5000) sorts without stack overflow", {
  n <- 5000L
  ids <- as.character(seq_len(n))
  sire <- c(NA_character_, utils::head(ids, -1L))
  dam <- rep(NA_character_, n)

  result <- hsquared:::hs_topological_pedigree(ids, sire, dam)

  # Already in parent-before-offspring order, so the sort is a no-op reorder.
  expect_equal(result$data$id, ids)
  expect_equal(nrow(result$data), n)
  # Each individual's sire sits one position earlier in the sorted order.
  expect_equal(result$parent_index$sire, c(0L, seq_len(n - 1L)))
  expect_equal(result$parent_index$dam, rep(0L, n))
})

test_that("topological sort places parents strictly before offspring", {
  ids <- c("c", "a", "b")
  sire <- c("a", NA, NA)
  dam <- c("b", NA, NA)

  result <- hsquared:::hs_topological_pedigree(ids, sire, dam)

  expect_equal(result$data$id, c("a", "b", "c"))
  expect_equal(result$parent_index$sire, c(0L, 0L, 1L))
  expect_equal(result$parent_index$dam, c(0L, 0L, 2L))
})

test_that("a true parent-offspring cycle is still detected", {
  expect_error(
    hsquared:::hs_topological_pedigree(
      c("a", "b"),
      c("b", "a"),
      c(NA_character_, NA_character_)
    ),
    "parent-offspring cycle involving ID `a`",
    fixed = TRUE
  )
})
