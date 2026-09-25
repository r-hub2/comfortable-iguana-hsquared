# Opt-in, experimental single-step model: a single primary effect whose
# relationship is a user-supplied single-step relationship inverse (Hinv).
# Surfaces fit_ai_reml on an Hinv-based animal_model_spec. REML only. This
# reuses the supplied-relationship-inverse primary path (cf. genomic).

hs_test_hinv <- function(ids) {
  n <- length(ids)
  h <- diag(n)
  for (i in seq_len(n - 1L)) {
    h[i, i + 1L] <- h[i + 1L, i] <- 0.15
  }
  hinv <- solve(h)
  dimnames(hinv) <- list(ids, ids)
  hinv
}

test_that("the parser accepts single_step(1 | id, Hinv = Hinv) as a primary", {
  ids <- paste0("a", 1:4)
  Hinv <- hs_test_hinv(ids)
  dat <- data.frame(y = c(1, 2, 3, 4), id = ids)

  spec <- hsquared:::hs_build_model_spec(
    y ~ single_step(1 | id, Hinv = Hinv),
    data = dat,
    family = stats::gaussian(),
    REML = TRUE
  )

  expect_named(spec$random, "single_step")
  expect_equal(spec$random$single_step$type, "single_step")
  expect_equal(spec$random$single_step$relationship, "single_step")
  expect_equal(spec$random$single_step$ids, ids)
  expect_match(spec$bridge$target, "Hinv", fixed = TRUE)
})

test_that("single_step() requires an Hinv argument", {
  dat <- data.frame(y = c(1, 2), id = paste0("a", 1:2))
  expect_error(
    hsquared:::hs_build_model_spec(
      y ~ single_step(1 | id),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "requires a `Hinv`",
    fixed = TRUE
  )
})

test_that("single_step() ids must be in the Hinv dimnames", {
  ids <- paste0("a", 1:3)
  Hinv <- hs_test_hinv(ids)
  dat <- data.frame(y = c(1, 2), id = c("a1", "ghost"))
  expect_error(
    hsquared:::hs_build_model_spec(
      y ~ single_step(1 | id, Hinv = Hinv),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "not in the `Hinv`",
    fixed = TRUE
  )
})

test_that("single_step is a valid opt-in julia target", {
  expect_equal(
    hsquared:::hs_validate_julia_target("single_step"),
    "single_step"
  )
})

test_that("the default engine = \"fit\" rejects a single_step() formula", {
  ids <- paste0("a", 1:3)
  Hinv <- hs_test_hinv(ids)
  dat <- data.frame(y = c(1, 2, 3), id = ids)
  expect_error(
    hsquared(
      y ~ single_step(1 | id, Hinv = Hinv),
      data = dat,
      family = stats::gaussian()
    ),
    "experimental and opt-in",
    fixed = TRUE
  )
})

test_that("the default single-step error names both opt-in routes", {
  ids <- paste0("a", 1:3)
  Hinv <- hs_test_hinv(ids)
  dat <- data.frame(y = c(1, 2, 3), id = ids)

  cnd <- expect_error(
    hsquared(
      y ~ single_step(1 | id, Hinv = Hinv),
      data = dat,
      family = stats::gaussian()
    ),
    class = "hsquared_unsupported_syntax"
  )
  expect_match(
    conditionMessage(cnd),
    "`single_step()` is not on the default `engine = \"fit\"` path",
    fixed = TRUE
  )
  expect_match(
    conditionMessage(cnd),
    "target = \"single_step\"",
    fixed = TRUE
  )
  expect_match(
    conditionMessage(cnd),
    "target = \"single_step_construct\"",
    fixed = TRUE
  )
  expect_match(
    conditionMessage(cnd),
    "Neither route is a covered default R formula path.",
    fixed = TRUE
  )
  expect_match(
    conditionMessage(cnd),
    "Ordinary/default single-step stays held",
    fixed = TRUE
  )
  expect_match(
    conditionMessage(cnd),
    "V2-SSHINV",
    fixed = TRUE
  )
  expect_match(
    conditionMessage(cnd),
    "public_covered_count stays 7",
    fixed = TRUE
  )
})

test_that("model_spec() errors clearly on a single_step formula", {
  ids <- paste0("a", 1:2)
  Hinv <- hs_test_hinv(ids)
  dat <- data.frame(y = c(1, 2), id = ids)
  expect_error(
    model_spec(y ~ single_step(1 | id, Hinv = Hinv), data = dat),
    "previews the pedigree animal-model grammar only",
    fixed = TRUE
  )
})

test_that("hsquared fits the opt-in single-step model", {
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for live single-step."
  )

  set.seed(9)
  na <- 8
  ids <- paste0("a", seq_len(na))
  Hinv <- hs_test_hinv(ids)

  n <- 24
  rec <- rep(ids, length.out = n)
  dat <- data.frame(
    y = 3 + stats::rnorm(n, 0, 1),
    id = rec
  )

  fit <- hsquared(
    y ~ single_step(1 | id, Hinv = Hinv),
    data = dat,
    family = stats::gaussian(),
    control = hs_control(
      engine = "julia",
      engine_control = list(target = "single_step")
    )
  )

  expect_s3_class(fit, "hsquared_fit")
  expect_equal(fit$spec$target, "single_step")
  vc <- variance_components(fit)
  expect_equal(vc$component, c("single_step", "residual"))
  expect_true(all(is.finite(vc$estimate)) && all(vc$estimate > 0))
  expect_equal(nrow(breeding_values(fit)), na)
  expect_equal(
    fit_diagnostics(fit)$value[
      fit_diagnostics(fit)$metric == "variance_components_source"
    ],
    "estimated_single_step_ai_reml"
  )
})
