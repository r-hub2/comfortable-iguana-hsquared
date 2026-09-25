test_that("model_spec previews the parsed v0.1 animal contract", {
  ped <- data.frame(
    id = c("a", "b", "c", "d"),
    sire = c(NA, NA, "a", "a"),
    dam = c(NA, NA, "b", "c")
  )
  dat <- data.frame(
    y = c(1, 2, 3),
    sex = c("f", "m", "f"),
    age = c(1, 2, 3),
    id = c("a", "c", "d")
  )

  spec <- model_spec(
    y ~ sex + age + animal(1 | id, pedigree = ped),
    data = dat,
    family = stats::gaussian(),
    REML = TRUE
  )

  expect_s3_class(spec, "hs_model_spec")
  expect_s3_class(spec$payload, "hs_bridge_payload")
  expect_equal(spec$response, "y")
  expect_equal(spec$method, "REML")
  expect_equal(spec$dimensions$observations, 3L)
  expect_equal(spec$dimensions$fixed_columns, 3L)
  expect_equal(spec$dimensions$animal_ids, 4L)
  expect_equal(spec$dimensions$random_design_nonzeros, 3L)
  expect_equal(spec$dimensions$pedigree_founders, 2L)
  expect_equal(spec$fixed$columns, c("(Intercept)", "sexm", "age"))
  expect_equal(spec$animal$group, "id")
  expect_equal(spec$animal$ids, c("a", "b", "c", "d"))
  expect_equal(spec$julia$ainv_status, "build_in_julia")
  expect_match(spec$julia$fit_target, "fit_animal_model", fixed = TRUE)
})

test_that("model_spec summaries are compact and do not fit models", {
  ped <- data.frame(
    id = c("a", "b"),
    sire = c(NA, NA),
    dam = c(NA, NA)
  )
  dat <- data.frame(y = c(1, 2), id = c("a", "b"))

  spec <- model_spec(y ~ animal(1 | id, pedigree = ped), data = dat)
  out <- summary(spec)

  expect_s3_class(out, "summary_hs_model_spec")
  expect_equal(out$model$response, "y")
  expect_equal(out$model$family, "gaussian")
  expect_equal(out$model$method, "REML")
  expect_equal(out$animal$relationship, "pedigree")
  expect_equal(out$julia$stage, c("Ainv", "model specification", "fit"))
  expect_equal(
    out$julia$target[[3L]],
    "HSquared.fit_animal_model(y, X, Z, Ainv; method = :REML)"
  )
  expect_null(spec$result)
})

test_that("model_spec uses the same honest planned-marker errors", {
  dat <- data.frame(y = c(1, 2), id = c("a", "b"))

  expect_error(
    model_spec(y ~ marker_scan(M, map = map), data = dat),
    "`marker_scan()` is planned, not implemented.",
    fixed = TRUE
  )

  # `genomic()` is parsed (opt-in primary effect) but `model_spec()` previews the
  # pedigree animal grammar only; it errors with a clear pointer to the opt-in fit.
  ids <- c("a", "b")
  Ginv <- diag(2)
  dimnames(Ginv) <- list(ids, ids)
  gdat <- data.frame(y = c(1, 2), id = ids)
  expect_error(
    model_spec(y ~ genomic(1 | id, Ginv = Ginv), data = gdat),
    "previews the pedigree animal-model grammar only",
    fixed = TRUE
  )
  expect_error(
    model_spec(
      y ~ animal(1 | id, pedigree = ped),
      data = dat,
      extra = TRUE
    ),
    "`...` is reserved for future `model_spec()` options.",
    fixed = TRUE
  )
})

test_that("model_spec can read phenotypes and pedigree from hs_data", {
  ped <- data.frame(
    id = c("a", "b", "c"),
    sire = c(NA, NA, "a"),
    dam = c(NA, NA, "b")
  )
  dat <- data.frame(
    y = c(1, 2),
    age = c(4, 5),
    id = c("a", "c")
  )
  bundle <- hs_data(
    phenotypes = dat,
    pedigree = ped
  )

  spec <- model_spec(
    y ~ age + animal(1 | id, pedigree = pedigree),
    data = bundle
  )

  expect_equal(spec$dimensions$observations, 2L)
  expect_equal(spec$animal$ids, c("a", "b", "c"))
  expect_equal(spec$animal$observed_id_index, c(1L, 3L))
  expect_equal(spec$payload$pedigree$sire_index, c(0L, 0L, 1L))
  expect_equal(spec$payload$pedigree$dam_index, c(0L, 0L, 2L))
})

test_that("model_spec uses the hs_data pedigree when animal omits pedigree", {
  ped <- data.frame(
    id = c("a", "b", "c"),
    sire = c(NA, NA, "a"),
    dam = c(NA, NA, "b")
  )
  dat <- data.frame(
    y = c(1, 2),
    age = c(4, 5),
    id = c("a", "c")
  )
  bundle <- hs_data(phenotypes = dat, pedigree = ped)

  spec <- model_spec(y ~ age + animal(1 | id), data = bundle)
  out <- summary(spec)

  expect_equal(spec$animal$pedigree_source, "hs_data")
  expect_equal(out$animal$pedigree_source, "hs_data")
  expect_equal(spec$animal$ids, c("a", "b", "c"))
  expect_equal(spec$payload$pedigree$sire_index, c(0L, 0L, 1L))
  expect_equal(spec$payload$pedigree$dam_index, c(0L, 0L, 2L))
})
