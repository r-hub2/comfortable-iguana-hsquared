# Experimental heritability_interval extractor (engine row V1-HERIT-CI, partial).
# These tests exercise the R-side extractor + normalizer with fixtures; the live
# engine path is exercised opportunistically by the bridge when Julia is present.

test_that("heritability_interval() returns the interval when present", {
  result <- list(
    heritability_interval = data.frame(
      estimate = 0.42,
      lower = 0.21,
      upper = 0.66,
      level = 0.95,
      se = 0.11,
      method = "delta",
      stringsAsFactors = FALSE
    )
  )
  fit <- hsquared:::hs_new_fit(
    call = quote(hsquared(y ~ animal(1 | id, pedigree = ped), data = dat)),
    spec = list(method = "REML", family = list(family = "gaussian")),
    payload = list(y = seq_len(10)),
    result = result
  )

  expect_s3_class(fit, "hsquared_fit")
  expect_equal(heritability_interval(fit), result$heritability_interval)
})

test_that("heritability_interval() errors clearly without the field or object", {
  expect_error(
    heritability_interval(seq_len(10)),
    "requires an `hsquared_fit` object"
  )

  fit_no_ci <- hsquared:::hs_new_fit(
    call = quote(hsquared(y ~ animal(1 | id, pedigree = ped), data = dat)),
    spec = list(method = "REML", family = list(family = "gaussian")),
    payload = list(y = seq_len(10)),
    result = list(heritability = data.frame(term = "animal", estimate = 0.4))
  )
  expect_error(
    heritability_interval(fit_no_ci),
    "experimental heritability confidence interval"
  )
})

test_that("hs_normalize_heritability_interval() handles delta and profile shapes", {
  delta <- hsquared:::hs_normalize_heritability_interval(
    list(
      heritability = 0.42,
      lower = 0.21,
      upper = 0.66,
      level = 0.95,
      se = 0.11,
      method = "delta"
    )
  )
  expect_equal(nrow(delta), 1L)
  expect_equal(delta$estimate, 0.42)
  expect_true(delta$lower <= delta$estimate && delta$estimate <= delta$upper)
  expect_equal(delta$se, 0.11)
  expect_equal(delta$method, "delta")

  # The profile method returns no `se`; it must normalize to NA, not error.
  profile <- hsquared:::hs_normalize_heritability_interval(
    list(
      heritability = 0.40,
      lower = 0.20,
      upper = 0.70,
      level = 0.95,
      method = "profile"
    )
  )
  expect_true(is.na(profile$se))
  expect_equal(profile$method, "profile")
})
