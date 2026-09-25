# Experimental multivariate covariance standard errors (engine row V4-MV-REML,
# partial; unstructured-only). R-side extractor tested with fixtures; the live
# engine path is exercised opportunistically by the multivariate bridge.

test_that("covariance_standard_errors() returns the SE list when present", {
  se <- list(
    genetic_covariance = matrix(c(0.1, 0.05, 0.05, 0.12), 2, 2),
    residual_covariance = matrix(c(0.2, 0.06, 0.06, 0.15), 2, 2),
    genetic_correlation = matrix(c(0, 0.08, 0.08, 0), 2, 2),
    residual_correlation = matrix(c(0, 0.09, 0.09, 0), 2, 2),
    heritability = c(t1 = 0.07, t2 = 0.06)
  )
  fit <- hsquared:::hs_new_fit(
    call = quote(hsquared(
      cbind(t1, t2) ~ animal(1 | id, pedigree = ped),
      data = dat
    )),
    spec = list(method = "REML", family = list(family = "gaussian")),
    payload = list(Y = matrix(0, 4, 2)),
    result = list(covariance_standard_errors = se)
  )
  expect_s3_class(fit, "hsquared_fit")
  expect_equal(covariance_standard_errors(fit), se)
})

test_that("covariance_standard_errors() errors clearly without the field or object", {
  expect_error(
    covariance_standard_errors(seq_len(10)),
    "requires an `hsquared_fit` object"
  )
  fit_no_se <- hsquared:::hs_new_fit(
    call = quote(hsquared(
      cbind(t1, t2) ~ animal(1 | id, pedigree = ped),
      data = dat
    )),
    spec = list(method = "REML", family = list(family = "gaussian")),
    payload = list(Y = matrix(0, 4, 2)),
    result = list(genetic_covariance = diag(2))
  )
  expect_error(
    covariance_standard_errors(fit_no_se),
    "experimental multivariate covariance standard errors"
  )
})
