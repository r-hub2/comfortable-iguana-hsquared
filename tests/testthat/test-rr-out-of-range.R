# hsquared issue #213: the R-side normalized-Legendre basis silently clamped
# a standardized covariate to [-1, 1] instead of erroring like the Julia twin
# (`HSquared.legendre_basis`, tolerance |t| <= 1 + 1e-10). These tests pin the
# error-not-clamp behaviour on every layer of the blast radius: the basis
# helper itself, the shared `at`-standardization helper, and a public
# extractor. No Julia is needed -- the fixture mirrors
# test-random-regression.R's hand-built `hsquared_fit`.

test_that("hs_legendre_basis errors beyond tolerance, not within it", {
  # Clearly out of range: errors.
  expect_error(
    hsquared:::hs_legendre_basis(1 + 1e-6, 2L),
    class = "hsquared_error"
  )
  expect_error(
    hsquared:::hs_legendre_basis(-1 - 1e-6, 2L),
    class = "hsquared_error"
  )

  # Floating-point overshoot within the Julia tolerance (1e-10): no error.
  expect_no_error(hsquared:::hs_legendre_basis(1 + 1e-11, 2L))
  expect_no_error(hsquared:::hs_legendre_basis(-1 - 1e-11, 2L))

  # In-range values are unaffected.
  expect_equal(
    hsquared:::hs_legendre_basis(1, 2L),
    hsquared:::hs_legendre_basis(1 + 1e-11, 2L),
    tolerance = 1e-8
  )
})

# Hand-built `rr` metadata, mirroring the shape `hs_rr_eval_points()` reads
# off `object$result$random_regression` (see test-random-regression.R).
hs_rr_out_of_range_fit <- function() {
  payload <- structure(
    list(
      y = rep(c(1, 2, 3), 4),
      X = matrix(1, nrow = 12L, ncol = 1L),
      ids = c("a", "b", "c", "d"),
      family = "gaussian",
      random_regression = list(
        covariate = "age",
        order = 2L,
        lower = 1,
        upper = 5
      ),
      metadata = list(
        fixed_colnames = "(Intercept)",
        random_regression = list(
          covariate = "age",
          order = 2L,
          lower = 1,
          upper = 5
        )
      )
    ),
    class = c("hs_bridge_payload", "list")
  )
  raw <- list(
    K_g = matrix(c(0.8, 0.1, 0.1, 0.3), 2L),
    sigma_e2 = 0.5,
    beta = 2.1,
    coef_ids = payload$ids,
    coef_values = matrix(seq(0.1, 0.8, length.out = 8L), nrow = 4L),
    loglik = -10.2,
    converged = TRUE,
    iterations = 33L,
    ncoef = 2L
  )
  result <- hsquared:::hs_normalize_random_regression_result(raw, payload)
  hsquared:::hs_new_fit(
    spec = list(
      method = "REML",
      family = list(family = "gaussian"),
      target = "random_regression"
    ),
    payload = payload,
    result = result
  )
}

test_that("hs_rr_eval_points rejects `at` outside the fitted range", {
  fit <- hs_rr_out_of_range_fit()
  rr <- fit$result$random_regression
  expect_equal(rr$lower, 1)
  expect_equal(rr$upper, 5)

  err <- tryCatch(
    hsquared:::hs_rr_eval_points(fit, at = c(-40, 1e6)),
    error = function(e) e
  )
  expect_s3_class(err, "hsquared_error")
  expect_match(err$message, "1", fixed = TRUE)
  expect_match(err$message, "5", fixed = TRUE)

  # In-range `at` still works.
  pts <- hsquared:::hs_rr_eval_points(fit, at = c(1, 3, 5))
  expect_equal(pts$at, c(1, 3, 5))
})

test_that("rr_heritability() errors on an out-of-range `at`, same class", {
  fit <- hs_rr_out_of_range_fit()

  # In range: works.
  expect_no_error(rr_heritability(fit, at = 3))

  # Out of range: errors with the same structured class as the extractor
  # helper, naming the fitted range.
  err <- tryCatch(
    rr_heritability(fit, at = 1e6),
    error = function(e) e
  )
  expect_s3_class(err, "hsquared_error")
  expect_match(err$message, "1", fixed = TRUE)
  expect_match(err$message, "5", fixed = TRUE)
})
