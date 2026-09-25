test_that("R REML reference does not throw a raw chol error near the h2 -> 1 boundary", {
  # Regression: when residual variance is driven to ~0 (h2 -> 1), an
  # overshooting Nelder-Mead step can make V numerically singular / non-PD, and
  # the reference Gaussian log-likelihood used to propagate a raw chol() error
  # ("the leading minor of order 2 is not positive") straight out of optim().
  # The optimizer must instead retreat from the bad region or report a non-zero
  # convergence code; it must never abort.
  set.seed(11)
  a <- 6L
  n <- 5L
  N <- a * n
  u <- rnorm(a, 0, sqrt(50))
  id <- rep(seq_len(a), each = n)
  y <- u[id] + rnorm(N, 0, sqrt(1e-6))
  X <- matrix(1, N, 1L)
  Z <- matrix(0, N, a)
  Z[cbind(seq_len(N), id)] <- 1

  ref <- NULL
  expect_no_error(
    ref <- hsquared:::hs_reml_estimate_reference(
      y,
      X,
      Z,
      diag(a),
      method = "REML"
    )
  )

  # The result is either a finite admissible estimate (optimizer retreated and
  # converged) or it carries a non-zero convergence code so callers can detect
  # failure — but never a silently "converged" non-finite estimate.
  admissible_finite <- all(is.finite(ref$estimate)) &&
    all(ref$estimate > 0) &&
    is.finite(ref$loglik)
  expect_true(admissible_finite || !identical(ref$convergence, 0L))
})

test_that("R REML reference still converges cleanly on well-conditioned data", {
  # Control: a non-degenerate replicated design must still report
  # convergence == 0 with a finite, positive estimate (unchanged behaviour).
  fixture <- hsquared:::hs_mrode_supplied_variance_validation_fixture()
  spec <- hsquared:::hs_build_model_spec(
    fixture$formula,
    data = fixture$data,
    family = stats::gaussian(),
    REML = TRUE
  )
  payload <- hsquared:::hs_build_bridge_payload(spec)
  ref <- hsquared:::hs_reml_estimate_reference(
    payload$y,
    payload$X,
    as.matrix(payload$Z),
    fixture$expected$Ainv,
    method = "REML"
  )

  expect_equal(ref$convergence, 0L)
  expect_true(all(is.finite(ref$estimate)) && all(ref$estimate > 0))
  expect_true(is.finite(ref$loglik))
})
