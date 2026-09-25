# G-matrix geometry / evolvability extractors (Hansen & Houle 2008). Computed in
# R from the fitted genetic covariance; correctness is pinned to hand-computed
# values, and a skip-guarded live test verifies parity with the engine's
# evolvability.jl definitions.

hs_make_g_fit <- function(G) {
  hsquared:::hs_new_fit(
    spec = list(
      method = "REML",
      family = list(family = "gaussian"),
      target = "multivariate"
    ),
    payload = list(Y = matrix(0, 4, nrow(G))),
    result = list(genetic_covariance = G)
  )
}

test_that("G-matrix geometry extractors match hand-computed values", {
  G <- matrix(
    c(2, 0.5, 0.5, 1),
    2,
    2,
    dimnames = list(
      c("t1", "t2"),
      c("t1", "t2")
    )
  )
  fit <- hs_make_g_fit(G)

  expect_equal(mean_evolvability(fit), 1.5) # tr(G)/2 = 3/2
  expect_equal(evolvability(fit, c(1, 0)), 2) # G[1,1]
  expect_equal(evolvability(fit, c(0, 1)), 1) # G[2,2]
  expect_equal(evolvability(fit, c(2, 0)), 2) # direction is normalised
  expect_equal(respondability(fit, c(1, 0)), sqrt(2^2 + 0.5^2))
  # variance_along_gradient: raw uses the gradient as given; normalized = evolvability.
  expect_equal(variance_along_gradient(fit, c(1, 0), normalize = FALSE), 2)
  expect_equal(variance_along_gradient(fit, c(2, 0), normalize = FALSE), 8)
  expect_equal(variance_along_gradient(fit, c(2, 0), normalize = TRUE), 2)
  expect_equal(conditional_evolvability(fit, c(1, 0)), 1.75) # det/G22 = 1.75/1
  expect_equal(autonomy(fit, c(1, 0)), 0.875) # 1.75 / 2

  eg <- eigen_G(fit)
  expect_equal(eg$values, c((3 + sqrt(2)) / 2, (3 - sqrt(2)) / 2))
  expect_equal(dim(eg$vectors), c(2L, 2L))
  # eigenvectors are sign-canonicalised: largest-magnitude element positive.
  expect_true(all(vapply(
    seq_len(ncol(eg$vectors)),
    function(j) eg$vectors[which.max(abs(eg$vectors[, j])), j] > 0,
    logical(1L)
  )))

  gm <- g_max(fit)
  expect_equal(gm$eigenvalue, (3 + sqrt(2)) / 2)
  expect_equal(length(gm$eigenvector), 2L)
})

test_that("G-matrix geometry guards inputs, singular G, and non-multivariate fits", {
  fit <- hs_make_g_fit(matrix(c(2, 0.5, 0.5, 1), 2, 2))
  expect_error(evolvability(fit, c(1, 2, 3)), "one entry per trait")
  expect_error(evolvability(fit, c(0, 0)), "nonzero")

  # A singular (reduced-rank) G has no conditional evolvability / autonomy.
  singular <- hs_make_g_fit(matrix(c(1, 1, 1, 1), 2, 2))
  expect_error(conditional_evolvability(singular, c(1, 0)), "positive-definite")
  expect_error(autonomy(singular, c(1, 0)), "positive-definite")
  # ...but the rotation-invariant eigenstructure is still defined.
  expect_equal(eigen_G(singular)$values, c(2, 0))

  # A non-multivariate fit has no genetic covariance matrix.
  univ <- hsquared:::hs_new_fit(
    spec = list(target = "ai_reml"),
    payload = list(y = 1:4),
    result = list()
  )
  expect_error(eigen_G(univ), "multivariate")
  expect_error(mean_evolvability(univ), "multivariate")

  # Asymmetric input is rejected.
  asym <- hs_make_g_fit(matrix(c(2, 0.5, 0.9, 1), 2, 2))
  expect_error(eigen_G(asym), "symmetric")
})

test_that("G-matrix geometry matches the engine (live parity)", {
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for the live evolvability parity check."
  )

  set.seed(7)
  A <- matrix(stats::rnorm(9), 3, 3)
  G <- crossprod(A) + diag(3) # a random symmetric positive-definite 3x3 G
  fit <- hs_make_g_fit(G)
  beta <- c(1, -2, 0.5)

  hsquared:::hs_julia_setup(hsquared:::hs_default_julia_project())
  JuliaCall::julia_assign("hsq_Gtest", G)
  JuliaCall::julia_assign("hsq_beta", beta)

  expect_equal(
    evolvability(fit, beta),
    JuliaCall::julia_eval("HSquared.evolvability(hsq_Gtest, hsq_beta)")
  )
  expect_equal(
    respondability(fit, beta),
    JuliaCall::julia_eval("HSquared.respondability(hsq_Gtest, hsq_beta)")
  )
  expect_equal(
    conditional_evolvability(fit, beta),
    JuliaCall::julia_eval(
      "HSquared.conditional_evolvability(hsq_Gtest, hsq_beta)"
    )
  )
  expect_equal(
    autonomy(fit, beta),
    JuliaCall::julia_eval("HSquared.autonomy(hsq_Gtest, hsq_beta)")
  )
  expect_equal(
    mean_evolvability(fit),
    JuliaCall::julia_eval("HSquared.mean_evolvability(hsq_Gtest)")
  )
  expect_equal(
    eigen_G(fit)$values,
    JuliaCall::julia_eval("HSquared.genetic_pca(hsq_Gtest).values")
  )
  expect_equal(
    g_max(fit)$eigenvalue,
    JuliaCall::julia_eval("HSquared.g_max(hsq_Gtest).eigenvalue")
  )
  expect_equal(
    variance_along_gradient(fit, beta),
    JuliaCall::julia_eval(
      "HSquared.variance_along_gradient(hsq_Gtest, hsq_beta)"
    )
  )
  expect_equal(
    variance_along_gradient(fit, beta, normalize = FALSE),
    JuliaCall::julia_eval(
      "HSquared.variance_along_gradient(hsq_Gtest, hsq_beta; normalize = false)"
    )
  )
  # normalize = TRUE equals evolvability().
  expect_equal(variance_along_gradient(fit, beta), evolvability(fit, beta))
})
