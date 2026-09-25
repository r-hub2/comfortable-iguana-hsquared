# Experimental covariance-structure LRT (diagonal vs unstructured genetic
# covariance). Computed R-side from two multivariate fits' stored `loglik` +
# `n_genetic_params` (the twin's #61 contract), so it is fully fixture-testable
# without a live engine. Engine row V4-MV-REML (partial).

make_mv_fit <- function(loglik, n_genetic_params, genetic_structure) {
  hsquared:::hs_new_fit(
    call = quote(hsquared(
      cbind(t1, t2) ~ animal(1 | id, pedigree = ped),
      data = dat
    )),
    spec = list(method = "REML", family = list(family = "gaussian")),
    payload = list(Y = matrix(0, 4, 2)),
    result = list(
      loglik = loglik,
      n_genetic_params = n_genetic_params,
      genetic_structure = genetic_structure
    )
  )
}

test_that("covariance_structure_lrt computes the diagonal-vs-unstructured test", {
  diag_fit <- make_mv_fit(-110, 2L, "diagonal")
  full_fit <- make_mv_fit(-108, 3L, "unstructured")
  lrt <- covariance_structure_lrt(diag_fit, full_fit)

  expect_equal(lrt$df, 1L) # t(t-1)/2 = 1 off-diagonal genetic covariance for t=2
  expect_equal(lrt$statistic, 2 * (-108 - -110)) # 2*(ll_full - ll_diag) = 4
  expect_false(lrt$boundary) # diagonal-in-unstructured is an interior null
  expect_equal(lrt$pvalue, stats::pchisq(4, 1, lower.tail = FALSE))
  expect_equal(lrt$constrained, "diagonal")
  expect_equal(lrt$full, "unstructured")
})

test_that("covariance_structure_lrt guards order, object class, and missing fields", {
  diag_fit <- make_mv_fit(-110, 2L, "diagonal")
  full_fit <- make_mv_fit(-108, 3L, "unstructured")

  # Wrong order (df <= 0).
  expect_error(
    covariance_structure_lrt(full_fit, diag_fit),
    "more genetic covariance parameters"
  )
  # Not hsquared_fit objects.
  expect_error(
    covariance_structure_lrt(seq_len(3), full_fit),
    "must both be"
  )
  # Missing loglik (e.g. a non-converged fit).
  no_ll <- hsquared:::hs_new_fit(
    call = quote(hsquared(
      cbind(t1, t2) ~ animal(1 | id, pedigree = ped),
      data = dat
    )),
    spec = list(method = "REML", family = list(family = "gaussian")),
    payload = list(Y = matrix(0, 4, 2)),
    result = list(n_genetic_params = 2L, genetic_structure = "diagonal")
  )
  expect_error(covariance_structure_lrt(no_ll, full_fit), "loglik")
})

hs_lrt_fixture_meta <- function(dir, key) {
  meta <- utils::read.csv(
    testthat::test_path("fixtures", dir, "expected_metadata.csv"),
    stringsAsFactors = FALSE
  )
  stats::setNames(meta$value, meta$key)[[key]]
}

test_that("covariance_structure_lrt runs end-to-end on the shared fixtures", {
  # The diagonal and unstructured targets are fitted on IDENTICAL inputs (the
  # `structured_covariance_parity` and `phase4_multitrait_parity` fixtures share
  # the same pedigree + phenotypes), so the two REML log-likelihoods form a
  # valid nested diagonal-vs-unstructured structure test.
  ll_diag <- as.numeric(
    hs_lrt_fixture_meta("structured_covariance_parity", "loglik")
  )
  ll_full <- as.numeric(
    hs_lrt_fixture_meta("phase4_multitrait_parity", "loglik")
  )
  # The diagonal genetic-parameter count is read from the fixture (it records
  # n_genetic_params = t); the unstructured count is the derived t(t+1)/2 = 3
  # for t = 2 (the phase4 fixture predates the n_genetic_params field).
  np_diag <- as.integer(
    hs_lrt_fixture_meta("structured_covariance_parity", "n_genetic_params")
  )
  expect_equal(np_diag, 2L)

  diag_fit <- make_mv_fit(ll_diag, np_diag, "diagonal")
  full_fit <- make_mv_fit(ll_full, 3L, "unstructured")
  lrt <- covariance_structure_lrt(diag_fit, full_fit)

  expect_equal(lrt$df, 1L) # t(t-1)/2 off-diagonal genetic covariances, t = 2
  expect_false(lrt$boundary) # interior null
  expect_equal(lrt$statistic, 2 * (ll_full - ll_diag))
  expect_gt(lrt$statistic, 0) # the unstructured fit cannot do worse
  expect_equal(
    lrt$pvalue,
    stats::pchisq(2 * (ll_full - ll_diag), df = 1, lower.tail = FALSE)
  )
  expect_equal(lrt$constrained, "diagonal")
  expect_equal(lrt$full, "unstructured")
})

test_that("covariance_structure_lrt flags a non-interior null and clamps a negative statistic", {
  # Any pairing other than diagonal-in-unstructured is boundary-conservative
  # (the naive chi-square is not valid at a variance boundary).
  diag_fit <- make_mv_fit(-110, 2L, "diagonal")
  lowrank_fit <- make_mv_fit(-108, 3L, "lowrank")
  expect_true(covariance_structure_lrt(diag_fit, lowrank_fit)$boundary)

  # A marginally-negative 2*Δloglik (optimizer noise) is clamped so the p-value
  # stays a valid probability rather than exceeding 1.
  worse_full <- make_mv_fit(-110.0001, 3L, "unstructured")
  noisy <- covariance_structure_lrt(diag_fit, worse_full)
  expect_lt(noisy$statistic, 0)
  expect_equal(noisy$pvalue, stats::pchisq(0, df = 1, lower.tail = FALSE))
})
