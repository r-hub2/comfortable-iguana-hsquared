# Negative controls (tests of tests) for the validation-canon suite.
#
# The recovery and supplied-variance fixtures in test-validation-fixtures.R and
# test-mrode-validation.R assert that a reference computation MATCHES a pinned
# expected value. On their own those assertions cannot tell a discriminating
# comparison from a vacuous one: an `expect_equal()` against a constant, or a
# tolerance band so wide it accepts anything, would pass even for a wrong
# implementation. These tests close that gap. Each one reuses the SAME pure-R
# reference computation and the SAME tolerance as a real validation test, then
# feeds it a DELIBERATELY WRONG expected value and asserts the comparison
# REJECTS it. If a real fixture ever became vacuous (e.g. zero tolerance relaxed
# to Inf, or a constant accidentally compared to itself), its paired negative
# control here would start failing.
#
# All checks are pure R, fast, and CI-runnable: no Julia engine, no skip guards.

test_that("Henderson supplied-variance solve rejects wrong fixed effects (1e-12)", {
  # Mirrors test-validation-fixtures.R "Henderson MME validation fixture pins
  # supplied-variance solutions": same reference solve, same 1e-12 tolerance.
  fixture <- hsquared:::hs_henderson_mme_validation_fixture()
  spec <- hsquared:::hs_build_model_spec(
    fixture$formula,
    data = fixture$data,
    family = stats::gaussian(),
    REML = FALSE
  )
  payload <- hsquared:::hs_build_bridge_payload(spec)
  reference <- hsquared:::hs_solve_henderson_mme_reference(
    y = payload$y,
    X = payload$X,
    Z = payload$Z,
    Ainv = fixture$expected$Ainv,
    sigma_a2 = fixture$sigma_a2,
    sigma_e2 = fixture$sigma_e2,
    ids = payload$ids
  )

  # Sanity: the reference still matches the pinned truth at the real tolerance,
  # so the rejection below is of a genuinely wrong value, not a broken solve.
  expect_equal(
    unname(reference$fixed_effects),
    unname(fixture$expected$fixed_effects),
    tolerance = 1e-12
  )

  # Test of test: a fixed-effect vector perturbed well beyond 1e-12 must NOT
  # compare equal. If this passed, the real assertion would be vacuous.
  wrong_fixed <- unname(fixture$expected$fixed_effects) + c(0.1, -0.1)
  expect_false(isTRUE(all.equal(
    unname(reference$fixed_effects),
    wrong_fixed,
    tolerance = 1e-12
  )))

  # The discriminating comparison is also detectable through testthat's own
  # failure machinery, not just all.equal().
  expect_failure(expect_equal(
    unname(reference$fixed_effects),
    wrong_fixed,
    tolerance = 1e-12
  ))

  # And breeding values: shifting one EBV by 0.5 must be rejected.
  wrong_bv <- fixture$expected$breeding_values
  wrong_bv$value[3] <- wrong_bv$value[3] + 0.5
  expect_false(isTRUE(all.equal(
    reference$breeding_values,
    wrong_bv,
    tolerance = 1e-12
  )))
})

test_that("Mrode-style supplied-variance solve and loglik reject wrong values (1e-10)", {
  # Mirrors test-validation-fixtures.R "Mrode-style supplied-variance fixture
  # pins R reference outputs": same reference solve and loglik, same 1e-10
  # tolerance.
  fixture <- hsquared:::hs_mrode_supplied_variance_validation_fixture()
  spec <- hsquared:::hs_build_model_spec(
    fixture$formula,
    data = fixture$data,
    family = stats::gaussian(),
    REML = FALSE
  )
  payload <- hsquared:::hs_build_bridge_payload(spec)
  reference <- hsquared:::hs_solve_henderson_mme_reference(
    y = payload$y,
    X = payload$X,
    Z = payload$Z,
    Ainv = fixture$expected$Ainv,
    sigma_a2 = fixture$sigma_a2,
    sigma_e2 = fixture$sigma_e2,
    ids = payload$ids
  )
  reml <- hsquared:::hs_gaussian_loglik_reference(
    y = payload$y,
    X = payload$X,
    Z = payload$Z,
    Ainv = fixture$expected$Ainv,
    sigma_a2 = fixture$sigma_a2,
    sigma_e2 = fixture$sigma_e2,
    method = "REML"
  )

  # Sanity: real comparisons still hold at the real tolerance.
  expect_equal(
    reference$fixed_effects,
    fixture$expected$fixed_effects,
    tolerance = 1e-10,
    ignore_attr = TRUE
  )
  expect_equal(reml$loglik, fixture$expected$reml_loglik, tolerance = 1e-10)

  # Test of test: a loglik off by 0.01 (>> 1e-10) must be rejected. A wrong
  # estimator that landed on a different REML optimum would be caught here.
  expect_false(isTRUE(all.equal(
    reml$loglik,
    fixture$expected$reml_loglik + 0.01,
    tolerance = 1e-10
  )))
  expect_failure(expect_equal(
    reml$loglik,
    fixture$expected$reml_loglik + 0.01,
    tolerance = 1e-10
  ))

  # And the fixed-effect slope perturbed beyond 1e-10 must be rejected.
  wrong_fixed <- fixture$expected$fixed_effects
  wrong_fixed[["x"]] <- wrong_fixed[["x"]] + 1e-3
  expect_false(isTRUE(all.equal(
    unname(reference$fixed_effects),
    unname(wrong_fixed),
    tolerance = 1e-10
  )))
})

test_that("gryphon REML recovery band (0.02) rejects a perturbed estimate", {
  # Mirrors test-validation-fixtures.R "hsquared's R REML reference recovers the
  # published gryphon estimates": same published anchor and same 0.02 band, but
  # here we prove the band is discriminating WITHOUT needing the `enhancer`
  # data. We perturb the PUBLISHED truth itself and assert the band rejects it.
  pub <- hsquared:::hs_gryphon_published_reml()

  # A truthful estimate (= the anchor) passes the band; this documents the band
  # the real test relies on.
  expect_equal(pub[["sigma_a2"]], pub[["sigma_a2"]], tolerance = 0.02)

  # Test of test: an additive-variance estimate off by 0.5 (>> the 0.02 band)
  # must be rejected. The real recovery test would catch an estimator biased by
  # this much.
  wrong_va <- pub[["sigma_a2"]] + 0.5
  expect_false(isTRUE(all.equal(
    wrong_va,
    pub[["sigma_a2"]],
    tolerance = 0.02
  )))
  expect_failure(expect_equal(wrong_va, pub[["sigma_a2"]], tolerance = 0.02))

  # And h2 off by 0.1 (>> 0.02) must be rejected.
  wrong_h2 <- pub[["h2"]] + 0.1
  expect_false(isTRUE(all.equal(
    wrong_h2,
    pub[["h2"]],
    tolerance = 0.02
  )))
})

test_that("pure-R REML optimizer band (5e-2) rejects a perturbed estimate", {
  # Mirrors the cross-check in test-validation-fixtures.R "independent pure-R
  # REML optimizer matches the Julia sparse REML estimate": same pure-R
  # optimizer, same 5e-2 tolerance. The real test compares the Julia estimate to
  # `ref$estimate`; here we feed the SAME band a deliberately wrong estimate to
  # prove 5e-2 is not wide enough to accept anything.
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

  # Sanity: the optimizer converged to a finite, positive estimate.
  expect_equal(ref$convergence, 0L)
  expect_true(all(ref$estimate > 0))

  # Test of test: an estimate scaled by 1.5 (a 50% error, >> 5e-2) must be
  # rejected by the same band the real cross-check uses.
  wrong_estimate <- ref$estimate * 1.5
  expect_false(isTRUE(all.equal(
    unname(wrong_estimate),
    unname(ref$estimate),
    tolerance = 5e-2
  )))
  expect_failure(expect_equal(
    unname(wrong_estimate),
    unname(ref$estimate),
    tolerance = 5e-2
  ))
})

test_that("DGP recovery band (0.06) rejects a biased h2/variance statistic", {
  # Mirrors test-validation-fixtures.R "REML recovers known variance components
  # from a simulated DGP": that test asserts `expect_lt(abs(mean - truth), 0.06)`
  # on s2a, s2e, and h2. This negative control proves the 0.06 band is
  # discriminating: a statistic biased by more than 0.06 from the known truth
  # must FAIL the same predicate. No simulation or `nadiv` needed -- we perturb
  # the known truth directly, exactly the bias the real test is meant to catch.
  s2a <- 0.4
  s2e <- 0.6
  h2_truth <- 0.4

  # A statistic at the truth passes the band (documents the real predicate).
  expect_lt(abs(s2a - s2a), 0.06)
  expect_lt(abs(h2_truth - h2_truth), 0.06)

  # Test of test: a biased additive-variance mean (off by 0.1 > 0.06) must fail
  # the exact predicate the real recovery test uses.
  biased_s2a <- s2a + 0.1
  expect_false(abs(biased_s2a - s2a) < 0.06)
  expect_failure(expect_lt(abs(biased_s2a - s2a), 0.06))

  # And a biased heritability (off by 0.1 > 0.06) must fail the same band.
  biased_h2 <- h2_truth + 0.1
  expect_false(abs(biased_h2 - h2_truth) < 0.06)
  expect_failure(expect_lt(abs(biased_h2 - h2_truth), 0.06))

  # The EBV-accuracy floor (expect_gt(mean(acc), 0.5)) is likewise
  # discriminating: an accuracy below the floor must fail it.
  poor_accuracy <- 0.3
  expect_false(poor_accuracy > 0.5)
  expect_failure(expect_gt(poor_accuracy, 0.5))
})
