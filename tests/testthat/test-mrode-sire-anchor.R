# Published external-canon anchor #2: Mrode (2014), 3rd ed., Example 3.2 (p.48) —
# the SIRE model on the same WWG data as Example 3.1, with related sires.
#
# This extends the published-EBV canon (test-mrode-published-anchor.R covers the
# animal model, 3.1) to a SECOND model class. Inputs + published solutions are
# confirmed against the masuday BLUPF90 tutorial (p.48) and the austin-putz
# chapter-3 reproduction, and independently re-solved (~1e-7). See
# hs_mrode_example_3_2_sire_fixture() for full provenance.
#
# The sire model is not a parsed hsquared path (sire models are framed as
# planned), so X / Z are built here and the general reference Henderson solver is
# pinned to the published digits. Pure R, CI-runnable: no Julia, no skip guards.

test_that("the Henderson MME solver reproduces the published Mrode Example 3.2 sire solutions", {
  fixture <- hsquared:::hs_mrode_example_3_2_sire_fixture()
  d <- fixture$data
  sire_ids <- fixture$sire_ids

  # No-intercept sex design (two identifiable means) + sire incidence.
  X <- stats::model.matrix(~ sex - 1, data = d)
  Z <- matrix(0, nrow(d), length(sire_ids))
  Z[cbind(seq_len(nrow(d)), match(d$sire, sire_ids))] <- 1

  expect_equal(dim(X), c(5L, 2L))
  expect_equal(dim(Z), c(5L, 3L))

  reference <- hsquared:::hs_solve_henderson_mme_reference(
    y = d$WWG,
    X = X,
    Z = Z,
    Ainv = fixture$expected$Ainv,
    sigma_a2 = fixture$sigma_s2,
    sigma_e2 = fixture$sigma_e2,
    ids = sire_ids
  )

  # Load-bearing: the solver's sire solutions equal the PUBLISHED p.48 digits.
  published_sire <- fixture$expected$sire_solutions[
    reference$breeding_values$id
  ]
  expect_equal(
    reference$breeding_values$value,
    unname(published_sire),
    tolerance = 1e-6
  )
  expect_equal(reference$breeding_values$id, sire_ids)

  # The sex block is full rank under the no-intercept design, so both published
  # sex means are directly reproducible.
  fe <- reference$fixed_effects
  expect_equal(
    unname(fe["sexmale"]),
    fixture$expected$sex_male,
    tolerance = 1e-6
  )
  expect_equal(
    unname(fe["sexfemale"]),
    fixture$expected$sex_female,
    tolerance = 1e-6
  )
  expect_equal(
    unname(fe["sexmale"] - fe["sexfemale"]),
    fixture$expected$sex_contrast_male_minus_female,
    tolerance = 1e-6
  )
})

test_that("the published Mrode 3.2 sire anchor rejects perturbed solutions (test of test)", {
  fixture <- hsquared:::hs_mrode_example_3_2_sire_fixture()
  d <- fixture$data
  sire_ids <- fixture$sire_ids
  X <- stats::model.matrix(~ sex - 1, data = d)
  Z <- matrix(0, nrow(d), length(sire_ids))
  Z[cbind(seq_len(nrow(d)), match(d$sire, sire_ids))] <- 1

  reference <- hsquared:::hs_solve_henderson_mme_reference(
    y = d$WWG,
    X = X,
    Z = Z,
    Ainv = fixture$expected$Ainv,
    sigma_a2 = fixture$sigma_s2,
    sigma_e2 = fixture$sigma_e2,
    ids = sire_ids
  )
  published <- unname(
    fixture$expected$sire_solutions[reference$breeding_values$id]
  )
  expect_equal(reference$breeding_values$value, published, tolerance = 1e-6)
  expect_failure(expect_equal(
    reference$breeding_values$value,
    published + 0.1,
    tolerance = 1e-6
  ))
})
