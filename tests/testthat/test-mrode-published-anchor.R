# Published external-canon anchor: Mrode (2014), "Linear Models for the
# Prediction of Animal Breeding Values", 3rd ed., Example 3.1 (p.39).
#
# The other pedigree anchors (test-pedigree-mme-anchor.R) prove the package
# solver agrees with an INDEPENDENT hand solve of the SAME MME -- closing the
# self-generated-number circularity. This file closes the remaining gap: it pins
# the package solver against the PUBLISHED textbook EBV digits. A solver that is
# internally self-consistent but systematically wrong (e.g. a mis-scaled lambda,
# a wrong Ainv convention) would pass the hand-solve anchors yet fail here.
#
# The Example 3.1 inputs and published solutions are confirmed against three
# independent citable sources (the masuday BLUPF90 tutorial citing Mrode 2014
# p.39, the austin-putz Mrode chapter-3 R reproduction, and the Bioconductor
# GeneticsPed Mrode3.1 dataset) and were independently re-solved from the stated
# inputs (alpha = sigma_e2/sigma_a2 = 2). See hs_mrode_example_3_1_fixture() in
# R/validation-fixtures.R for the full provenance.
#
# X / Z / y / ids come from the payload builders (not the solver under test);
# Ainv is built by the tabular numerator-relationship method in pure base R
# (not nadiv, not the solver). All checks are pure R, fast, CI-runnable: no Julia
# engine, no skip guards.

test_that("the Henderson MME solver reproduces the published Mrode Example 3.1 EBVs", {
  fixture <- hsquared:::hs_mrode_example_3_1_fixture()
  spec <- hsquared:::hs_build_model_spec(
    fixture$formula,
    data = fixture$data,
    family = stats::gaussian(),
    REML = TRUE
  )
  payload <- hsquared:::hs_build_bridge_payload(spec)

  X <- as.matrix(payload$X)
  Z <- as.matrix(payload$Z)
  y <- payload$y
  ids <- payload$ids

  # Confirm the Example 3.1 design: 5 records on animals 4-8, 8 pedigree animals,
  # a 5x8 record->animal incidence Z, and the natural 1..8 animal order so the
  # published EBV vector aligns by id.
  expect_equal(ids, as.character(1:8))
  expect_equal(dim(X), c(5L, 2L))
  expect_equal(dim(Z), c(5L, 8L))
  expect_equal(length(y), 5L)

  reference <- hsquared:::hs_solve_henderson_mme_reference(
    y = y,
    X = X,
    Z = Z,
    Ainv = fixture$expected$Ainv,
    sigma_a2 = fixture$sigma_a2,
    sigma_e2 = fixture$sigma_e2,
    ids = ids
  )

  # The load-bearing check: the solver's EBVs must equal the PUBLISHED textbook
  # digits (aligned by id), not a value the solver itself produced.
  published <- fixture$expected$breeding_values[reference$breeding_values$id]
  expect_equal(
    reference$breeding_values$value,
    unname(published),
    tolerance = 1e-6
  )
  expect_equal(reference$breeding_values$id, ids)

  # The two sex solutions are parameterization-dependent (the sex block is
  # rank-deficient; the textbook uses a generalized inverse, the package uses
  # intercept + treatment contrast). The published invariant is the male - female
  # contrast, computed parameterization-free from a male and a female record's
  # fixed-effect prediction (the intercept cancels).
  male_row <- X[which(fixture$data$sex == "male")[1L], ]
  female_row <- X[which(fixture$data$sex == "female")[1L], ]
  contrast <- as.numeric((male_row - female_row) %*% reference$fixed_effects)
  expect_equal(
    contrast,
    fixture$expected$sex_contrast_male_minus_female,
    tolerance = 1e-6
  )
})

test_that("the published Mrode 3.1 anchor rejects perturbed EBVs (test of test)", {
  fixture <- hsquared:::hs_mrode_example_3_1_fixture()
  spec <- hsquared:::hs_build_model_spec(
    fixture$formula,
    data = fixture$data,
    family = stats::gaussian(),
    REML = TRUE
  )
  payload <- hsquared:::hs_build_bridge_payload(spec)

  reference <- hsquared:::hs_solve_henderson_mme_reference(
    y = payload$y,
    X = as.matrix(payload$X),
    Z = as.matrix(payload$Z),
    Ainv = fixture$expected$Ainv,
    sigma_a2 = fixture$sigma_a2,
    sigma_e2 = fixture$sigma_e2,
    ids = payload$ids
  )
  published <- unname(
    fixture$expected$breeding_values[reference$breeding_values$id]
  )

  # Sanity: the solver genuinely matches the published digits at the real
  # tolerance, so the rejection below is of a deliberately wrong value.
  expect_equal(reference$breeding_values$value, published, tolerance = 1e-6)

  # Test of test: published EBVs perturbed by +0.1 (far beyond 1e-6) must NOT
  # compare equal. If they did, the anchor above would be vacuous.
  expect_failure(expect_equal(
    reference$breeding_values$value,
    published + 0.1,
    tolerance = 1e-6
  ))
})
