# The twin's serialized Julia-native fitted univariate animal-model target
# (HSquared.jl #46): the engine fitting its own REML model on a 20-animal
# pedigree. The fixture is internally checked julia-free (the heritability and
# reliability formulas on the stored values), and a skip-guarded live test
# confirms the engine reproduces the serialized REML estimates end-to-end.

hs_fitted_target_csv <- function(file) {
  utils::read.csv(
    testthat::test_path("fixtures", "animal_model_fitted_target", file),
    stringsAsFactors = FALSE
  )
}

hs_fitted_target_meta <- function() {
  m <- hs_fitted_target_csv("expected_metadata.csv")
  stats::setNames(m$value, m$key)
}

test_that("the fitted-target fixture is internally consistent (h2 + reliability formulas)", {
  vc <- hs_fitted_target_csv("expected_variance_components.csv")
  sigma_a2 <- vc$value[vc$name == "sigma_a2"]
  sigma_e2 <- vc$value[vc$name == "sigma_e2"]
  meta <- hs_fitted_target_meta()

  # h2 = sigma_a2 / (sigma_a2 + sigma_e2).
  expect_equal(
    as.numeric(meta[["h2"]]),
    sigma_a2 / (sigma_a2 + sigma_e2),
    tolerance = 1e-10
  )

  # Animal-model reliability = 1 - PEV / (sigma_a2 * A_ii), where A_ii = 1 + F_i
  # is the diagonal of the relationship matrix (1 for non-inbred animals). PEV
  # and reliability are therefore mutually consistent only if the implied
  # A_ii = (PEV / sigma_a2) / (1 - reliability) are all >= 1 (inbreeding F >= 0).
  rel <- hs_fitted_target_csv("expected_reliability.csv")
  expect_true(all(rel$reliability >= 0 & rel$reliability <= 1))
  expect_true(all(rel$pev > 0))
  implied_a_ii <- (rel$pev / sigma_a2) / (1 - rel$reliability)
  expect_true(all(implied_a_ii >= 1 - 1e-8))
  # The least-related animals (the pedigree founders) are non-inbred, so the
  # smallest implied A_ii is exactly 1.
  expect_equal(min(implied_a_ii), 1, tolerance = 1e-8)
})

test_that("the live engine reproduces the serialized fitted-target REML estimates", {
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for the live fitted-target reproduction."
  )

  ped <- hs_fitted_target_csv("pedigree.csv")
  names(ped)[names(ped) == "animal"] <- "id"
  ped$sire[ped$sire == "0"] <- NA
  ped$dam[ped$dam == "0"] <- NA
  pheno <- hs_fitted_target_csv("phenotypes.csv")

  vc <- hs_fitted_target_csv("expected_variance_components.csv")
  sigma_a2 <- vc$value[vc$name == "sigma_a2"]
  sigma_e2 <- vc$value[vc$name == "sigma_e2"]
  beta <- hs_fitted_target_csv("expected_beta.csv")
  ebv <- hs_fitted_target_csv("expected_ebv.csv")
  rel <- hs_fitted_target_csv("expected_reliability.csv")
  meta <- hs_fitted_target_meta()

  fit <- hsquared(
    y ~ x + animal(1 | animal, pedigree = ped),
    data = pheno,
    family = stats::gaussian(),
    REML = TRUE
  )

  fvc <- variance_components(fit)
  expect_equal(
    fvc$estimate[fvc$component == "animal"],
    sigma_a2,
    tolerance = 1e-4
  )
  expect_equal(
    fvc$estimate[fvc$component == "residual"],
    sigma_e2,
    tolerance = 1e-4
  )
  expect_equal(
    heritability(fit)$estimate,
    as.numeric(meta[["h2"]]),
    tolerance = 1e-4
  )
  expect_equal(
    as.numeric(stats::logLik(fit)),
    as.numeric(meta[["loglik"]]),
    tolerance = 1e-3
  )

  # Fixed effects (Intercept, x) in model-matrix order.
  expect_equal(as.numeric(stats::coef(fit)), beta$value, tolerance = 1e-4)

  # Breeding values + reliabilities, aligned by animal id.
  bv <- merge(breeding_values(fit), ebv, by = "id")
  expect_equal(nrow(bv), nrow(ebv))
  expect_equal(bv$value.x, bv$value.y, tolerance = 1e-3)

  r <- merge(reliability(fit), rel, by = "id")
  expect_equal(nrow(r), nrow(rel))
  expect_equal(r$value, r$reliability, tolerance = 1e-3)
})
