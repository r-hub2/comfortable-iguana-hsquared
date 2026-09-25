# A15 engine = "julia" smoke harness (B4 bridge prep).
# Consolidates live smokes S2–S3 from ~/local-scratch/h2-b4-bridge-plan.md.
# S1 (default ai_reml) remains in test-julia-bridge.R; S4 gryphon in
# test-validation-fixtures.R.

hs_engine_julia_smoke_skip <- function() {
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for engine=julia smoke."
  )
}

hs_engine_julia_mrode_fixture <- function() {
  hsquared:::hs_mrode_supplied_variance_validation_fixture()
}

test_that("S2: engine = julia dense fit_animal_model differs from default ai_reml on Mrode fixture", {
  hs_engine_julia_smoke_skip()

  fixture <- hs_engine_julia_mrode_fixture()
  fit_args <- list(
    formula = fixture$formula,
    data = fixture$data,
    family = stats::gaussian(),
    REML = TRUE
  )

  default_fit <- do.call(hsquared, fit_args)

  dense_fit <- do.call(
    hsquared,
    c(
      fit_args,
      list(
        control = hs_control(
          engine = "julia",
          engine_control = list(
            initial = c(sigma_a2 = 1, sigma_e2 = 1)
          )
        )
      )
    )
  )

  expect_s3_class(default_fit, "hsquared_fit")
  expect_s3_class(dense_fit, "hsquared_fit")
  expect_equal(default_fit$spec$target, "ai_reml")

  default_diag <- fit_diagnostics(default_fit)
  dense_diag <- fit_diagnostics(dense_fit)
  expect_equal(default_diag$value[default_diag$metric == "target"], "ai_reml")
  expect_equal(
    dense_diag$value[dense_diag$metric == "dense_validation_path"],
    "TRUE"
  )
  expect_false(
    identical(
      dense_diag$value[dense_diag$metric == "target"],
      "ai_reml"
    )
  )

  default_vc <- variance_components(default_fit)$estimate
  dense_vc <- variance_components(dense_fit)$estimate
  expect_equal(variance_components(default_fit)$component, c("animal", "residual"))
  expect_equal(variance_components(dense_fit)$component, c("animal", "residual"))
  expect_true(all(is.finite(default_vc)) && all(default_vc > 0))
  expect_true(all(is.finite(dense_vc)) && all(dense_vc > 0))
  expect_true(is.finite(stats::logLik(default_fit)))
  expect_true(is.finite(stats::logLik(dense_fit)))

  default_h2 <- heritability(default_fit)$estimate
  dense_h2 <- heritability(dense_fit)$estimate
  expect_true(default_h2 > 0 && default_h2 < 1)
  expect_true(dense_h2 > 0 && dense_h2 < 1)

  # The two paths are different estimators; that is asserted structurally above
  # (dense_validation_path, target != ai_reml), never by requiring the numbers
  # to disagree — a sharper future optimizer must not turn this test red.
  # Bounds carry ~5x headroom over the 2026-09-01 measurement on this fixture
  # (max relative VC delta 9.9e-4, h2 5.2e-5, logLik 1.3e-8); see the A15
  # check-log entry and tolerance_rule in bridge-parity-smoke-status.tsv.
  expect_lt(max(abs(default_vc - dense_vc) / abs(default_vc)), 5e-3)
  expect_lt(abs(default_h2 - dense_h2), 1e-3)
  expect_lt(
    abs(as.numeric(stats::logLik(default_fit)) - as.numeric(stats::logLik(dense_fit))),
    1e-5
  )
})

test_that("S3: engine = julia target = ai_reml matches the default fit path on Mrode fixture", {
  hs_engine_julia_smoke_skip()

  fixture <- hs_engine_julia_mrode_fixture()
  fit_args <- list(
    formula = fixture$formula,
    data = fixture$data,
    family = stats::gaussian(),
    REML = TRUE
  )

  default_fit <- do.call(hsquared, fit_args)
  explicit_fit <- do.call(
    hsquared,
    c(
      fit_args,
      list(
        control = hs_control(
          engine = "julia",
          engine_control = list(
            target = "ai_reml",
            initial = c(sigma_a2 = 1, sigma_e2 = 1),
            iterations = 200L
          )
        )
      )
    )
  )

  expect_s3_class(default_fit, "hsquared_fit")
  expect_s3_class(explicit_fit, "hsquared_fit")
  expect_equal(default_fit$spec$target, "ai_reml")
  expect_equal(explicit_fit$spec$target, "ai_reml")

  default_vc <- variance_components(default_fit)$estimate
  explicit_vc <- variance_components(explicit_fit)$estimate
  expect_equal(explicit_vc, default_vc, tolerance = 1e-5)
  expect_equal(
    heritability(explicit_fit)$estimate,
    heritability(default_fit)$estimate,
    tolerance = 1e-5
  )
  expect_equal(
    stats::logLik(explicit_fit),
    stats::logLik(default_fit),
    tolerance = 1e-5
  )
})
