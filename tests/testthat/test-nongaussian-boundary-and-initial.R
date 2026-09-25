# hsquared#222: the generic non-Gaussian bridge wrapper reads only `converged`
# off the Julia v0.9 envelope; it must also consume `boundary` (HSquared.jl#342
# / HSquared.jl#327), matching what `converged` already gets.
#
# hsquared#225: `target = "nongaussian"` never forwards `engine_control$initial`
# or `engine_control$restart_check` to `HSquared.fit_laplace_reml()`, so every
# live fit starts its search at the Julia-side hard-coded `sigma_a2 = 1.0` and a
# user cannot follow HSquared.jl#342's own retry advice ("retry with
# restart_check = true or a different initial").

test_that("hsquared#225: nongaussian honours `initial` and `restart_check`, not `max_dense_cells`", {
  honoured <- hsquared:::hs_engine_control_honoured_keys[["nongaussian"]]
  expect_true("initial" %in% honoured)
  expect_true("restart_check" %in% honoured)

  control <- hs_control(
    engine = "julia",
    engine_control = list(
      target = "nongaussian",
      initial = list(sigma_a2 = 0.05),
      restart_check = TRUE
    )
  )
  forwarded <- hsquared:::hs_engine_control_forwarding(control, "nongaussian")
  expect_equal(forwarded$initial, list(sigma_a2 = 0.05))
  expect_true(isTRUE(forwarded$restart_check))

  # `max_dense_cells` is not in the nongaussian honoured-keys row (checked
  # above), so supplying it must still error naming the target (hsquared#212).
  bad_control <- hs_control(
    engine = "julia",
    engine_control = list(target = "nongaussian", max_dense_cells = 1e5)
  )
  expect_error(
    hsquared:::hs_engine_control_forwarding(bad_control, "nongaussian"),
    class = "hsquared_unsupported_syntax"
  )
  cnd <- tryCatch(
    hsquared:::hs_engine_control_forwarding(bad_control, "nongaussian"),
    error = function(e) e
  )
  expect_match(conditionMessage(cnd), "max_dense_cells", fixed = TRUE)
  expect_match(conditionMessage(cnd), "nongaussian", fixed = TRUE)
})

test_that("hsquared#225: the v0.9 Julia command forwards `initial` and `restart_check` conditionally", {
  baseline <- hsquared:::hs_nongaussian_three_field_julia_command(
    family_symbol = "bernoulli",
    marginal = "laplace"
  )
  # An unsupplied `initial`/`restart_check` must emit no keyword at all, so the
  # `fit_laplace_reml(...)` CALL is byte-identical to the pre-fix one and the
  # Julia-side defaults (`sigma_a2 = 1.0`, `restart_check = false`) apply. The
  # full command asserted below is NOT byte-identical to pre-fix: the trailing
  # Dict deliberately gains "boundary"/"restart_estimate" (hsquared#222).
  expect_false(grepl("initial =", baseline, fixed = TRUE))
  expect_false(grepl("restart_check", baseline, fixed = TRUE))
  expect_identical(
    baseline,
    paste0(
      "hsq_ped = HSquared.normalize_pedigree(hsq_id, hsq_sire, hsq_dam); ",
      "hsq_Ainv = HSquared.pedigree_inverse(hsq_ped); ",
      "hsq_fit = HSquared.fit_laplace_reml(",
      "hsq_y, hsq_X, hsq_Z, hsq_Ainv; ",
      "family = Symbol(hsq_family), marginal = Symbol(hsq_marginal), ",
      "ids = hsq_ped.ids, iterations = hsq_iterations); ",
      "hsq_result = HSquared.nongaussian_three_field_payload(",
      "hsq_fit; predictor_variance = 0.0, response_length = length(hsq_y)); ",
      "hsq_ng_raw = Dict(",
      "\"schema\" => hsq_result.schema, ",
      "\"family\" => hsq_result.family, ",
      "\"method\" => hsq_result.method, ",
      "\"loglik\" => hsq_result.loglik, ",
      "\"components\" => Dict(\"names\" => [\"V_A\", \"V_RE\", \"V_O\"], ",
      "\"values\" => [hsq_result.components.V_A, hsq_result.components.V_RE, hsq_result.components.V_O]), ",
      "\"fixed_effects\" => Dict(\"names\" => hsq_result.fixed_effects.names, ",
      "\"values\" => hsq_result.fixed_effects.values), ",
      "\"breeding_ids\" => string.(collect(hsq_fit.ids)), ",
      "\"breeding_values\" => collect(Float64, hsq_fit.breeding_values), ",
      "\"h2_latent\" => hsq_result.h2_latent, ",
      "\"h2_liability\" => hsq_result.h2_liability, ",
      "\"h2_observation\" => hsq_result.h2_observation, ",
      "\"h2_observation_undefined_reason\" => hsq_result.h2_observation_undefined_reason, ",
      "\"n_trials\" => hsq_result.n_trials, ",
      "\"converged\" => hsq_fit.converged, ",
      "\"boundary\" => hsq_fit.boundary, ",
      "\"restart_estimate\" => hsq_fit.restart_estimate);"
    )
  )

  with_initial <- hsquared:::hs_nongaussian_three_field_julia_command(
    family_symbol = "bernoulli",
    marginal = "laplace",
    initial = 0.05
  )
  expect_match(with_initial, "initial = \\(sigma_a2 = 0\\.05,\\)", fixed = FALSE)
  expect_false(grepl("restart_check", with_initial, fixed = TRUE))

  with_restart <- hsquared:::hs_nongaussian_three_field_julia_command(
    family_symbol = "bernoulli",
    marginal = "laplace",
    restart_check = TRUE
  )
  expect_match(with_restart, "restart_check = true", fixed = TRUE)
  expect_false(grepl("initial =", with_restart, fixed = TRUE))

  with_both <- hsquared:::hs_nongaussian_three_field_julia_command(
    family_symbol = "bernoulli",
    marginal = "laplace",
    initial = 0.05,
    restart_check = TRUE
  )
  expect_match(with_both, "initial = \\(sigma_a2 = 0\\.05,\\)", fixed = FALSE)
  expect_match(with_both, "restart_check = true", fixed = TRUE)
})

test_that("hsquared#225: `initial` and `restart_check` validators reject malformed input", {
  expect_null(hsquared:::hs_validate_nongaussian_initial(NULL))
  expect_equal(
    hsquared:::hs_validate_nongaussian_initial(list(sigma_a2 = 0.05)),
    0.05
  )
  expect_error(
    hsquared:::hs_validate_nongaussian_initial(list(sigma_e2 = 0.05)),
    "sigma_a2"
  )
  expect_error(
    hsquared:::hs_validate_nongaussian_initial(list(sigma_a2 = -1)),
    "positive"
  )
  expect_error(
    hsquared:::hs_validate_nongaussian_initial(list(sigma_a2 = c(0.1, 0.2))),
    "single"
  )

  expect_false(hsquared:::hs_validate_restart_check(FALSE))
  expect_true(hsquared:::hs_validate_restart_check(TRUE))
  expect_error(
    hsquared:::hs_validate_restart_check(NA),
    "TRUE/FALSE"
  )
  expect_error(
    hsquared:::hs_validate_restart_check(c(TRUE, FALSE)),
    "TRUE/FALSE"
  )
})

# hsquared#222's `hs_ng09_boundary()`/normalizer tests live in
# test-nongaussian-three-field-v09.R, next to `ng09_raw()` and the rest of the
# v0.9 envelope contract tests they depend on.
