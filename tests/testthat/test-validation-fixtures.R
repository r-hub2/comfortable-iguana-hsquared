test_that("tiny animal validation fixture pins R payload ordering", {
  fixture <- hsquared:::hs_tiny_animal_validation_fixture()
  spec <- hsquared:::hs_build_model_spec(
    fixture$formula,
    data = fixture$data,
    family = stats::gaussian(),
    REML = TRUE
  )
  payload <- hsquared:::hs_build_bridge_payload(spec)

  expect_equal(payload$ids, fixture$expected$ids)
  expect_equal(payload$pedigree$id, fixture$expected$ids)
  expect_equal(payload$pedigree$sire_index, fixture$expected$sire_index)
  expect_equal(payload$pedigree$dam_index, fixture$expected$dam_index)
  expect_equal(payload$metadata$observed_ids, fixture$data$id)
  expect_equal(payload$metadata$observed_id_index, c(1L, 2L, 3L))
  expect_equal(unname(as.matrix(payload$Z)), fixture$expected$Z)
})

test_that("tiny animal validation fixture matches Julia Ainv when available", {
  fixture <- hsquared:::hs_tiny_animal_validation_fixture()
  spec <- hsquared:::hs_build_model_spec(
    fixture$formula,
    data = fixture$data,
    family = stats::gaussian(),
    REML = TRUE
  )
  payload <- hsquared:::hs_build_bridge_payload(spec)

  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for live Ainv validation."
  )

  hs_skip_live_julia()
  hsquared:::hs_julia_setup(hsquared:::hs_default_julia_project())
  JuliaCall::julia_assign("hsq_val_id", payload$pedigree$id)
  JuliaCall::julia_assign(
    "hsq_val_sire",
    hsquared:::hs_parent_for_julia(payload$pedigree$sire)
  )
  JuliaCall::julia_assign(
    "hsq_val_dam",
    hsquared:::hs_parent_for_julia(payload$pedigree$dam)
  )
  JuliaCall::julia_command(paste(
    "hsq_val_ped = HSquared.normalize_pedigree(",
    "hsq_val_id, hsq_val_sire, hsq_val_dam);",
    "hsq_val_Ainv = Matrix(HSquared.pedigree_inverse(hsq_val_ped));"
  ))

  observed <- JuliaCall::julia_eval("hsq_val_Ainv")
  dimnames(observed) <- dimnames(fixture$expected$Ainv)

  expect_equal(observed, fixture$expected$Ainv, tolerance = 1e-12)
})

test_that("Henderson MME validation fixture pins supplied-variance solutions", {
  fixture <- hsquared:::hs_henderson_mme_validation_fixture()
  spec <- hsquared:::hs_build_model_spec(
    fixture$formula,
    data = fixture$data,
    family = stats::gaussian(),
    REML = FALSE
  )
  payload <- hsquared:::hs_build_bridge_payload(spec)

  expect_equal(payload$ids, fixture$expected$ids)
  expect_equal(payload$pedigree$id, fixture$expected$ids)
  expect_equal(
    as.matrix(fixture$expected$Ainv),
    matrix(
      c(
        2,
        1,
        -1,
        -1,
        0,
        1,
        2,
        -1,
        -1,
        0,
        -1,
        -1,
        2.5,
        0.5,
        -1,
        -1,
        -1,
        0.5,
        2.5,
        -1,
        0,
        0,
        -1,
        -1,
        2
      ),
      nrow = 5,
      byrow = TRUE,
      dimnames = list(fixture$expected$ids, fixture$expected$ids)
    )
  )

  reference <- hsquared:::hs_solve_henderson_mme_reference(
    y = payload$y,
    X = payload$X,
    Z = payload$Z,
    Ainv = fixture$expected$Ainv,
    sigma_a2 = fixture$sigma_a2,
    sigma_e2 = fixture$sigma_e2,
    ids = payload$ids
  )

  expect_equal(
    unname(reference$fixed_effects),
    unname(fixture$expected$fixed_effects),
    tolerance = 1e-12
  )
  expect_equal(
    reference$breeding_values,
    fixture$expected$breeding_values,
    tolerance = 1e-12
  )
  expect_equal(reference$fitted, fixture$expected$fitted, tolerance = 1e-12)
  expect_equal(reference$heritability, fixture$expected$heritability)
})

test_that("Henderson MME fixture matches Julia henderson_mme when available", {
  fixture <- hsquared:::hs_henderson_mme_validation_fixture()
  spec <- hsquared:::hs_build_model_spec(
    fixture$formula,
    data = fixture$data,
    family = stats::gaussian(),
    REML = FALSE
  )
  payload <- hsquared:::hs_build_bridge_payload(spec)

  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for live Henderson MME validation."
  )

  hs_skip_live_julia()
  hsquared:::hs_julia_setup(hsquared:::hs_default_julia_project())
  hsquared:::hs_julia_assign_payload(
    payload,
    initial = c(
      sigma_a2 = fixture$sigma_a2,
      sigma_e2 = fixture$sigma_e2
    )
  )
  JuliaCall::julia_assign("hsq_mme_sigma_a2", fixture$sigma_a2)
  JuliaCall::julia_assign("hsq_mme_sigma_e2", fixture$sigma_e2)
  JuliaCall::julia_command(paste(
    "hsq_mme_ped = HSquared.normalize_pedigree(hsq_id, hsq_sire, hsq_dam);",
    "hsq_mme_Ainv = HSquared.pedigree_inverse(hsq_mme_ped);",
    "hsq_mme_spec = HSquared.animal_model_spec(",
    "hsq_y, hsq_X, hsq_Z, hsq_mme_Ainv;",
    "ids = hsq_mme_ped.ids, method = :ML);",
    "hsq_mme = HSquared.henderson_mme(",
    "hsq_mme_spec, hsq_mme_sigma_a2, hsq_mme_sigma_e2);",
    "hsq_mme_bv = HSquared.breeding_values(hsq_mme);",
    "hsq_mme_payload = Dict(",
    "\"fixed_effects\" => HSquared.fixed_effects(hsq_mme),",
    "\"animal_ids\" => hsq_mme_bv.ids,",
    "\"animal_effects\" => hsq_mme_bv.values,",
    "\"fitted\" => HSquared.fitted_values(hsq_mme)",
    ");"
  ))

  observed <- JuliaCall::julia_eval("hsq_mme_payload")

  expect_equal(
    as.numeric(observed$fixed_effects),
    unname(fixture$expected$fixed_effects),
    tolerance = 1e-10
  )
  expect_equal(
    as.character(observed$animal_ids),
    fixture$expected$breeding_values$id
  )
  expect_equal(
    as.numeric(observed$animal_effects),
    fixture$expected$breeding_values$value,
    tolerance = 1e-10
  )
  expect_equal(
    as.numeric(observed$fitted),
    fixture$expected$fitted$.fitted,
    tolerance = 1e-10
  )
})

test_that("sparse REML likelihood fixture pins closed-form targets", {
  fixture <- hsquared:::hs_reml_likelihood_validation_fixture()
  spec <- hsquared:::hs_build_model_spec(
    fixture$formula,
    data = fixture$data,
    family = stats::gaussian(),
    REML = TRUE
  )
  payload <- hsquared:::hs_build_bridge_payload(spec)

  expect_equal(payload$ids, fixture$expected$ids)
  expect_equal(payload$pedigree$id, fixture$expected$ids)
  expect_equal(payload$pedigree$sire_index, fixture$expected$sire_index)
  expect_equal(payload$pedigree$dam_index, fixture$expected$dam_index)
  expect_equal(as.numeric(payload$X), rep(1, 3))
  expect_equal(unname(as.matrix(payload$Z)), fixture$expected$Z)
  expect_equal(fixture$expected$Ainv, diag(3))
  expect_equal(unname(fixture$expected$fixed_effects), 2)
  expect_equal(
    fixture$expected$ml_loglik,
    -0.5 * (3 * log(2 * pi) + 3 * log(2) + 1)
  )
  expect_equal(
    fixture$expected$reml_loglik,
    -0.5 * (2 * log(2 * pi) + 3 * log(2) + log(1.5) + 1)
  )
})

test_that("sparse REML likelihood fixture matches Julia dense REML when available", {
  fixture <- hsquared:::hs_reml_likelihood_validation_fixture()
  spec <- hsquared:::hs_build_model_spec(
    fixture$formula,
    data = fixture$data,
    family = stats::gaussian(),
    REML = TRUE
  )
  payload <- hsquared:::hs_build_bridge_payload(spec)

  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    paste(
      "JuliaCall, Julia, and local HSquared.jl are required for live sparse",
      "REML likelihood validation."
    )
  )

  hs_skip_live_julia()
  hsquared:::hs_julia_setup(hsquared:::hs_default_julia_project())
  hsquared:::hs_julia_assign_payload(
    payload,
    initial = c(
      sigma_a2 = fixture$sigma_a2,
      sigma_e2 = fixture$sigma_e2
    )
  )
  JuliaCall::julia_assign("hsq_reml_sigma_a2", fixture$sigma_a2)
  JuliaCall::julia_assign("hsq_reml_sigma_e2", fixture$sigma_e2)
  JuliaCall::julia_command(paste(
    "hsq_reml_ped = HSquared.normalize_pedigree(hsq_id, hsq_sire, hsq_dam);",
    "hsq_reml_Ainv = HSquared.pedigree_inverse(hsq_reml_ped);",
    "hsq_reml_spec = HSquared.animal_model_spec(",
    "hsq_y, hsq_X, hsq_Z, hsq_reml_Ainv;",
    "ids = hsq_reml_ped.ids, method = :REML);",
    "hsq_dense_reml = HSquared.gaussian_loglik(",
    "hsq_reml_spec, hsq_reml_sigma_a2, hsq_reml_sigma_e2;",
    "method = :REML);",
    "hsq_sparse_reml = HSquared.sparse_reml_loglik(",
    "hsq_reml_spec, hsq_reml_sigma_a2, hsq_reml_sigma_e2);",
    "hsq_dense_ml = HSquared.gaussian_loglik(",
    "hsq_reml_spec, hsq_reml_sigma_a2, hsq_reml_sigma_e2;",
    "method = :ML);",
    "hsq_reml_payload = Dict(",
    "\"dense_reml_loglik\" => hsq_dense_reml.loglik,",
    "\"sparse_reml_loglik\" => hsq_sparse_reml.loglik,",
    "\"dense_ml_loglik\" => hsq_dense_ml.loglik,",
    "\"dense_beta\" => hsq_dense_reml.beta,",
    "\"sparse_beta\" => hsq_sparse_reml.beta",
    ");"
  ))

  observed <- JuliaCall::julia_eval("hsq_reml_payload")

  expect_equal(
    as.numeric(observed$dense_reml_loglik),
    fixture$expected$reml_loglik,
    tolerance = 1e-10
  )
  expect_equal(
    as.numeric(observed$sparse_reml_loglik),
    as.numeric(observed$dense_reml_loglik),
    tolerance = 1e-10
  )
  expect_equal(
    as.numeric(observed$dense_ml_loglik),
    fixture$expected$ml_loglik,
    tolerance = 1e-10
  )
  expect_equal(
    as.numeric(observed$dense_beta),
    unname(fixture$expected$fixed_effects),
    tolerance = 1e-10
  )
  expect_equal(
    as.numeric(observed$sparse_beta),
    unname(fixture$expected$fixed_effects),
    tolerance = 1e-10
  )
})

test_that("Mrode-style supplied-variance fixture pins R reference outputs", {
  fixture <- hsquared:::hs_mrode_supplied_variance_validation_fixture()
  spec <- hsquared:::hs_build_model_spec(
    fixture$formula,
    data = fixture$data,
    family = stats::gaussian(),
    REML = FALSE
  )
  payload <- hsquared:::hs_build_bridge_payload(spec)

  expect_equal(fixture$name, "mrode_style_supplied_variance_outputs")
  expect_equal(payload$ids, fixture$expected$ids)
  expect_equal(payload$pedigree$id, fixture$expected$ids)
  expect_equal(as.numeric(payload$X[, 1]), rep(1, 12))
  expect_equal(as.numeric(payload$X[, 2]), fixture$data$x)
  expect_equal(unname(as.matrix(payload$Z)), diag(12))
  expect_equal(fixture$expected$Ainv, t(fixture$expected$Ainv))

  reference <- hsquared:::hs_solve_henderson_mme_reference(
    y = payload$y,
    X = payload$X,
    Z = payload$Z,
    Ainv = fixture$expected$Ainv,
    sigma_a2 = fixture$sigma_a2,
    sigma_e2 = fixture$sigma_e2,
    ids = payload$ids
  )
  ml <- hsquared:::hs_gaussian_loglik_reference(
    y = payload$y,
    X = payload$X,
    Z = payload$Z,
    Ainv = fixture$expected$Ainv,
    sigma_a2 = fixture$sigma_a2,
    sigma_e2 = fixture$sigma_e2,
    method = "ML"
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

  expect_equal(
    reference$fixed_effects,
    fixture$expected$fixed_effects,
    tolerance = 1e-10,
    ignore_attr = TRUE
  )
  expect_equal(
    reference$breeding_values,
    fixture$expected$breeding_values,
    tolerance = 1e-10
  )
  expect_equal(reference$fitted, fixture$expected$fitted, tolerance = 1e-10)
  expect_equal(
    reference$prediction_error_variance,
    fixture$expected$prediction_error_variance,
    tolerance = 1e-10
  )
  expect_equal(
    reference$reliability,
    fixture$expected$reliability,
    tolerance = 1e-10
  )
  expect_equal(reference$heritability, fixture$expected$heritability)
  expect_equal(ml$loglik, fixture$expected$ml_loglik, tolerance = 1e-10)
  expect_equal(reml$loglik, fixture$expected$reml_loglik, tolerance = 1e-10)
  expect_equal(
    reml$beta,
    unname(fixture$expected$fixed_effects),
    tolerance = 1e-10
  )
})

test_that("Mrode-style supplied-variance fixture matches Julia when available", {
  fixture <- hsquared:::hs_mrode_supplied_variance_validation_fixture()
  spec <- hsquared:::hs_build_model_spec(
    fixture$formula,
    data = fixture$data,
    family = stats::gaussian(),
    REML = FALSE
  )
  payload <- hsquared:::hs_build_bridge_payload(spec)

  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    paste(
      "JuliaCall, Julia, and local HSquared.jl are required for live",
      "Mrode-style supplied-variance validation."
    )
  )

  hs_skip_live_julia()
  hsquared:::hs_julia_setup(hsquared:::hs_default_julia_project())
  hsquared:::hs_julia_assign_payload(
    payload,
    initial = c(
      sigma_a2 = fixture$sigma_a2,
      sigma_e2 = fixture$sigma_e2
    )
  )
  JuliaCall::julia_assign("hsq_ms_sigma_a2", fixture$sigma_a2)
  JuliaCall::julia_assign("hsq_ms_sigma_e2", fixture$sigma_e2)
  JuliaCall::julia_command(paste(
    "hsq_ms_ped = HSquared.normalize_pedigree(hsq_id, hsq_sire, hsq_dam);",
    "hsq_ms_Ainv = HSquared.pedigree_inverse(hsq_ms_ped);",
    "hsq_ms_spec_ml = HSquared.animal_model_spec(",
    "hsq_y, hsq_X, hsq_Z, hsq_ms_Ainv;",
    "ids = hsq_ms_ped.ids, method = :ML);",
    "hsq_ms_spec_reml = HSquared.animal_model_spec(",
    "hsq_y, hsq_X, hsq_Z, hsq_ms_Ainv;",
    "ids = hsq_ms_ped.ids, method = :REML);",
    "hsq_ms_ml = HSquared.gaussian_loglik(",
    "hsq_ms_spec_ml, hsq_ms_sigma_a2, hsq_ms_sigma_e2;",
    "method = :ML);",
    "hsq_ms_reml = HSquared.gaussian_loglik(",
    "hsq_ms_spec_reml, hsq_ms_sigma_a2, hsq_ms_sigma_e2;",
    "method = :REML);",
    "hsq_ms_sparse_reml = HSquared.sparse_reml_loglik(",
    "hsq_ms_spec_reml, hsq_ms_sigma_a2, hsq_ms_sigma_e2);",
    "hsq_ms_mme = HSquared.henderson_mme(",
    "hsq_ms_spec_ml, hsq_ms_sigma_a2, hsq_ms_sigma_e2);",
    "hsq_ms_bv = HSquared.breeding_values(hsq_ms_mme);",
    "hsq_ms_payload = Dict(",
    "\"ids\" => hsq_ms_ped.ids,",
    "\"Ainv\" => Matrix(hsq_ms_Ainv),",
    "\"ml_loglik\" => hsq_ms_ml.loglik,",
    "\"dense_reml_loglik\" => hsq_ms_reml.loglik,",
    "\"sparse_reml_loglik\" => hsq_ms_sparse_reml.loglik,",
    "\"sparse_beta\" => hsq_ms_sparse_reml.beta,",
    "\"fixed_effects\" => HSquared.fixed_effects(hsq_ms_mme),",
    "\"animal_ids\" => hsq_ms_bv.ids,",
    "\"animal_effects\" => hsq_ms_bv.values,",
    "\"fitted\" => HSquared.fitted_values(hsq_ms_mme),",
    "\"heritability\" => HSquared.heritability(hsq_ms_mme)",
    ");",
    "if isdefined(HSquared, :prediction_error_variance) &&",
    "isdefined(HSquared, :reliability) &&",
    "applicable(HSquared.prediction_error_variance, hsq_ms_mme) &&",
    "applicable(HSquared.reliability, hsq_ms_mme);",
    "hsq_ms_pev = HSquared.prediction_error_variance(hsq_ms_mme);",
    "hsq_ms_rel = HSquared.reliability(hsq_ms_mme);",
    "hsq_ms_payload[\"pev_ids\"] = hsq_ms_pev.ids;",
    "hsq_ms_payload[\"pev\"] = hsq_ms_pev.values;",
    "hsq_ms_payload[\"reliability_ids\"] = hsq_ms_rel.ids;",
    "hsq_ms_payload[\"reliability\"] = hsq_ms_rel.values;",
    "end;"
  ))

  observed <- JuliaCall::julia_eval("hsq_ms_payload")

  expect_equal(as.character(observed$ids), fixture$expected$ids)
  dimnames(observed$Ainv) <- list(fixture$expected$ids, fixture$expected$ids)
  expect_equal(observed$Ainv, fixture$expected$Ainv, tolerance = 1e-10)
  expect_equal(
    as.numeric(observed$fixed_effects),
    unname(fixture$expected$fixed_effects),
    tolerance = 1e-10
  )
  expect_equal(as.character(observed$animal_ids), fixture$expected$ids)
  expect_equal(
    as.numeric(observed$animal_effects),
    fixture$expected$breeding_values$value,
    tolerance = 1e-10
  )
  expect_equal(
    as.numeric(observed$fitted),
    fixture$expected$fitted$.fitted,
    tolerance = 1e-10
  )
  expect_equal(
    as.numeric(observed$ml_loglik),
    fixture$expected$ml_loglik,
    tolerance = 1e-10
  )
  expect_equal(
    as.numeric(observed$dense_reml_loglik),
    fixture$expected$reml_loglik,
    tolerance = 1e-10
  )
  expect_equal(
    as.numeric(observed$sparse_reml_loglik),
    fixture$expected$reml_loglik,
    tolerance = 1e-10
  )
  expect_equal(
    as.numeric(observed$sparse_beta),
    unname(fixture$expected$fixed_effects),
    tolerance = 1e-10
  )
  expect_equal(
    as.numeric(observed$heritability),
    fixture$expected$heritability,
    tolerance = 1e-12
  )

  if (!is.null(observed$pev) && !is.null(observed$reliability)) {
    expect_equal(as.character(observed$pev_ids), fixture$expected$ids)
    expect_equal(as.character(observed$reliability_ids), fixture$expected$ids)
    expect_equal(
      as.numeric(observed$pev),
      fixture$expected$prediction_error_variance$value,
      tolerance = 1e-10
    )
    expect_equal(
      as.numeric(observed$reliability),
      fixture$expected$reliability$value,
      tolerance = 1e-10
    )
  }
})

test_that("sparse REML optimizer reaches the same REML optimum from different starts", {
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for live sparse REML estimate-recovery validation."
  )

  fixture <- hsquared:::hs_mrode_supplied_variance_validation_fixture()
  fit_from <- function(init) {
    hsquared(
      fixture$formula,
      data = fixture$data,
      family = stats::gaussian(),
      REML = TRUE,
      control = hs_control(
        engine = "julia",
        engine_control = list(
          target = "sparse_reml",
          initial = init,
          iterations = 1000L
        )
      )
    )
  }

  # Estimate-recovery discipline (see 04-validation-canon comparator rule): the
  # REML objective optimum is start-independent. This validates the experimental
  # optimizer, NOT data-generating recovery, supplied-truth recovery, or ASReml
  # parity. Both fits optimize the SAME estimand (the REML objective).
  fit_a <- fit_from(c(sigma_a2 = 0.5, sigma_e2 = 0.5))
  fit_b <- fit_from(c(sigma_a2 = 3.0, sigma_e2 = 1.5))

  la <- fit_a$result$loglik
  lb <- fit_b$result$loglik
  expect_true(is.finite(la) && is.finite(lb))
  # Same REML optimum from different starts (loglik is flat near the optimum, so
  # compare it tightly and the variance estimates a little more loosely).
  expect_equal(la, lb, tolerance = 1e-3)
  expect_equal(
    variance_components(fit_a)$estimate,
    variance_components(fit_b)$estimate,
    tolerance = 1e-2
  )
  expect_true(all(variance_components(fit_a)$estimate > 0))
  expect_equal(
    fit_diagnostics(fit_a)$value[
      fit_diagnostics(fit_a)$metric == "variance_components_source"
    ],
    "estimated_sparse_reml"
  )
})

test_that("sparse and dense REML optimizers reach the same REML optimum", {
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for live sparse-vs-dense REML optimizer validation."
  )

  fixture <- hsquared:::hs_mrode_supplied_variance_validation_fixture()
  init <- c(sigma_a2 = 1, sigma_e2 = 1)
  fit_with <- function(extra) {
    hsquared(
      fixture$formula,
      data = fixture$data,
      family = stats::gaussian(),
      REML = TRUE,
      control = hs_control(
        engine = "julia",
        engine_control = c(list(initial = init), extra)
      )
    )
  }

  # Internal comparator: the dense REML optimizer (fit_variance_components) and
  # the sparse REML optimizer (fit_sparse_reml) maximize the SAME REML objective
  # via different linear algebra, so on the same data they must reach the same
  # optimum. This cross-validates the sparse optimizer against the dense one; it
  # is NOT an external comparator, DGP recovery, or production-fitting claim.
  # hsquared#212: the implicit `target = "fit_animal_model"` dense path has no
  # `iterations` control (it never did -- `hs_fit_julia_payload()` has no such
  # formal), so `iterations` is supplied only on the sparse_reml side, which
  # does honour it; this changes nothing about the dense fit's actual
  # behaviour, only whether the call passes the new #212 key-forwarding gate.
  dense <- fit_with(list())
  sparse <- fit_with(list(target = "sparse_reml", iterations = 1000L))

  expect_equal(dense$spec$target %||% "fit_animal_model", "fit_animal_model")
  expect_equal(sparse$spec$target, "sparse_reml")
  expect_true(is.finite(dense$result$loglik) && is.finite(sparse$result$loglik))
  # Same REML optimum (loglik is the shared objective; compare tightly).
  expect_equal(dense$result$loglik, sparse$result$loglik, tolerance = 1e-3)
  # Same variance estimates (flat near the optimum; compare a little loosely).
  expect_equal(
    variance_components(dense)$estimate,
    variance_components(sparse)$estimate,
    tolerance = 5e-2
  )
})

test_that("AI-REML and sparse REML optimizers reach the same REML optimum", {
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for live AI-REML-vs-sparse REML optimizer validation."
  )

  fixture <- hsquared:::hs_mrode_supplied_variance_validation_fixture()
  init <- c(sigma_a2 = 1, sigma_e2 = 1)
  fit_target <- function(target) {
    hsquared(
      fixture$formula,
      data = fixture$data,
      family = stats::gaussian(),
      REML = TRUE,
      control = hs_control(
        engine = "julia",
        engine_control = list(
          target = target,
          initial = init,
          iterations = 1000L
        )
      )
    )
  }

  # Internal comparator: average-information REML (fit_ai_reml) and the sparse
  # NelderMead REML optimizer (fit_sparse_reml) maximize the SAME REML objective
  # via different algorithms (an AI/Newton step vs derivative-free search), so on
  # the same data they must reach the same optimum. This cross-validates the
  # AI-REML estimator against the sparse one; it is NOT an external comparator,
  # DGP recovery, or production-fitting claim.
  ai <- fit_target("ai_reml")
  sparse <- fit_target("sparse_reml")

  expect_equal(ai$spec$target, "ai_reml")
  expect_equal(sparse$spec$target, "sparse_reml")
  expect_true(is.finite(ai$result$loglik) && is.finite(sparse$result$loglik))
  # Same REML optimum (loglik is the shared objective; compare tightly).
  expect_equal(ai$result$loglik, sparse$result$loglik, tolerance = 1e-3)
  # Same variance estimates (flat near the optimum; compare a little loosely).
  expect_equal(
    variance_components(ai)$estimate,
    variance_components(sparse)$estimate,
    tolerance = 5e-2
  )
})

test_that("independent pure-R REML optimizer matches the Julia sparse REML estimate", {
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

  # Independent pure-R REML optimization (no Julia) — always runs, including CI.
  expect_equal(ref$convergence, 0L)
  expect_true(all(is.finite(ref$estimate)) && all(ref$estimate > 0))
  expect_true(is.finite(ref$loglik))

  # Cross-check: the Julia sparse REML estimate matches the independent pure-R
  # optimum (same estimand, fully independent implementation). Skip-guarded.
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for the live sparse REML cross-check."
  )
  fit <- hsquared(
    fixture$formula,
    data = fixture$data,
    family = stats::gaussian(),
    REML = TRUE,
    control = hs_control(
      engine = "julia",
      engine_control = list(
        target = "sparse_reml",
        initial = c(sigma_a2 = 1, sigma_e2 = 1),
        iterations = 1000L
      )
    )
  )
  expect_equal(
    unname(variance_components(fit)$estimate),
    unname(ref$estimate),
    tolerance = 5e-2
  )
})

test_that("hsquared's REML solution is at least as good as the pedigreemm comparator", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("pedigreemm")
  testthat::skip_if_not_installed("withr")

  fixture <- hsquared:::hs_replicated_animal_comparator_fixture()
  lab <- as.character(fixture$pedigree$id)

  # External comparator: pedigreemm (lme4-based pedigree animal model, REML).
  # pedigreemm Depends on lme4, which must be attached for its internal lmer().
  withr::local_package("pedigreemm")
  ped <- pedigreemm::pedigree(
    sire = match(as.character(fixture$pedigree$sire), lab),
    dam = match(as.character(fixture$pedigree$dam), lab),
    label = lab
  )
  dat <- fixture$data
  dat$id <- factor(as.character(dat$id), levels = lab)
  m <- pedigreemm::pedigreemm(
    y ~ x + (1 | id),
    data = dat,
    pedigree = list(id = ped)
  )
  vc <- as.data.frame(lme4::VarCorr(m))
  ext <- c(
    sigma_a2 = vc$vcov[vc$grp == "id"],
    sigma_e2 = vc$vcov[vc$grp == "Residual"]
  )

  # hsquared's REML estimate via the independent pure-R reference (which equals
  # the Julia sparse_reml estimate; see the pure-R cross-check test above).
  spec <- hsquared:::hs_build_model_spec(
    fixture$formula,
    data = fixture$data,
    family = stats::gaussian(),
    REML = TRUE
  )
  payload <- hsquared:::hs_build_bridge_payload(spec)
  Z <- as.matrix(payload$Z)
  ref <- hsquared:::hs_reml_estimate_reference(
    payload$y,
    payload$X,
    Z,
    fixture$expected$Ainv,
    method = "REML"
  )

  ll <- function(theta) {
    hsquared:::hs_gaussian_loglik_reference(
      payload$y,
      payload$X,
      Z,
      fixture$expected$Ainv,
      theta[["sigma_a2"]],
      theta[["sigma_e2"]],
      "REML"
    )$loglik
  }
  h2 <- function(theta) {
    theta[["sigma_a2"]] / sum(theta[c("sigma_a2", "sigma_e2")])
  }

  # Core claim: under the same verified REML objective, hsquared's solution is at
  # least as good as the established external package (it reaches the true
  # optimum; pedigreemm's optimizer lands slightly off on pedigree models).
  expect_gte(ll(ref$estimate), ll(ext) - 1e-6)
  # Heritabilities agree within a sane band (flat REML surfaces; not DGP recovery
  # or ASReml parity).
  expect_equal(h2(ref$estimate), h2(ext), tolerance = 0.1)
})

test_that("hsquared's R REML reference recovers the published gryphon estimates", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("enhancer")

  # External published-estimate anchor: the gryphon birth-weight univariate
  # animal model (BWT ~ 1 + animal). This validates hsquared's INDEPENDENT pure-R
  # REML reference optimizer against a named published REML estimate (Wilson et
  # al. 2010), plus an optional agreement check against the external sommer
  # package. It is the first externally-anchored REML-recovery atom in the R
  # lane. It does NOT exercise the production fit path and does NOT satisfy the
  # twin-owned V1-MRODE-FIT gate row; the gryphon population is teaching/simulated
  # data shipped in `enhancer`.
  e <- new.env()
  utils::data("DT_gryphon", package = "enhancer", envir = e)
  DT <- get("DT_gryphon", e)
  A <- get("A_gryphon", e)
  ids <- rownames(A)
  if (is.null(ids)) {
    ids <- as.character(seq_len(nrow(A)))
  }
  dat <- DT[!is.na(DT$BWT), ]

  pub <- hsquared:::hs_gryphon_published_reml()

  # hsquared's independent pure-R REML reference on the same y, X, Z, Ainv.
  y <- dat$BWT
  X <- matrix(1, length(y), 1L)
  j <- match(as.character(dat$ANIMAL), ids)
  Z <- matrix(0, length(y), length(ids))
  Z[cbind(seq_along(j), j)] <- 1
  Ainv <- solve(A)
  ref <- hsquared:::hs_reml_estimate_reference(
    y,
    X,
    Z,
    Ainv,
    method = "REML",
    initial = c(sigma_a2 = 3, sigma_e2 = 4)
  )
  va <- ref$estimate[["sigma_a2"]]
  ve <- ref$estimate[["sigma_e2"]]
  h2 <- va / (va + ve)

  expect_equal(ref$convergence, 0L)
  expect_equal(va, pub[["sigma_a2"]], tolerance = 0.02)
  expect_equal(ve, pub[["sigma_e2"]], tolerance = 0.02)
  expect_equal(h2, pub[["h2"]], tolerance = 0.02)

  # Optional two-sided agreement against the external sommer package, robust to
  # sommer API churn (skip the agreement leg if its API differs, do not fail).
  if (requireNamespace("sommer", quietly = TRUE)) {
    m <- tryCatch(
      {
        d2 <- dat
        d2$ANIMAL <- factor(as.character(d2$ANIMAL))
        sommer::mmes(
          BWT ~ 1,
          random = ~ sommer::vsm(sommer::ism(ANIMAL), Gu = A),
          data = d2,
          verbose = FALSE
        )
      },
      error = function(e) NULL
    )
    if (!is.null(m)) {
      # Maintainer-signed-off V1-COMPARATORS band (2026-06-13): variance
      # components within ~1-2% relative, h2 within ~0.01-0.02 absolute.
      sv <- sort(as.numeric(unlist(m$theta)))
      expect_equal(sv[1], min(va, ve), tolerance = 0.02)
      expect_equal(sv[2], max(va, ve), tolerance = 0.02)
      expect_lt(abs(sv[1] / sum(sv) - va / (va + ve)), 0.02)
    }
  }
})

test_that("REML recovers known variance components from a simulated DGP", {
  # Known-truth recovery (ADEMP; Morris/White/Crowther 2019). Data are simulated
  # from a univariate Gaussian animal model with KNOWN variance components over a
  # clean simulated pedigree; the estimator must recover the generating values
  # (near-unbiased) and produce EBVs that track the true breeding values. This is
  # the statistical-correctness check the v0.1 promotion predicate (item 3) calls
  # for, distinct from optimizer reproducibility. This CI leg uses hsquared's
  # independent pure-R REML reference; the full engine study lives in
  # data-raw/dgp-recovery-study.R. It does NOT flip the twin-owned estimator gate
  # row.
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("nadiv")

  ped <- hs_sim_pedigree(n_founder = 40, n_per_gen = 80, n_gen = 2, seed = 1)
  A <- as.matrix(nadiv::makeA(ped[, c("id", "sire", "dam")]))[ped$id, ped$id]
  U <- chol(A)
  Ainv <- solve(A)
  n <- nrow(ped)
  s2a <- 0.4
  s2e <- 0.6
  set.seed(20240613L)
  seeds <- sample.int(.Machine$integer.max, 25L)

  out <- t(vapply(
    seeds,
    function(sd) {
      set.seed(sd)
      sim <- hs_sim_animal_phenotypes(U, s2a, s2e, mu = 5)
      ref <- hsquared:::hs_reml_estimate_reference(
        sim$y,
        matrix(1, n, 1L),
        diag(n),
        Ainv,
        method = "REML",
        initial = c(sigma_a2 = 0.5, sigma_e2 = 0.5)
      )
      va <- ref$estimate[["sigma_a2"]]
      ve <- ref$estimate[["sigma_e2"]]
      ebv <- hs_sim_blup_ebv(sim$y, Ainv, va, ve)
      c(va, ve, va / (va + ve), ref$convergence, stats::cor(ebv, sim$u))
    },
    numeric(5)
  ))
  colnames(out) <- c("s2a", "s2e", "h2", "conv", "acc")

  # All replicates converged.
  expect_true(all(out[, "conv"] == 0L))
  # Near-unbiased recovery of the known truth (loose absolute band, well above
  # the Monte Carlo SE for 25 replicates).
  expect_lt(abs(mean(out[, "s2a"]) - s2a), 0.06)
  expect_lt(abs(mean(out[, "s2e"]) - s2e), 0.06)
  expect_lt(abs(mean(out[, "h2"]) - 0.4), 0.06)
  # EBVs track the true breeding values (correlation floor for h2 = 0.4).
  expect_gt(mean(out[, "acc"]), 0.5)
})

test_that("the Julia engine recovers the published gryphon estimates via supplied A", {
  # V1-MRODE-FIT engine evidence: both REML optimizers must recover the published
  # gryphon REML estimate (Wilson et al. 2010) when given the published
  # relationship matrix directly. The raw gryphon pedigree is pathological
  # (ancestral loops) and the engine correctly rejects it, so the signed-off
  # V1-MRODE-FIT anchor uses supplied A_gryphon. Tolerance = the signed-off
  # comparator band (~1-2% / h2 ~0.01-0.02).
  hs_skip_live_julia()
  testthat::skip_if_not_installed("enhancer")
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for the live engine gryphon recovery."
  )

  e <- new.env()
  utils::data("DT_gryphon", package = "enhancer", envir = e)
  DT <- get("DT_gryphon", e)
  A <- get("A_gryphon", e)
  ids <- rownames(A)
  dat <- DT[!is.na(DT$BWT), ]
  pub <- hsquared:::hs_gryphon_published_reml()

  y <- dat$BWT
  X <- matrix(1, length(y), 1L)
  j <- match(as.character(dat$ANIMAL), ids)
  Z <- methods::as(
    Matrix::sparseMatrix(
      i = seq_along(j),
      j = j,
      x = 1,
      dims = c(length(y), length(ids))
    ),
    "CsparseMatrix"
  )
  Ainv <- solve(A)

  hs_skip_live_julia()
  hsquared:::hs_julia_setup(hsquared:::hs_default_julia_project())
  JuliaCall::julia_assign("gy", y)
  JuliaCall::julia_assign("gX", X)
  hsquared:::hs_julia_assign_sparse_csc("gZ", Z)
  JuliaCall::julia_assign("gAinv", as.matrix(Ainv))
  JuliaCall::julia_assign("gids", ids)
  JuliaCall::julia_command(paste(
    "gspec = HSquared.animal_model_spec(gy, gX, gZ,",
    "SparseArrays.sparse(gAinv); ids = string.(gids), method = :REML);"
  ))

  for (tg in c("fit_sparse_reml", "fit_ai_reml")) {
    JuliaCall::julia_command(sprintf(
      "gfit = HSquared.%s(gspec; initial = (sigma_a2 = 3.0, sigma_e2 = 4.0), iterations = 2000);",
      tg
    ))
    vc <- JuliaCall::julia_eval(
      "[gfit.variance_components.sigma_a2, gfit.variance_components.sigma_e2]"
    )
    expect_equal(vc[1], pub[["sigma_a2"]], tolerance = 0.02)
    expect_equal(vc[2], pub[["sigma_e2"]], tolerance = 0.02)
    expect_lt(abs(vc[1] / sum(vc) - pub[["h2"]]), 0.02)
  }
})
