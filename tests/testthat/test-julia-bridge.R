test_that("Julia bridge availability check is conservative", {
  expect_false(hsquared:::hs_julia_bridge_available(project = tempfile()))
})

test_that("Julia result normalizer accepts optional PEV and reliability fields", {
  raw <- list(
    variance_components = list(sigma_a2 = 0.4, sigma_e2 = 0.6),
    heritability = 0.4,
    breeding_values = list(
      ids = c("a", "b"),
      values = c(0.1, -0.1)
    ),
    fixed_effects = c(1.2),
    random_effects = list(
      animal = list(ids = c("a", "b"), values = c(0.1, -0.1))
    ),
    loglik = -12.5,
    df = 3L,
    nobs = 2L,
    predictions = c(1.1, 1.2),
    diagnostics = list(optimizer_status = "test"),
    converged = TRUE,
    prediction_error_variance = list(
      ids = c("a", "b"),
      values = c(0.2, 0.25)
    ),
    reliability = list(ids = c("a", "b"), values = c(0.5, 0.375))
  )
  payload <- list(
    metadata = list(fixed_colnames = "(Intercept)"),
    y = c(1, 2)
  )

  result <- hsquared:::hs_normalize_julia_result(raw, payload)

  expect_equal(
    result$prediction_error_variance,
    data.frame(id = c("a", "b"), value = c(0.2, 0.25))
  )
  expect_equal(
    result$reliability,
    data.frame(id = c("a", "b"), value = c(0.5, 0.375))
  )
})

test_that("hs_data components marshal to Julia HSData [live]", {
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for live HSData marshalling."
  )

  ped <- data.frame(
    id = c("a", "b", "c"),
    sire = c(NA, NA, "a"),
    dam = c(NA, NA, "b")
  )
  dat <- data.frame(
    id = c("a", "b", "c"),
    y = c(1, 2, 3)
  )
  genotypes <- data.frame(
    id = c("a", "c"),
    m1 = c(0, 2),
    m2 = c(1, 1)
  )

  bundle <- hs_data(
    phenotypes = dat,
    pedigree = ped,
    genotypes = genotypes
  )
  status <- data_status(bundle)
  expect_equal(
    status$id_overlap$count[
      status$id_overlap$metric == "phenotypes_without_genotypes"
    ],
    1L
  )

  hsquared:::hs_julia_setup(Sys.getenv("HSQUARED_JULIA_PROJECT"))
  JuliaCall::julia_assign("hsq_hsdata_pheno", bundle$phenotypes)
  JuliaCall::julia_assign("hsq_hsdata_ped", bundle$pedigree)
  JuliaCall::julia_assign("hsq_hsdata_geno", bundle$genotypes)
  JuliaCall::julia_command(paste(
    "hsq_hsdata = HSquared.HSData(",
    "hsq_hsdata_pheno;",
    "id = :id,",
    "pedigree = hsq_hsdata_ped,",
    "pedigree_id = :id,",
    "genotypes = hsq_hsdata_geno,",
    "genotype_id = :id",
    ");"
  ))

  expect_equal(
    JuliaCall::julia_eval(
      "collect(String, string.(HSquared.id_map(hsq_hsdata).phenotype_ids))"
    ),
    c("a", "b", "c")
  )
  expect_equal(
    JuliaCall::julia_eval(
      "collect(String, string.(HSquared.id_map(hsq_hsdata).genotype_ids))"
    ),
    c("a", "c")
  )
  expect_equal(
    JuliaCall::julia_eval(
      "collect(String, string.(HSquared.id_map(hsq_hsdata).phenotypes_without_genotypes))"
    ),
    "b"
  )
  expect_equal(
    JuliaCall::julia_eval(
      "collect(String, string.(HSquared.data_status(hsq_hsdata).components))"
    ),
    c("phenotypes", "pedigree", "genotypes")
  )
})

test_that("experimental Julia bridge smoke fits the tiny payload", {
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for live bridge smoke."
  )

  ped <- data.frame(
    id = c("sire", "dam", "calf"),
    sire = c(NA, NA, "sire"),
    dam = c(NA, NA, "dam")
  )
  dat <- data.frame(
    y = c(1, 2.5, 4),
    id = c("sire", "dam", "calf")
  )
  spec <- hsquared:::hs_build_model_spec(
    y ~ animal(1 | id, pedigree = ped),
    data = dat,
    family = stats::gaussian(),
    REML = FALSE
  )
  payload <- hsquared:::hs_build_bridge_payload(spec)

  fit <- hsquared:::hs_fit_julia_payload(
    payload,
    initial = c(sigma_a2 = 0.8, sigma_e2 = 0.4)
  )

  expect_s3_class(fit, "hsquared_fit")
  expect_true(fit$result$converged)
  expect_equal(fit$result$nobs, 3L)
  expect_true(is.finite(fit$result$loglik))
  expect_equal(
    fit$result$variance_components$component,
    c("animal", "residual")
  )
  expect_true(all(fit$result$variance_components$estimate >= 0))
  expect_equal(breeding_values(fit)$id, c("sire", "dam", "calf"))
  expect_equal(names(fixef(fit)), "(Intercept)")
  expect_true(is.finite(heritability(fit)$estimate))
  expect_equal(prediction_error_variance(fit)$id, c("sire", "dam", "calf"))
  expect_equal(reliability(fit)$id, c("sire", "dam", "calf"))
  expect_true(all(is.finite(prediction_error_variance(fit)$value)))
  expect_true(all(is.finite(reliability(fit)$value)))
  expect_s3_class(stats::logLik(fit), "logLik")
})

test_that("hsquared can use the opt-in experimental Julia engine", {
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for live bridge smoke."
  )

  ped <- data.frame(
    id = c("sire", "dam", "calf"),
    sire = c(NA, NA, "sire"),
    dam = c(NA, NA, "dam")
  )
  dat <- data.frame(
    y = c(1, 2.5, 4),
    id = c("sire", "dam", "calf")
  )

  fit <- hsquared(
    y ~ animal(1 | id, pedigree = ped),
    data = dat,
    family = stats::gaussian(),
    REML = TRUE,
    control = hs_control(
      engine = "julia",
      engine_control = list(
        initial = c(sigma_a2 = 0.8, sigma_e2 = 0.4)
      )
    )
  )

  expect_s3_class(fit, "hsquared_fit")
  expect_equal(fit$engine, "HSquared.jl")
  expect_equal(fit$result$nobs, 3L)
  expect_equal(breeding_values(fit)$id, c("sire", "dam", "calf"))
  expect_true(is.finite(heritability(fit)$estimate))
  expect_equal(prediction_error_variance(fit)$id, c("sire", "dam", "calf"))
  expect_equal(reliability(fit)$id, c("sire", "dam", "calf"))
})

test_that("metafounder target is supplied-variance and formula-gated", {
  ped <- data.frame(
    id = c("sire", "dam", "calf"),
    sire = c(NA, NA, "sire"),
    dam = c(NA, NA, "dam")
  )
  dat <- data.frame(y = c(1, 2.5, 4), id = c("sire", "dam", "calf"))
  mf_group <- c(sire = "base", dam = "base", calf = "")
  Gamma <- matrix(0.25, nrow = 1, dimnames = list("base", "base"))

  expect_error(
    hsquared(
      y ~ metafounder(1 | id, pedigree = ped, group = mf_group, Gamma = Gamma),
      data = dat,
      control = hs_control(
        engine = "julia",
        engine_control = list(target = "metafounder")
      )
    ),
    "`engine_control$variance_components` is required when `target = \"metafounder\"`.",
    fixed = TRUE
  )

  expect_error(
    hsquared(
      y ~ animal(1 | id, pedigree = ped),
      data = dat,
      control = hs_control(
        engine = "julia",
        engine_control = list(
          target = "metafounder",
          variance_components = c(sigma_a2 = 0.8, sigma_e2 = 0.4)
        )
      )
    ),
    "`target = \"metafounder\"` requires",
    fixed = TRUE
  )
})

test_that("metafounder supplied-variance bridge matches Gamma-zero Henderson path [live]", {
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for live metafounder bridge validation."
  )

  ped <- data.frame(
    id = c("sire", "dam", "calf", "sib"),
    sire = c(NA, NA, "sire", "sire"),
    dam = c(NA, NA, "dam", "dam")
  )
  dat <- data.frame(
    y = c(1, 2.5, 4, 3.2),
    id = c("sire", "dam", "calf", "sib")
  )
  mf_group <- c(sire = "base", dam = "base", calf = "", sib = "")
  Gamma0 <- matrix(0, nrow = 1, dimnames = list("base", "base"))
  Gammag <- matrix(0.25, nrow = 1, dimnames = list("base", "base"))
  vc <- c(sigma_a2 = 0.8, sigma_e2 = 0.4)

  ordinary <- hsquared(
    y ~ animal(1 | id, pedigree = ped),
    data = dat,
    control = hs_control(
      engine = "julia",
      engine_control = list(
        target = "henderson_mme",
        variance_components = vc
      )
    )
  )
  mf0 <- hsquared(
    y ~ metafounder(1 | id, pedigree = ped, group = mf_group, Gamma = Gamma0),
    data = dat,
    control = hs_control(
      engine = "julia",
      engine_control = list(
        target = "metafounder",
        variance_components = vc
      )
    )
  )
  mfg <- hsquared(
    y ~ metafounder(1 | id, pedigree = ped, group = mf_group, Gamma = Gammag),
    data = dat,
    control = hs_control(
      engine = "julia",
      engine_control = list(
        target = "metafounder",
        variance_components = vc
      )
    )
  )

  expect_equal(mf0$spec$target, "metafounder")
  expect_equal(
    fit_diagnostics(mf0)$value[
      fit_diagnostics(mf0)$metric == "variance_components_source"
    ],
    "supplied_metafounder"
  )
  expect_equal(
    fit_diagnostics(mf0)$value[
      fit_diagnostics(mf0)$metric == "gamma_source"
    ],
    "supplied"
  )
  expect_equal(variance_components(mf0)$component, c("metafounder", "residual"))
  expect_equal(heritability(mf0)$term, "metafounder")
  expect_equal(gamma_matrix(mf0), Gamma0)
  expect_equal(gamma_matrix(mfg), Gammag)
  expect_equal(
    metafounder_groups(mfg),
    data.frame(
      id = c("sire", "dam", "calf", "sib"),
      metafounder_group = c("base", "base", NA, NA),
      is_metafounder = c(TRUE, TRUE, FALSE, FALSE),
      stringsAsFactors = FALSE
    )
  )
  expect_equal(
    breeding_values(mf0)$value,
    breeding_values(ordinary)$value,
    tolerance = 1e-8
  )
  expect_equal(prediction_error_variance(mf0)$id, ped$id)
  expect_equal(reliability(mf0)$id, ped$id)

  bv0 <- breeding_values(mf0)
  bvg <- breeding_values(mfg)
  expect_equal(bvg$id, bv0$id)
  expect_true(all(is.finite(bvg$value)))
  expect_gt(mean(abs(bvg$value - bv0$value)), 1e-8)
})

test_that("the opt-in engine = \"julia\" estimation path rejects REML = FALSE (ML)", {
  ped <- data.frame(
    id = c("sire", "dam", "calf"),
    sire = c(NA, NA, "sire"),
    dam = c(NA, NA, "dam")
  )
  dat <- data.frame(y = c(1, 2.5, 4), id = c("sire", "dam", "calf"))

  # ML is not implemented in v0.1. The estimation targets (the default
  # fit_animal_model, and the REML-only sparse_reml/ai_reml) must reject
  # REML = FALSE rather than run the ML optimizer or silently ignore the
  # request. This is a pure request-validity error: it fires before any Julia
  # engine call, so it needs no live bridge.
  for (tgt in list(NULL, "ai_reml", "sparse_reml")) {
    ec <- if (is.null(tgt)) list() else list(target = tgt)
    expect_error(
      hsquared(
        y ~ animal(1 | id, pedigree = ped),
        data = dat,
        family = stats::gaussian(),
        REML = FALSE,
        control = hs_control(engine = "julia", engine_control = ec)
      ),
      "ML estimation",
      fixed = TRUE
    )
  }
})

test_that("sparse CSC slot helper exposes R Matrix slots", {
  Z <- Matrix::sparseMatrix(
    i = c(1, 3, 2),
    j = c(1, 2, 3),
    x = c(1, 2, 3),
    dims = c(3, 3)
  )

  slots <- hsquared:::hs_sparse_csc_slots(Z)

  expect_equal(slots$nrow, 3L)
  expect_equal(slots$ncol, 3L)
  expect_equal(slots$colptr, Z@p)
  expect_equal(slots$rowval, Z@i)
  expect_equal(slots$nzval, Z@x)

  expect_error(
    hsquared:::hs_sparse_csc_slots(matrix(1)),
    "`x` must be a `Matrix::dgCMatrix` object.",
    fixed = TRUE
  )
})

test_that("experimental Julia bridge uses sparse Z marshalling", {
  ped <- data.frame(
    id = c("sire", "dam", "calf"),
    sire = c(NA, NA, "sire"),
    dam = c(NA, NA, "dam")
  )
  dat <- data.frame(y = c(1, 2.5, 4), id = c("sire", "dam", "calf"))
  spec <- hsquared:::hs_build_model_spec(
    y ~ animal(1 | id, pedigree = ped),
    data = dat,
    family = stats::gaussian(),
    REML = FALSE
  )
  payload <- hsquared:::hs_build_bridge_payload(spec)

  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for live bridge smoke."
  )

  fit <- hsquared:::hs_fit_julia_payload(
    payload,
    initial = c(sigma_a2 = 0.8, sigma_e2 = 0.4)
  )

  expect_s3_class(fit, "hsquared_fit")
  expect_equal(fit$result$nobs, 3L)
  expect_equal(prediction_error_variance(fit)$id, c("sire", "dam", "calf"))
  expect_equal(reliability(fit)$id, c("sire", "dam", "calf"))
})

test_that("Julia bridge validates target control", {
  expect_equal(
    hsquared:::hs_validate_julia_target("fit_animal_model"),
    "fit_animal_model"
  )
  expect_equal(
    hsquared:::hs_validate_julia_target("henderson_mme"),
    "henderson_mme"
  )
  expect_equal(
    hsquared:::hs_validate_julia_target("sparse_reml"),
    "sparse_reml"
  )
  expect_error(
    hsquared:::hs_validate_julia_target("AI_REML"),
    "`engine_control\\$target` must be one of",
    fixed = FALSE
  )
})

test_that("Julia bridge validates iterations control", {
  expect_equal(hsquared:::hs_validate_iterations(500), 500L)
  expect_equal(hsquared:::hs_validate_iterations(1000L), 1000L)
  expect_error(
    hsquared:::hs_validate_iterations(0),
    "must be a single positive integer",
    fixed = TRUE
  )
  expect_error(
    hsquared:::hs_validate_iterations(c(10, 20)),
    "must be a single positive integer",
    fixed = TRUE
  )
})

test_that("Julia bridge validates em_warmup control", {
  # 0 (off / default) is valid, unlike iterations which must be positive.
  expect_equal(hsquared:::hs_validate_em_warmup(0), 0L)
  expect_equal(hsquared:::hs_validate_em_warmup(3), 3L)
  expect_equal(hsquared:::hs_validate_em_warmup(5L), 5L)
  expect_error(
    hsquared:::hs_validate_em_warmup(-1),
    "must be a single non-negative integer",
    fixed = TRUE
  )
  expect_error(
    hsquared:::hs_validate_em_warmup(c(1, 2)),
    "must be a single non-negative integer",
    fixed = TRUE
  )
})

test_that("sparse REML payload requires an internal bridge payload", {
  expect_error(
    hsquared:::hs_fit_julia_sparse_reml_payload(list()),
    "must be an internal `hs_bridge_payload`",
    fixed = TRUE
  )
})

test_that("Julia bridge accepts the ai_reml target", {
  expect_equal(
    hsquared:::hs_validate_julia_target("ai_reml"),
    "ai_reml"
  )
})

test_that("AI-REML payload requires an internal bridge payload", {
  expect_error(
    hsquared:::hs_fit_julia_ai_reml_payload(list()),
    "must be an internal `hs_bridge_payload`",
    fixed = TRUE
  )
})

test_that("Julia Henderson MME bridge requires supplied variance components", {
  expect_error(
    hsquared:::hs_validate_supplied_variances(NULL),
    "`engine_control$variance_components` is required",
    fixed = TRUE
  )
  expect_error(
    hsquared:::hs_validate_supplied_variances(c(sigma_a2 = 1)),
    "must include `sigma_a2` and `sigma_e2`",
    fixed = TRUE
  )
  expect_error(
    hsquared:::hs_validate_supplied_variances(c(sigma_a2 = 1, sigma_e2 = 0)),
    "must be positive and finite",
    fixed = TRUE
  )
  expect_equal(
    hsquared:::hs_validate_supplied_variances(c(
      sigma_a2 = 1.2,
      sigma_e2 = 0.8
    )),
    c(sigma_a2 = 1.2, sigma_e2 = 0.8)
  )
})

test_that("Julia Henderson MME normalizer accepts optional PEV and reliability fields", {
  raw <- list(
    fixed_effects = c(3.2),
    animal_ids = c("sire", "dam", "calf"),
    animal_effects = c(-0.1, 0.2, 0.3),
    fitted = c(3.1, 3.4, 3.5),
    nobs = 3L,
    prediction_error_variance = list(
      ids = c("sire", "dam", "calf"),
      values = c(0.2, 0.25, 0.3)
    ),
    reliability = list(
      ids = c("sire", "dam", "calf"),
      values = c(0.8, 0.75, 0.7)
    )
  )
  payload <- list(metadata = list(fixed_colnames = "(Intercept)"))

  result <- hsquared:::hs_normalize_julia_henderson_mme_result(
    raw,
    payload,
    variance_components = c(sigma_a2 = 1.2, sigma_e2 = 0.8)
  )

  expect_equal(
    result$prediction_error_variance,
    data.frame(id = c("sire", "dam", "calf"), value = c(0.2, 0.25, 0.3))
  )
  expect_equal(
    result$reliability,
    data.frame(id = c("sire", "dam", "calf"), value = c(0.8, 0.75, 0.7))
  )
})

test_that("experimental Julia Henderson MME bridge matches validation fixture", {
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for live Henderson MME bridge validation."
  )

  fixture <- hsquared:::hs_henderson_mme_validation_fixture()
  spec <- hsquared:::hs_build_model_spec(
    fixture$formula,
    data = fixture$data,
    family = stats::gaussian(),
    REML = FALSE
  )
  payload <- hsquared:::hs_build_bridge_payload(spec)
  fit <- hsquared:::hs_fit_julia_henderson_mme_payload(
    payload,
    variance_components = c(
      sigma_a2 = fixture$sigma_a2,
      sigma_e2 = fixture$sigma_e2
    )
  )

  expect_s3_class(fit, "hsquared_fit")
  expect_true(fit$result$converged)
  expect_equal(fit$result$diagnostics$target, "henderson_mme")
  expect_equal(
    variance_components(fit)$estimate,
    c(fixture$sigma_a2, fixture$sigma_e2),
    tolerance = 1e-12
  )
  expect_equal(
    fixef(fit),
    fixture$expected$fixed_effects,
    tolerance = 1e-10,
    ignore_attr = TRUE
  )
  expect_equal(
    breeding_values(fit),
    fixture$expected$breeding_values,
    tolerance = 1e-10
  )
  expect_equal(
    predict(fit),
    fixture$expected$fitted,
    tolerance = 1e-10
  )
  expect_equal(heritability(fit)$estimate, fixture$expected$heritability)
  # PEV/reliability are now attached unconditionally on the Henderson MME path
  # (dense, validation-scale): assert they are present, not merely "if present".
  expect_false(is.null(fit$result$prediction_error_variance))
  expect_false(is.null(fit$result$reliability))
  expect_equal(
    prediction_error_variance(fit)$id,
    fixture$expected$breeding_values$id
  )
  expect_equal(reliability(fit)$id, fixture$expected$breeding_values$id)
  expect_true(all(is.finite(prediction_error_variance(fit)$value)))
  expect_true(all(is.finite(reliability(fit)$value)))
  expect_error(
    stats::logLik(fit),
    "does not contain log-likelihood",
    fixed = TRUE
  )
})

test_that("hsquared can use the opt-in supplied-variance Henderson MME bridge", {
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for live Henderson MME bridge validation."
  )

  fixture <- hsquared:::hs_henderson_mme_validation_fixture()
  fit <- hsquared(
    fixture$formula,
    data = fixture$data,
    family = stats::gaussian(),
    REML = FALSE,
    control = hs_control(
      engine = "julia",
      engine_control = list(
        target = "henderson_mme",
        variance_components = c(
          sigma_a2 = fixture$sigma_a2,
          sigma_e2 = fixture$sigma_e2
        )
      )
    )
  )

  expect_s3_class(fit, "hsquared_fit")
  expect_equal(fit$spec$target, "henderson_mme")
  expect_equal(breeding_values(fit), fixture$expected$breeding_values)
})

test_that("hsquared can use the opt-in experimental sparse REML estimator bridge", {
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for live sparse REML estimator bridge validation."
  )

  fixture <- hsquared:::hs_mrode_supplied_variance_validation_fixture()
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
        iterations = 500L
      )
    )
  )

  # Experimental, opt-in, Julia-owned REML optimizer that R only surfaces.
  # Honest behaviour checks only: not DGP recovery, not ASReml parity.
  expect_s3_class(fit, "hsquared_fit")
  expect_equal(fit$engine, "HSquared.jl")
  expect_equal(fit$spec$target, "sparse_reml")
  est <- variance_components(fit)$estimate
  expect_equal(variance_components(fit)$component, c("animal", "residual"))
  expect_true(all(is.finite(est)) && all(est > 0))
  expect_true(is.finite(fit$result$loglik))
  h2 <- heritability(fit)$estimate
  expect_true(is.finite(h2) && h2 > 0 && h2 < 1)

  # B3: estimated-vs-supplied provenance is explicit and distinct.
  diag <- fit_diagnostics(fit)
  expect_equal(diag$value[diag$metric == "target"], "sparse_reml")
  expect_equal(
    diag$value[diag$metric == "variance_components_source"],
    "estimated_sparse_reml"
  )
})

test_that("hsquared can use the opt-in experimental AI-REML estimator bridge", {
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for live AI-REML estimator bridge validation."
  )

  fixture <- hsquared:::hs_mrode_supplied_variance_validation_fixture()
  fit <- hsquared(
    fixture$formula,
    data = fixture$data,
    family = stats::gaussian(),
    REML = TRUE,
    control = hs_control(
      engine = "julia",
      engine_control = list(
        target = "ai_reml",
        initial = c(sigma_a2 = 1, sigma_e2 = 1),
        iterations = 200L
      )
    )
  )

  # Experimental, opt-in, Julia-owned average-information REML optimizer that R
  # only surfaces. Honest behaviour checks only: not DGP recovery, not ASReml
  # parity.
  expect_s3_class(fit, "hsquared_fit")
  expect_equal(fit$engine, "HSquared.jl")
  expect_equal(fit$spec$target, "ai_reml")
  est <- variance_components(fit)$estimate
  expect_equal(variance_components(fit)$component, c("animal", "residual"))
  expect_true(all(is.finite(est)) && all(est > 0))
  expect_true(is.finite(fit$result$loglik))
  h2 <- heritability(fit)$estimate
  expect_true(is.finite(h2) && h2 > 0 && h2 < 1)

  # Estimated-vs-supplied provenance is explicit and distinct from sparse_reml.
  diag <- fit_diagnostics(fit)
  expect_equal(diag$value[diag$metric == "target"], "ai_reml")
  expect_equal(
    diag$value[diag$metric == "variance_components_source"],
    "estimated_ai_reml"
  )
})

test_that("em_warmup forwards through the bridge and is optimum-invariant", {
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for live em_warmup bridge validation."
  )

  fixture <- hsquared:::hs_mrode_supplied_variance_validation_fixture()
  fit_one <- function(em) {
    hsquared(
      fixture$formula,
      data = fixture$data,
      family = stats::gaussian(),
      REML = TRUE,
      control = hs_control(
        engine = "julia",
        engine_control = list(
          target = "ai_reml",
          initial = c(sigma_a2 = 1, sigma_e2 = 1),
          iterations = 200L,
          em_warmup = em
        )
      )
    )
  }

  # Forwarding: a non-zero `em_warmup` must reach the engine as a well-formed
  # `fit_ai_reml(...; em_warmup = k)` call. A malformed/mis-named kwarg would make
  # the Julia call error, so a successful fit here is the forwarding check.
  warm_fit <- fit_one(3L)
  expect_s3_class(warm_fit, "hsquared_fit")
  expect_equal(warm_fit$spec$target, "ai_reml")

  # Optimum-invariance (matches the engine A finding, V1-AI-REML): the EM warm-start
  # only changes the AI iterations' starting point, so on this identified fixture
  # it reaches the SAME REML optimum as em_warmup = 0. NOTE: the bad-start
  # convergence *rescue* is an engine property validated in HSquared.jl, not here.
  base <- variance_components(fit_one(0L))$estimate
  warm <- variance_components(warm_fit)$estimate
  expect_equal(warm, base, tolerance = 1e-6)
})

test_that("the default engine fits the v0.1 contract via ai_reml", {
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for the default-fit path."
  )

  fixture <- hsquared:::hs_mrode_supplied_variance_validation_fixture()
  # Default control => engine = "fit": the default call now fits the v0.1
  # Gaussian animal model by REML (average-information) through the engine.
  fit <- hsquared(
    fixture$formula,
    data = fixture$data,
    family = stats::gaussian(),
    REML = TRUE
  )

  expect_s3_class(fit, "hsquared_fit")
  expect_equal(fit$engine, "HSquared.jl")
  expect_equal(fit$spec$target, "ai_reml")
  est <- variance_components(fit)$estimate
  expect_true(all(is.finite(est)) && all(est > 0))
  h2 <- heritability(fit)$estimate
  expect_true(is.finite(h2) && h2 > 0 && h2 < 1)
  diag <- fit_diagnostics(fit)
  expect_equal(
    diag$value[diag$metric == "variance_components_source"],
    "estimated_ai_reml"
  )
})
