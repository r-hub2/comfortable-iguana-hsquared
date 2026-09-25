test_that("multivariate cbind response builds Y payload and preserves NA cells", {
  ped <- data.frame(
    id = c("sire", "dam", "calf1", "calf2"),
    sire = c(NA, NA, "sire", "sire"),
    dam = c(NA, NA, "dam", "dam")
  )
  dat <- data.frame(
    y1 = c(1, 2, 3, 4),
    y2 = c(1.5, NA, 3.5, 4.5),
    sex = c("m", "f", "f", "m"),
    id = ped$id
  )

  spec <- hsquared:::hs_build_model_spec(
    cbind(y1, y2) ~ sex + animal(1 | id, pedigree = ped),
    data = dat,
    family = stats::gaussian(),
    REML = TRUE
  )
  payload <- hsquared:::hs_build_bridge_payload(spec)

  expect_true(spec$response$multivariate)
  expect_equal(spec$response$trait_names, c("y1", "y2"))
  expect_match(spec$bridge$target, "fit_multivariate_reml", fixed = TRUE)
  expect_null(payload$y)
  expect_equal(dim(payload$Y), c(4L, 2L))
  expect_true(is.na(payload$Y[2, 2]))
  expect_equal(payload$metadata$response_type, "multivariate")
  expect_equal(payload$metadata$trait_names, c("y1", "y2"))
  expect_s4_class(payload$Z, "dgCMatrix")
})

test_that("multivariate cbind response requires unique non-empty trait names", {
  ped <- data.frame(
    id = c("sire", "dam", "calf1", "calf2"),
    sire = c(NA, NA, "sire", "sire"),
    dam = c(NA, NA, "dam", "dam")
  )
  dat <- data.frame(
    y1 = c(1, 2, 3, 4),
    y2 = c(1.5, 2.5, 3.5, 4.5),
    id = ped$id
  )

  expect_error(
    hsquared:::hs_build_model_spec(
      cbind(y1, y1) ~ animal(1 | id, pedigree = ped),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "duplicate names: y1",
    fixed = TRUE
  )

  expect_error(
    hsquared:::hs_validate_multivariate_trait_names(c("y1", "")),
    "empty or missing names",
    fixed = TRUE
  )
})

test_that("multivariate parser rejects fixed-effect NA and rank-deficient X", {
  ped <- data.frame(
    id = c("a", "b", "c", "d"),
    sire = c(NA, NA, "a", "a"),
    dam = c(NA, NA, "b", "b")
  )
  dat <- data.frame(
    y1 = c(1, 2, 3, 4),
    y2 = c(1.5, 2.5, 3.5, 4.5),
    x = c(0, 1, 0, 1),
    x_dup = c(0, 1, 0, 1),
    id = ped$id
  )
  dat_na <- dat
  dat_na$x[2] <- NA

  expect_error(
    hsquared:::hs_build_model_spec(
      cbind(y1, y2) ~ x + animal(1 | id, pedigree = ped),
      data = dat_na,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "Missing values in fixed-effect variables",
    fixed = TRUE
  )
  expect_error(
    hsquared:::hs_build_model_spec(
      cbind(y1, y2) ~ x + x_dup + animal(1 | id, pedigree = ped),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "rank deficient",
    fixed = TRUE
  )
})

test_that("multivariate cbind auto-routes on the default path and stays cbind-only", {
  ped <- data.frame(
    id = c("sire", "dam", "calf"),
    sire = c(NA, NA, "sire"),
    dam = c(NA, NA, "dam")
  )
  dat <- data.frame(
    y1 = c(1, 2, 3),
    y2 = c(1.5, 2.5, 3.5),
    id = ped$id
  )

  # MV-4 (doc 38): the multivariate cbind + gaussian response auto-routes to the
  # multivariate REML target on the default path (grammar confirmed Julia-free).
  spec_default <- suppressMessages(hsquared(
    cbind(y1, y2) ~ animal(1 | id, pedigree = ped),
    data = dat,
    control = hs_control(engine = "validate")
  ))
  expect_true(grepl(
    "fit_multivariate_reml",
    spec_default$bridge$target,
    fixed = TRUE
  ))

  # The auto-route (and the engine = "julia" no-target auto-select, §H2) reach the
  # fitter dispatch: without the Julia bridge neither still raises the old
  # "experimental and opt-in" / "requires the opt-in target" aborts.
  if (!hsquared:::hs_julia_bridge_available()) {
    err_default <- tryCatch(
      hsquared(
        cbind(y1, y2) ~ animal(1 | id, pedigree = ped),
        data = dat
      ),
      error = function(e) conditionMessage(e)
    )
    expect_false(grepl("experimental and opt-in", err_default, fixed = TRUE))

    err_julia <- tryCatch(
      hsquared(
        cbind(y1, y2) ~ animal(1 | id, pedigree = ped),
        data = dat,
        control = hs_control(engine = "julia")
      ),
      error = function(e) conditionMessage(e)
    )
    expect_false(grepl("requires the opt-in", err_julia, fixed = TRUE))
  }
  expect_no_warning(suppressMessages(hsquared(
    cbind(y1, y2) ~ animal(1 | id, pedigree = ped),
    data = dat,
    control = hs_control(engine = "validate")
  )))
  expect_error(
    hsquared(
      y1 ~ animal(1 | id, pedigree = ped),
      data = dat,
      control = hs_control(
        engine = "julia",
        engine_control = list(target = "multivariate")
      )
    ),
    "requires a `cbind",
    fixed = TRUE
  )
  expect_equal(
    hsquared:::hs_validate_julia_target("multivariate"),
    "multivariate"
  )
})

test_that("default-path cbind() warns once that multivariate is experimental", {
  hsquared:::hs_reset_session_flags()
  withr::defer(hsquared:::hs_reset_session_flags())
  expect_warning(
    hsquared:::hs_warn_cbind_experimental_once(),
    "This cbind() model fitted; multivariate is covered at validation scale (experimental).",
    fixed = TRUE
  )
  expect_no_warning(hsquared:::hs_warn_cbind_experimental_once())
  hsquared:::hs_reset_session_flags()
  expect_warning(
    hsquared:::hs_warn_cbind_experimental_once(),
    "Report point estimates for t=2 unstructured G0/R0 only; intervals are not coverage-calibrated.",
    fixed = TRUE
  )
})

test_that("multivariate genetic_structure control is fenced", {
  ped <- data.frame(
    id = c("sire", "dam", "calf"),
    sire = c(NA, NA, "sire"),
    dam = c(NA, NA, "dam")
  )
  dat <- data.frame(
    y1 = c(1, 2, 3),
    y2 = c(1.5, 2.5, 3.5),
    id = ped$id
  )

  expect_equal(
    hsquared:::hs_validate_genetic_structure_control(
      hs_control(
        engine = "julia",
        engine_control = list(
          target = "multivariate",
          genetic_structure = "unstructured"
        )
      ),
      "multivariate"
    ),
    "unstructured"
  )
  expect_error(
    hsquared:::hs_validate_genetic_structure_control(
      hs_control(
        engine = "julia",
        engine_control = list(
          target = "multivariate",
          genetic_structure = c("unstructured", "diagonal")
        )
      ),
      "multivariate"
    ),
    "must be a single string",
    fixed = TRUE
  )
  expect_error(
    hsquared:::hs_validate_genetic_structure_control(
      hs_control(
        engine = "julia",
        engine_control = list(
          target = "multivariate",
          genetic_structure = "toeplitz"
        )
      ),
      "multivariate"
    ),
    "must be one of",
    fixed = TRUE
  )
  expect_error(
    hsquared:::hs_validate_genetic_structure_control(
      hs_control(
        engine = "julia",
        engine_control = list(
          target = "ai_reml",
          genetic_structure = "unstructured"
        )
      ),
      "ai_reml"
    ),
    "only planned for the `target = \"multivariate\"` bridge",
    fixed = TRUE
  )
  # "diagonal" is now accepted (no loadings / no rotation ambiguity).
  expect_equal(
    hsquared:::hs_validate_genetic_structure_control(
      hs_control(
        engine = "julia",
        engine_control = list(
          target = "multivariate",
          genetic_structure = "diagonal"
        )
      ),
      "multivariate"
    ),
    "diagonal"
  )
  # "lowrank" / "factor_analytic" stay rejected; FA names the engine-covered
  # twin without becoming an R-public fit.
  for (gs in c("lowrank", "factor_analytic")) {
    expect_error(
      hsquared:::hs_validate_genetic_structure_control(
        hs_control(
          engine = "julia",
          engine_control = list(
            target = "multivariate",
            genetic_structure = gs
          )
        ),
        "multivariate"
      ),
      "not activated on the R bridge",
      fixed = TRUE
    )
  }
  expect_error(
    hsquared:::hs_validate_genetic_structure_control(
      hs_control(
        engine = "julia",
        engine_control = list(
          target = "multivariate",
          genetic_structure = "factor_analytic"
        )
      ),
      "multivariate"
    ),
    "60895208",
    fixed = TRUE
  )
  expect_error(
    hsquared:::hs_validate_genetic_structure_control(
      hs_control(
        engine = "julia",
        engine_control = list(
          target = "multivariate",
          genetic_structure = "factor_analytic"
        )
      ),
      "multivariate"
    ),
    "planned on the R surface",
    fixed = TRUE
  )
  expect_error(
    hsquared:::hs_validate_genetic_structure_control(
      hs_control(
        engine = "julia",
        engine_control = list(
          target = "multivariate",
          genetic_structure = "unstructured",
          rank = 0
        )
      ),
      "multivariate"
    ),
    "must be a single positive integer",
    fixed = TRUE
  )
  expect_error(
    hsquared(
      cbind(y1, y2) ~ animal(1 | id, pedigree = ped),
      data = dat,
      control = hs_control(
        engine = "julia",
        engine_control = list(
          target = "multivariate",
          genetic_structure = "unstructured",
          rank = 1L
        )
      )
    ),
    "reserved for future `lowrank` and `factor_analytic`",
    fixed = TRUE
  )
})

test_that("multivariate initial values are named covariance matrices", {
  expect_equal(
    hsquared:::hs_validate_multivariate_initial(NULL, 2L),
    list(G0 = diag(1, 2L), R0 = diag(1, 2L))
  )
  expect_error(
    hsquared:::hs_validate_multivariate_initial(list(diag(2), diag(2)), 2L),
    "named list",
    fixed = TRUE
  )
  expect_error(
    hsquared:::hs_validate_multivariate_initial(
      list(G0 = diag(2), R0 = matrix(1, 2, 2)),
      2L
    ),
    "positive definite",
    fixed = TRUE
  )
  expect_error(
    hsquared:::hs_fit_julia_multivariate_payload(list()),
    "must be an internal `hs_bridge_payload`",
    fixed = TRUE
  )
})

test_that("multivariate result normalizer exposes G, R, h2, and cross-trait EBVs", {
  payload <- list(
    Y = matrix(c(1, 2, 3, 4, 1.5, NA, 3.5, 4.5), nrow = 4),
    X = matrix(c(1, 1, 1, 1, 0, 1, 0, 1), nrow = 4),
    ids = c("sire", "dam", "calf1", "calf2"),
    family = "gaussian",
    metadata = list(
      fixed_colnames = c("(Intercept)", "sexm"),
      trait_names = c("y1", "y2")
    )
  )
  raw <- list(
    genetic_covariance = matrix(c(1.0, 0.2, 0.2, 1.5), 2),
    residual_covariance = matrix(c(2.0, 0.1, 0.1, 2.5), 2),
    genetic_correlation = matrix(c(1.0, 0.1633, 0.1633, 1.0), 2),
    residual_correlation = matrix(c(1.0, 0.0447, 0.0447, 1.0), 2),
    heritability = c(1 / 3, 1.5 / 4),
    beta = matrix(c(1, 0.5, 2, 0.7), nrow = 2),
    breeding_ids = payload$ids,
    breeding_traits = c("y1", "y2"),
    breeding_values = matrix(seq(0.1, 0.8, length.out = 8), nrow = 4),
    loglik = -22.5,
    converged = TRUE,
    iterations = 18L,
    traits = c("y1", "y2"),
    genetic_structure = "unstructured"
  )

  result <- hsquared:::hs_normalize_multivariate_result(raw, payload)
  fit <- hsquared:::hs_new_fit(
    spec = list(
      method = "REML",
      family = list(family = "gaussian"),
      target = "multivariate"
    ),
    payload = payload,
    result = result
  )

  expect_equal(
    dimnames(genetic_covariance(fit)),
    list(c("y1", "y2"), c("y1", "y2"))
  )
  expect_equal(residual_covariance(fit), result$residual_covariance)
  expect_equal(genetic_correlation(fit), result$genetic_correlation)
  expect_equal(residual_correlation(fit), result$residual_correlation)
  expect_equal(nrow(heritability(fit)), 2L)
  expect_equal(nrow(breeding_values(fit)), 8L)
  expect_equal(stats::nobs(fit), 7L)
  expect_s3_class(stats::logLik(fit), "logLik")
})

test_that("non-converged multivariate fits do not expose logLik or AIC", {
  payload <- list(
    Y = matrix(1:4, nrow = 2),
    X = matrix(1, nrow = 2, ncol = 1),
    ids = c("a", "b"),
    metadata = list(fixed_colnames = "(Intercept)", trait_names = c("y1", "y2"))
  )
  raw <- list(
    genetic_covariance = diag(2),
    residual_covariance = diag(2),
    genetic_correlation = diag(2),
    residual_correlation = diag(2),
    heritability = c(0.5, 0.5),
    beta = matrix(c(1, 2), nrow = 1),
    breeding_ids = c("a", "b"),
    breeding_traits = c("y1", "y2"),
    breeding_values = matrix(0, 2, 2),
    loglik = -1,
    converged = FALSE,
    iterations = 10L,
    traits = c("y1", "y2"),
    genetic_structure = "unstructured"
  )
  fit <- hsquared:::hs_new_fit(
    spec = list(
      method = "REML",
      family = list(family = "gaussian"),
      target = "multivariate"
    ),
    payload = payload,
    result = hsquared:::hs_normalize_multivariate_result(raw, payload)
  )

  expect_error(stats::logLik(fit), "did not converge", fixed = TRUE)
  expect_error(AIC(fit), "did not converge", fixed = TRUE)
})

test_that("R multivariate payload keeps NA; assign copy becomes NaN", {
  ped <- data.frame(
    id = c("sire", "dam", "calf1", "calf2"),
    sire = c(NA, NA, "sire", "sire"),
    dam = c(NA, NA, "dam", "dam")
  )
  dat <- data.frame(
    y1 = c(1, 2, 3, 4),
    y2 = c(1.5, NA, 3.5, 4.5),
    id = ped$id
  )
  spec <- hsquared:::hs_build_model_spec(
    cbind(y1, y2) ~ animal(1 | id, pedigree = ped),
    data = dat,
    family = stats::gaussian(),
    REML = TRUE
  )
  payload <- hsquared:::hs_build_bridge_payload(spec)
  expect_equal(sum(is.na(payload$Y)), 1L)
  expect_true(is.na(payload$Y[2L, 2L]) && !is.nan(payload$Y[2L, 2L]))

  Yj <- hsquared:::hs_y_matrix_for_julia(payload$Y)
  expect_equal(sum(is.nan(Yj)), 1L)
  expect_true(is.nan(Yj[2L, 2L]))
})

test_that("bridge marshals multivariate response NA cells as Julia NaN", {
  # Live: assign the NaN copy, then count via julia_command + Int eval.
  # Do NOT julia_eval(sum(isnan.(...))) on a raw-NA assign: JuliaCall 0.17.6
  # delivers Missing, and that eval path segfaults in Rcpp precious-preserve.
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for live bridge smoke."
  )

  ped <- data.frame(
    id = c("sire", "dam", "calf1", "calf2"),
    sire = c(NA, NA, "sire", "sire"),
    dam = c(NA, NA, "dam", "dam")
  )
  dat <- data.frame(
    y1 = c(1, 2, 3, 4),
    y2 = c(1.5, NA, 3.5, 4.5),
    id = ped$id
  )
  spec <- hsquared:::hs_build_model_spec(
    cbind(y1, y2) ~ animal(1 | id, pedigree = ped),
    data = dat,
    family = stats::gaussian(),
    REML = TRUE
  )
  payload <- hsquared:::hs_build_bridge_payload(spec)
  Yj <- hsquared:::hs_y_matrix_for_julia(payload$Y)

  hsquared:::hs_julia_setup(hsquared:::hs_default_julia_project())
  JuliaCall::julia_assign("hsq_Y_roundtrip", Yj)
  JuliaCall::julia_command(
    "hsq_nan_count = Int(sum(isnan.(hsq_Y_roundtrip)))"
  )
  expect_equal(JuliaCall::julia_eval("hsq_nan_count"), 1L)
})

test_that("hsquared can use the opt-in experimental multivariate REML bridge", {
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for live multivariate bridge smoke."
  )

  ped <- data.frame(
    id = c("s1", "d1", "s2", "d2", "a", "b", "c", "d"),
    sire = c(NA, NA, NA, NA, "s1", "s1", "s2", "s2"),
    dam = c(NA, NA, NA, NA, "d1", "d1", "d2", "d2")
  )
  dat <- data.frame(
    y1 = c(1.0, 1.8, 1.2, 2.0, 3.0, 3.4, 2.8, 3.2),
    y2 = c(2.1, 1.5, 2.2, 1.7, 3.1, NA, 3.0, 2.8),
    id = ped$id
  )

  # PINNED (hsquared#235): this fixture's optimum is flat/boundary, so the
  # engine's observed information is not positive definite and
  # multivariate_covariance_standard_errors() refuses. The bridge surfaces
  # that as one warning rather than a silent NA -- asserting it keeps the
  # surfacing from regressing to silence. One record per animal here, so the
  # repeated-records warning does NOT apply. If the engine later makes these
  # SEs computable this FAILS on purpose: update it deliberately.
  expect_warning(
    fit <- hsquared(
      cbind(y1, y2) ~ animal(1 | id, pedigree = ped),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE,
      control = hs_control(
        engine = "julia",
        engine_control = list(target = "multivariate", iterations = 400L)
      )
    ),
    "multivariate_covariance_standard_errors\\(\\).*observed information is not finite positive-definite"
  )

  expect_s3_class(fit, "hsquared_fit")
  expect_equal(fit$spec$target, "multivariate")
  expect_equal(dim(genetic_covariance(fit)), c(2L, 2L))
  expect_equal(dim(residual_covariance(fit)), c(2L, 2L))
  expect_equal(dim(genetic_correlation(fit)), c(2L, 2L))
  expect_equal(dim(residual_correlation(fit)), c(2L, 2L))
  expect_equal(heritability(fit)$trait, c("y1", "y2"))
  expect_equal(nrow(breeding_values(fit)), 16L)
  expect_equal(stats::nobs(fit), 15L)
  expect_false(is.null(fit$result$genetic_correlation_plot_data))
  expect_false(is.null(fit$result$genetic_pca_plot_data))
  expect_equal(
    fit$result$genetic_correlation_plot_data$traits,
    c("y1", "y2")
  )
  expect_equal(
    fit$result$genetic_correlation_plot_data$genetic_correlations,
    genetic_correlation(fit),
    tolerance = 1e-8
  )
  expect_true(isTRUE(
    fit$result$genetic_pca_plot_data$rotation_invariant
  ))
  expect_true(isTRUE(
    fit$result$genetic_pca_plot_data$is_eigenstructure_not_loadings
  ))
  expect_equal(
    fit_diagnostics(fit)$value[
      fit_diagnostics(fit)$metric == "variance_components_source"
    ],
    "estimated_multivariate_reml"
  )
})

hs_phase4_fixture_path <- function(file) {
  testthat::test_path("fixtures", "phase4_multitrait_parity", file)
}

hs_read_phase4_fixture <- function(file) {
  utils::read.csv(
    hs_phase4_fixture_path(file),
    stringsAsFactors = FALSE,
    na.strings = c("", "NA")
  )
}

hs_phase4_matrix <- function(file) {
  dat <- hs_read_phase4_fixture(file)
  out <- as.matrix(dat[-1])
  storage.mode(out) <- "double"
  rownames(out) <- dat[[1]]
  out
}

test_that("R consumes the shared Phase 4 multivariate parity fixture", {
  ped <- hs_read_phase4_fixture("pedigree.csv")
  names(ped)[names(ped) == "animal"] <- "id"
  ped$sire[ped$sire == "0"] <- NA
  ped$dam[ped$dam == "0"] <- NA

  pheno <- hs_read_phase4_fixture("phenotypes.csv")
  G0 <- hs_phase4_matrix("expected_genetic_covariance.csv")
  R0 <- hs_phase4_matrix("expected_residual_covariance.csv")
  beta <- hs_phase4_matrix("expected_beta.csv")
  h2 <- hs_read_phase4_fixture("expected_heritability.csv")
  ebv <- hs_read_phase4_fixture("expected_ebv.csv")
  metadata <- stats::setNames(
    hs_read_phase4_fixture("expected_metadata.csv")$value,
    hs_read_phase4_fixture("expected_metadata.csv")$key
  )

  spec <- hsquared:::hs_build_model_spec(
    cbind(trait1, trait2) ~ x + animal(1 | animal, pedigree = ped),
    data = pheno,
    family = stats::gaussian(),
    REML = TRUE
  )
  payload <- hsquared:::hs_build_bridge_payload(spec)

  expect_equal(payload$metadata$response_type, "multivariate")
  expect_equal(payload$metadata$trait_names, c("trait1", "trait2"))
  expect_equal(payload$metadata$fixed_colnames, c("(Intercept)", "x"))
  expect_equal(payload$metadata$observed_ids, pheno$animal)
  expect_equal(payload$ids, ped$id)
  expect_equal(dim(payload$Z), c(nrow(pheno), nrow(ped)))
  expect_equal(payload$Y, unname(as.matrix(pheno[c("trait1", "trait2")])))
  expect_equal(payload$X, unname(stats::model.matrix(~x, data = pheno)))

  raw <- list(
    genetic_covariance = G0,
    residual_covariance = R0,
    genetic_correlation = stats::cov2cor(G0),
    residual_correlation = stats::cov2cor(R0),
    heritability = h2$h2,
    beta = beta,
    breeding_ids = ebv$animal,
    breeding_traits = c("trait1", "trait2"),
    breeding_values = as.matrix(ebv[c("trait1", "trait2")]),
    loglik = as.numeric(metadata[["loglik"]]),
    converged = identical(tolower(metadata[["converged"]]), "true"),
    iterations = as.integer(metadata[["iterations"]]),
    traits = c("trait1", "trait2"),
    genetic_structure = "unstructured"
  )
  result <- hsquared:::hs_normalize_multivariate_result(raw, payload)
  fit <- hsquared:::hs_new_fit(
    spec = list(
      method = "REML",
      family = list(family = "gaussian"),
      target = "multivariate"
    ),
    payload = payload,
    result = result
  )

  expect_equal(genetic_covariance(fit), G0, tolerance = 1e-10)
  expect_equal(G_matrix(fit), G0, tolerance = 1e-10)
  expect_equal(residual_covariance(fit), R0, tolerance = 1e-10)
  expect_equal(R_matrix(fit), R0, tolerance = 1e-10)
  expect_equal(genetic_correlation(fit), stats::cov2cor(G0), tolerance = 1e-10)
  expect_equal(residual_correlation(fit), stats::cov2cor(R0), tolerance = 1e-10)
  expect_equal(heritability(fit)$estimate, h2$h2, tolerance = 1e-10)

  # MV-3 covered-flip identity gate (Standard-Tier gate, docs/dev-log/decisions.md;
  # locked citation docs/design/04-validation-canon.md -- Falconer & Mackay 1996;
  # Lynch & Walsh 1998 ch. 4, 21): each derived estimand equals its defining
  # function of the covered components G0, R0. The genetic-correlation identity
  # r_g == cov2cor(G0) is the extractor assertion above; the per-trait
  # heritability identity h2_k == G0[k,k]/(G0[k,k]+R0[k,k]) holds on the engine's
  # serialized values. Both off-diagonals are genuine (not 0/1), so the r_g
  # identity is a real test.
  expect_equal(
    h2$h2,
    unname(diag(G0) / (diag(G0) + diag(R0))),
    tolerance = 1e-5
  )
  expect_gt(abs(stats::cov2cor(G0)[1, 2]), 0.05)
  expect_lt(abs(stats::cov2cor(G0)[1, 2]), 0.99)

  fixed <- fixef(fit)
  expected_fixed <- data.frame(
    term = rep(c("(Intercept)", "x"), times = 2L),
    trait = rep(c("trait1", "trait2"), each = 2L),
    estimate = as.vector(beta),
    stringsAsFactors = FALSE
  )
  expect_equal(fixed, expected_fixed, tolerance = 1e-10)

  expected_ebv <- data.frame(
    id = rep(ebv$animal, times = 2L),
    trait = rep(c("trait1", "trait2"), each = nrow(ebv)),
    value = c(ebv$trait1, ebv$trait2),
    stringsAsFactors = FALSE
  )
  expect_equal(breeding_values(fit), expected_ebv, tolerance = 1e-10)
  expect_equal(stats::nobs(fit), nrow(pheno) * 2L)
  expect_equal(as.numeric(stats::logLik(fit)), as.numeric(metadata[["loglik"]]))
  expect_equal(
    attr(stats::logLik(fit), "df"),
    ncol(payload$X) * 2L + 2L * (2L + 1L)
  )
  expect_equal(
    fit_diagnostics(fit)$value[
      fit_diagnostics(fit)$metric == "genetic_structure"
    ],
    "unstructured"
  )
})

test_that("optional sommer comparator matches the Phase 4 diagonal-residual target", {
  testthat::skip_on_cran()
  # A26b: on NOT_CRAN / maintainer CI, missing Suggests must fail loudly — not
  # skip_if_not_installed silently. CRAN still skips via skip_on_cran above.
  hs_require_suggests("sommer", gate = "MV-1b")
  hs_require_suggests("nadiv", gate = "MV-1b")

  ped <- hs_read_phase4_fixture("pedigree.csv")
  pheno <- hs_read_phase4_fixture("phenotypes.csv")
  G0 <- hs_phase4_matrix("expected_genetic_covariance.csv")
  R0 <- hs_phase4_matrix("expected_residual_covariance.csv")
  h2 <- hs_read_phase4_fixture("expected_heritability.csv")

  ped_a <- ped
  names(ped_a)[names(ped_a) == "animal"] <- "id"
  ped_a$sire[ped_a$sire == "0"] <- NA
  ped_a$dam[ped_a$dam == "0"] <- NA
  A <- suppressWarnings(
    as.matrix(nadiv::makeA(ped_a[, c("id", "sire", "dam")]))
  )
  A <- A[ped$animal, ped$animal]

  long <- stats::reshape(
    pheno[c("record", "animal", "x", "trait1", "trait2")],
    varying = c("trait1", "trait2"),
    v.names = "value",
    timevar = "trait",
    times = c("trait1", "trait2"),
    idvar = "record",
    direction = "long"
  )
  long$trait <- factor(long$trait, levels = c("trait1", "trait2"))
  long$animal <- factor(long$animal, levels = ped$animal)
  long$.record_index <- match(long$record, pheno$record)
  long <- long[with(long, order(trait, .record_index)), ]
  row.names(long) <- NULL

  fit <- tryCatch(
    sommer::mmes(
      value ~ trait + trait:x - 1,
      random = ~ sommer::vsm(
        sommer::usm(trait),
        sommer::ism(animal),
        Gu = A
      ),
      rcov = ~ sommer::vsm(sommer::dsm(trait), sommer::ism(units)),
      data = long,
      verbose = FALSE,
      dateWarning = FALSE,
      nIters = 80L
    ),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    testthat::skip(paste(
      "sommer multivariate comparator API did not fit this fixture:",
      conditionMessage(fit)
    ))
  }
  expect_true(isTRUE(fit$convergence))

  genetic_i <- grep("animal", names(fit$theta), fixed = TRUE)
  residual_i <- match("units", names(fit$theta))
  if (length(genetic_i) != 1L || is.na(residual_i)) {
    testthat::skip(
      "sommer theta layout changed; comparator extraction needs review."
    )
  }

  Ghat <- as.matrix(fit$theta[[genetic_i]])
  Rhat <- as.matrix(fit$theta[[residual_i]])
  dimnames(Ghat) <- dimnames(G0)
  dimnames(Rhat) <- dimnames(R0)

  expect_equal(Ghat, G0, tolerance = 5e-4)
  expect_equal(diag(Rhat), diag(R0), tolerance = 5e-4)
  expect_equal(
    unname(diag(Ghat) / (diag(Ghat) + diag(Rhat))),
    unname(h2$h2),
    tolerance = 5e-4
  )
  expect_equal(Rhat[upper.tri(Rhat)], 0, tolerance = 1e-12)
})

test_that("optional sommer comparator matches the Phase 4 FULL-UNSTRUCTURED target", {
  # MV-1: the diagonal-residual check above (sommer::mmes + dsm(trait)) fixes the
  # residual off-diagonal to zero and so cannot confront the engine's off-diagonal
  # R0[2,1]. The classic sommer::mmer interface fits a FULL UNSTRUCTURED residual
  # (mmes raises an Armadillo out-of-bounds on the unstructured residual in this
  # records-within-animal layout), so this in-suite test reaches the off-diagonal
  # the diagonal check cannot -- promoting the reproducible comparator study
  # (data-raw/multivariate-comparator-study.R) into a CI-gated same-estimand check.
  testthat::skip_on_cran()
  # A26b: headline MV-1 comparator — fail loudly when NOT_CRAN and sommer/nadiv
  # missing; never silent-skip on the maintainer path.
  hs_require_suggests("sommer", gate = "MV-1")
  hs_require_suggests("nadiv", gate = "MV-1")

  ped <- hs_read_phase4_fixture("pedigree.csv")
  pheno <- hs_read_phase4_fixture("phenotypes.csv")
  G0 <- hs_phase4_matrix("expected_genetic_covariance.csv")
  R0 <- hs_phase4_matrix("expected_residual_covariance.csv")
  h2 <- hs_read_phase4_fixture("expected_heritability.csv")
  ebv_target <- hs_read_phase4_fixture("expected_ebv.csv")

  # Rebuild A from the pedigree (nadiv), NOT copied from the engine.
  pedn <- data.frame(
    id = ped$animal,
    sire = ped$sire,
    dam = ped$dam,
    stringsAsFactors = FALSE
  )
  pedn$sire[pedn$sire == "0"] <- NA
  pedn$dam[pedn$dam == "0"] <- NA
  A <- suppressWarnings(as.matrix(nadiv::makeA(pedn)))
  A <- A[ped$animal, ped$animal]
  pheno$animal <- factor(pheno$animal, levels = rownames(A))

  fit <- tryCatch(
    suppressWarnings(sommer::mmer(
      cbind(trait1, trait2) ~ x,
      random = ~ sommer::vsr(animal, Gu = A, Gtc = sommer::unsm(2)),
      rcov = ~ sommer::vsr(units, Gtc = sommer::unsm(2)),
      data = pheno,
      verbose = FALSE,
      dateWarning = FALSE
    )),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    testthat::skip(paste(
      "sommer mmer unstructured comparator did not fit this fixture:",
      conditionMessage(fit)
    ))
  }

  G0_hat <- fit$sigma[["u:animal"]]
  R0_hat <- fit$sigma[["u:units"]]
  if (is.null(G0_hat) || is.null(R0_hat)) {
    testthat::skip(
      "sommer sigma layout changed; comparator extraction needs review."
    )
  }
  dimnames(G0_hat) <- dimnames(G0)
  dimnames(R0_hat) <- dimnames(R0)

  # Full G0 and full R0 -- including the OFF-DIAGONAL residual the diagonal check
  # cannot reach -- match the serialized engine target.
  expect_equal(G0_hat, G0, tolerance = 5e-4)
  expect_equal(R0_hat, R0, tolerance = 5e-4)
  # The residual off-diagonal is genuinely estimated (not fixed to zero) and
  # agrees with the engine within optimiser tolerance.
  expect_gt(abs(R0[2, 1]), 0)
  expect_lt(abs(R0_hat[2, 1] - R0[2, 1]), 5e-4)

  expect_equal(
    unname(diag(G0_hat) / (diag(G0_hat) + diag(R0_hat))),
    unname(h2$h2),
    tolerance = 5e-4
  )

  # EBV confrontation against the serialized target (independent A + optimiser).
  ebv_hat <- fit$U[["u:animal"]]
  ebv_hat_mat <- cbind(
    trait1 = ebv_hat$trait1[ebv_target$animal],
    trait2 = ebv_hat$trait2[ebv_target$animal]
  )
  expect_equal(
    unname(ebv_hat_mat),
    unname(as.matrix(ebv_target[, c("trait1", "trait2")])),
    tolerance = 5e-3
  )
})
