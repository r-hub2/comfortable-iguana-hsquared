test_that("internal hsquared_fit object supports v0.1 extractors", {
  result <- list(
    variance_components = data.frame(
      component = c("animal", "residual"),
      estimate = c(0.4, 0.6)
    ),
    heritability = data.frame(term = "animal", estimate = 0.4),
    genetic_covariance = matrix(0.4, 1, 1),
    residual_covariance = matrix(0.6, 1, 1),
    genetic_correlation = matrix(1, 1, 1),
    residual_correlation = matrix(1, 1, 1),
    breeding_values = data.frame(id = c("a", "b"), value = c(0.1, -0.1)),
    prediction_error_variance = data.frame(
      id = c("a", "b"),
      value = c(0.2, 0.25)
    ),
    reliability = data.frame(id = c("a", "b"), value = c(0.8, 0.75)),
    marker_effects = data.frame(marker = c("m1", "m2"), effect = c(0.2, -0.1)),
    marker_variance_explained = data.frame(
      marker = c("m1", "m2"),
      proportion = c(0.05, 0.02)
    ),
    qtl_table = data.frame(marker = "m1", lod = 3.2),
    gwas_table = data.frame(marker = "m1", p_value = 0.01),
    eqtl_table = data.frame(gene = "g1", marker = "m1", p_value = 0.02),
    lod_scores = data.frame(position = c(10, 20), lod = c(2.1, 3.2)),
    fixed_effects = c("(Intercept)" = 1.2, sexm = -0.3),
    random_effects = list(animal = c(a = 0.1, b = -0.1)),
    loglik = -12.5,
    df = 4L,
    nobs = 10L,
    predictions = data.frame(.fitted = seq_len(10) + 0.5),
    diagnostics = list(
      optimizer_status = "converged",
      iterations = 7L,
      gradient_norm = 0.001
    ),
    converged = TRUE
  )
  Gamma <- matrix(
    c(0.25, 0.05, 0.05, 0.2),
    nrow = 2,
    dimnames = list(c("base1", "base2"), c("base1", "base2"))
  )

  fit <- hsquared:::hs_new_fit(
    call = quote(hsquared(y ~ animal(1 | id, pedigree = ped), data = dat)),
    spec = list(method = "REML", family = list(family = "gaussian")),
    payload = list(
      y = seq_len(10),
      ids = c("a", "b"),
      group_of = c(a = "base1", b = ""),
      Gamma = Gamma
    ),
    result = result
  )

  expect_s3_class(fit, "hsquared_fit")
  expect_equal(variance_components(fit), result$variance_components)
  expect_equal(heritability(fit), result$heritability)
  expect_equal(genetic_covariance(fit), result$genetic_covariance)
  expect_equal(G_matrix(fit), result$genetic_covariance)
  expect_equal(residual_covariance(fit), result$residual_covariance)
  expect_equal(R_matrix(fit), result$residual_covariance)
  expect_equal(genetic_correlation(fit), result$genetic_correlation)
  expect_equal(residual_correlation(fit), result$residual_correlation)
  expect_equal(gamma_matrix(fit), Gamma)
  expect_equal(
    metafounder_groups(fit),
    data.frame(
      id = c("a", "b"),
      metafounder_group = c("base1", NA),
      is_metafounder = c(TRUE, FALSE),
      stringsAsFactors = FALSE
    )
  )
  expect_equal(breeding_values(fit), result$breeding_values)
  expect_equal(EBV(fit), result$breeding_values)
  expect_equal(BLUP(fit), result$breeding_values)
  expect_equal(prediction_error_variance(fit), result$prediction_error_variance)
  expect_equal(reliability(fit), result$reliability)
  expect_equal(
    accuracy(fit),
    transform(result$reliability, value = sqrt(value))
  )
  expect_equal(marker_effects(fit), result$marker_effects)
  expect_equal(
    marker_variance_explained(fit),
    result$marker_variance_explained
  )
  expect_equal(qtl_table(fit), result$qtl_table)
  expect_equal(gwas_table(fit), result$gwas_table)
  expect_equal(eqtl_table(fit), result$eqtl_table)
  expect_equal(lod_scores(fit), result$lod_scores)
  expect_equal(stats::coef(fit), result$fixed_effects)
  expect_equal(fixef(fit), result$fixed_effects)
  expect_equal(ranef(fit), result$random_effects)
  expect_equal(predict(fit), result$predictions)
  expect_equal(fitted(fit), result$predictions$.fitted)
  expect_equal(residuals(fit), seq_len(10) - result$predictions$.fitted)
  expect_equal(stats::nobs(fit), 10L)
  expect_equal(as.numeric(logLik(fit)), -12.5)
  expect_equal(attr(logLik(fit), "df"), 4L)
  expect_equal(attr(logLik(fit), "nobs"), 10L)
  expect_equal(AIC(fit), 2 * 4 - 2 * -12.5)
  diagnostics <- fit_diagnostics(fit)
  expect_s3_class(diagnostics, "hs_fit_diagnostics")
  expect_equal(
    diagnostics$metric,
    c(
      "engine",
      "method",
      "family",
      "target",
      "converged",
      "optimizer_status",
      "iterations",
      "loglik",
      "df",
      "nobs",
      "at_boundary",
      "gradient_norm"
    )
  )
  expect_equal(
    diagnostics$value[diagnostics$metric == "gradient_norm"],
    "0.001"
  )
  # Interior fit (h2 = 0.4) is not flagged at a variance-component boundary.
  expect_equal(
    diagnostics$value[diagnostics$metric == "at_boundary"],
    "FALSE"
  )
  expect_match(capture.output(print(diagnostics))[[1L]], "<hs_fit_diagnostics>")
  expect_s3_class(summary(fit), "summary_hsquared_fit")

  mf_result <- result
  mf_result$variance_components$component <- c("metafounder", "residual")
  mf_fit <- hsquared:::hs_new_fit(
    spec = list(method = "REML", family = list(family = "gaussian")),
    payload = list(
      y = seq_len(10),
      ids = c("a", "b"),
      group_of = c("base", "")
    ),
    result = mf_result
  )
  mf_diag <- fit_diagnostics(mf_fit)
  expect_equal(
    mf_diag$value[mf_diag$metric == "at_boundary"],
    "FALSE"
  )
})

test_that("extractor defaults do not imply fitted model support", {
  expect_error(
    variance_components(list()),
    "requires an `hsquared_fit` object",
    fixed = TRUE
  )
  expect_error(
    heritability(list()),
    "requires an `hsquared_fit` object",
    fixed = TRUE
  )
  expect_error(
    genetic_covariance(list()),
    "requires an `hsquared_fit` object",
    fixed = TRUE
  )
  expect_error(
    G_matrix(list()),
    "requires an `hsquared_fit` object",
    fixed = TRUE
  )
  expect_error(
    R_matrix(list()),
    "requires an `hsquared_fit` object",
    fixed = TRUE
  )
  expect_error(
    residual_correlation(list()),
    "requires an `hsquared_fit` object",
    fixed = TRUE
  )
  expect_error(
    gamma_matrix(list()),
    "requires an `hsquared_fit` object",
    fixed = TRUE
  )
  expect_error(
    metafounder_groups(list()),
    "requires an `hsquared_fit` object",
    fixed = TRUE
  )
  expect_error(
    metafounder_effects(list()),
    "requires an `hsquared_fit` object",
    fixed = TRUE
  )
  expect_error(
    breeding_values(list()),
    "requires an `hsquared_fit` object",
    fixed = TRUE
  )
  expect_error(
    EBV(list()),
    "requires an `hsquared_fit` object",
    fixed = TRUE
  )
  expect_error(
    BLUP(list()),
    "requires an `hsquared_fit` object",
    fixed = TRUE
  )
  expect_error(
    prediction_error_variance(list()),
    "requires an `hsquared_fit` object",
    fixed = TRUE
  )
  expect_error(
    reliability(list()),
    "requires an `hsquared_fit` object",
    fixed = TRUE
  )
  expect_error(
    accuracy(list()),
    "requires an `hsquared_fit` object",
    fixed = TRUE
  )
  expect_error(
    fit_diagnostics(list()),
    "requires an `hsquared_fit` object",
    fixed = TRUE
  )
  expect_error(
    marker_effects(list()),
    "requires an `hsquared_fit` object with marker",
    fixed = TRUE
  )
  expect_error(
    marker_variance_explained(list()),
    "populated only ",
    fixed = TRUE
  )
  expect_error(
    qtl_table(list()),
    "reserved for the planned map-annotated scan API",
    fixed = TRUE
  )
  # The reservation now points users to the post-fit gwas() verb and names the
  # uncalibrated-significance honesty gate, instead of staying silent.
  expect_error(
    gwas_table(list()),
    "call `gwas_table(scan)`",
    fixed = TRUE
  )
  expect_error(
    eqtl_table(list()),
    "future calibrated-threshold activation",
    fixed = TRUE
  )
  expect_error(
    specific_variance(list()),
    "reserves this extractor name",
    fixed = TRUE
  )
  expect_error(
    latent_breeding_values(list()),
    "reserves this extractor name",
    fixed = TRUE
  )
  expect_error(
    genetic_loadings(list()),
    "reserves this extractor name",
    fixed = TRUE
  )
  # eigen_G() is no longer reserved (rotation-invariant eigenstructure); its
  # behaviour is covered in test-evolvability.R.
})

test_that("reserved factor-analytic extractors fail with rotation-aware scope", {
  fit <- hsquared:::hs_new_fit(
    call = quote(hsquared(y ~ animal(1 | id, pedigree = ped), data = dat)),
    spec = list(method = "REML", family = list(family = "gaussian")),
    payload = list(y = c(1, 2, 3)),
    result = list(
      genetic_covariance = matrix(c(1, 0.2, 0.2, 1), 2, 2),
      genetic_correlation = matrix(c(1, 0.2, 0.2, 1), 2, 2),
      converged = TRUE
    )
  )

  expect_error(
    genetic_loadings(fit),
    "planned, not implemented.*rotation-nonunique",
    perl = TRUE
  )
  expect_error(
    genetic_loadings(fit, rotate = "varimax"),
    "rotation controls are planned, not implemented",
    fixed = TRUE
  )
  # specific_variance() (Psi) is rotation-INVARIANT and identified, unlike
  # the loading-based extractors above; its message says so and points at
  # the bridge-activation gate instead of rotation-nonuniqueness (#218).
  expect_error(
    specific_variance(fit),
    "planned, not implemented.*rotation-INVARIANT",
    perl = TRUE
  )
  expect_error(
    latent_breeding_values(fit),
    "planned, not implemented.*rotation-nonunique",
    perl = TRUE
  )
  # eigen_G() now returns the rotation-invariant eigenstructure (see
  # test-evolvability.R); only the loading-based extractors stay reserved.
  expect_equal(eigen_G(fit)$values, eigen(fit$result$genetic_covariance)$values)
})

test_that("unsupported inference helpers fail with explicit scope", {
  fit <- hsquared:::hs_new_fit(
    call = quote(hsquared(y ~ animal(1 | id, pedigree = ped), data = dat)),
    spec = list(method = "REML", family = list(family = "gaussian")),
    payload = list(y = c(1, 2, 3)),
    result = list(
      variance_components = data.frame(
        component = c("animal", "residual"),
        estimate = c(0.4, 0.6)
      ),
      heritability = data.frame(term = "animal", estimate = 0.4),
      loglik = -4.2,
      df = 3L,
      nobs = 3L,
      converged = TRUE
    )
  )

  expect_error(
    stats::confint(fit),
    "Validated confidence intervals for `hsquared_fit` quantities are not implemented",
    fixed = TRUE
  )
  expect_error(
    stats::confint(fit),
    "heritability_interval()",
    fixed = TRUE
  )
  expect_error(
    stats::vcov(fit),
    "variance-covariance matrix is not implemented",
    fixed = TRUE
  )
  expect_error(
    stats::vcov(fit),
    "heritability_interval()",
    fixed = TRUE
  )
  expect_error(
    stats::profile(fit),
    "Profile-likelihood intervals.*planned, not implemented",
    perl = TRUE
  )
  expect_error(
    stats::anova(fit),
    "Likelihood-ratio / ANOVA comparison.*planned, not implemented",
    perl = TRUE
  )
})

test_that("hsquared_fit extractors fail loudly when a result field is absent", {
  fit <- hsquared:::hs_new_fit(
    spec = list(method = "REML", family = list(family = "gaussian")),
    payload = list(y = 1:2),
    result = list(converged = TRUE)
  )

  expect_error(
    variance_components(fit),
    "does not contain variance components",
    fixed = TRUE
  )
  expect_error(
    reliability(fit),
    "does not contain reliability estimates",
    fixed = TRUE
  )
  expect_error(
    accuracy(fit),
    "does not contain reliability estimates",
    fixed = TRUE
  )
  expect_error(
    qtl_table(fit),
    "does not contain QTL table",
    fixed = TRUE
  )
  expect_error(
    qtl_table(fit),
    "gwas(fit, markers)",
    fixed = TRUE
  )
  expect_error(
    marker_variance_explained(fit),
    "does not contain marker variance explained",
    fixed = TRUE
  )
  expect_error(
    gamma_matrix(fit),
    "does not contain a supplied metafounder `Gamma` matrix",
    fixed = TRUE
  )
  expect_error(
    metafounder_groups(fit),
    "does not contain supplied metafounder group assignments",
    fixed = TRUE
  )
  expect_error(
    metafounder_effects(fit),
    "reserved but not implemented",
    fixed = TRUE
  )
  expect_error(
    metafounder_effects(fit),
    "gamma_matrix",
    fixed = TRUE
  )
  expect_error(
    fitted(fit),
    "does not contain predictions",
    fixed = TRUE
  )
  expect_error(
    residuals(fit),
    "does not contain predictions",
    fixed = TRUE
  )
})

test_that("hsquared_fit nobs falls back to response payload", {
  fit <- hsquared:::hs_new_fit(
    spec = list(method = "REML", family = list(family = "gaussian")),
    payload = list(y = 1:4),
    result = list(converged = TRUE)
  )

  expect_equal(stats::nobs(fit), 4L)

  missing <- hsquared:::hs_new_fit(
    spec = list(method = "REML", family = list(family = "gaussian")),
    payload = list(),
    result = list(converged = TRUE)
  )
  expect_error(
    stats::nobs(missing),
    "does not contain number-of-observations metadata",
    fixed = TRUE
  )
})

test_that("fit_diagnostics tolerates scalar diagnostics payloads", {
  fit <- hsquared:::hs_new_fit(
    spec = list(method = "REML", family = list(family = "gaussian")),
    payload = list(y = 1:2),
    result = list(
      diagnostics = "scalar-status",
      converged = FALSE
    )
  )

  diagnostics <- fit_diagnostics(fit)
  expect_equal(
    diagnostics$value[diagnostics$metric == "diagnostics"],
    "scalar-status"
  )
})

test_that("fit_diagnostics flags a variance-component boundary solution", {
  boundary <- hsquared:::hs_new_fit(
    spec = list(method = "REML", family = list(family = "gaussian")),
    payload = list(y = 1:3),
    result = list(
      variance_components = data.frame(
        component = c("animal", "residual"),
        estimate = c(0, 1.2)
      ),
      converged = TRUE
    )
  )
  diag <- fit_diagnostics(boundary)
  expect_equal(diag$value[diag$metric == "at_boundary"], "TRUE")

  interior <- hsquared:::hs_new_fit(
    spec = list(method = "REML", family = list(family = "gaussian")),
    payload = list(y = 1:3),
    result = list(
      variance_components = data.frame(
        component = c("animal", "residual"),
        estimate = c(0.5, 0.5)
      ),
      converged = TRUE
    )
  )
  expect_equal(
    fit_diagnostics(interior)$value[
      fit_diagnostics(interior)$metric == "at_boundary"
    ],
    "FALSE"
  )
})

# hsquared#230: `search_boundary` (the non-Gaussian bridge's optimizer-bracket
# flag) is a distinct concept from `at_boundary` (an estimated variance share
# at/near zero). It is keyed off `result$boundary`, which only the non-Gaussian
# bridge route populates -- a Gaussian `engine = "fit"` result never carries a
# `boundary` field, so it must get no `search_boundary` row at all (not even
# NA), unlike `at_boundary`, which IS derived from Gaussian variance estimates.
test_that("fit_diagnostics reports search_boundary for a non-Gaussian-shaped result", {
  at_rail <- hsquared:::hs_new_fit(
    spec = list(method = "laplace", family = list(family = "binomial")),
    payload = list(y = 1:3),
    result = list(
      converged = TRUE,
      boundary = TRUE
    )
  )
  diag <- fit_diagnostics(at_rail)
  expect_equal(diag$value[diag$metric == "search_boundary"], "TRUE")
  expect_match(
    diag$value[diag$metric == "search_boundary_condition"],
    "log(initial$sigma_a2)",
    fixed = TRUE
  )
  expect_match(
    diag$value[diag$metric == "search_boundary_condition"],
    "HSquared.jl#327",
    fixed = TRUE
  )

  interior_ng <- hsquared:::hs_new_fit(
    spec = list(method = "laplace", family = list(family = "binomial")),
    payload = list(y = 1:3),
    result = list(
      converged = TRUE,
      boundary = FALSE
    )
  )
  diag_interior <- fit_diagnostics(interior_ng)
  expect_equal(
    diag_interior$value[diag_interior$metric == "search_boundary"],
    "FALSE"
  )
  condition_value <- diag_interior$value[
    diag_interior$metric == "search_boundary_condition"
  ]
  expect_equal(length(condition_value), 1L)
  expect_equal(condition_value, "interior")
})

test_that("fit_diagnostics gains no search_boundary row for a Gaussian fit without `boundary`", {
  gaussian_fit <- hsquared:::hs_new_fit(
    spec = list(method = "REML", family = list(family = "gaussian")),
    payload = list(y = 1:3),
    result = list(
      variance_components = data.frame(
        component = c("animal", "residual"),
        estimate = c(0.5, 0.5)
      ),
      converged = TRUE
    ),
    engine = "fit"
  )
  diag <- fit_diagnostics(gaussian_fit)
  expect_false("search_boundary" %in% diag$metric)
  expect_false("search_boundary_condition" %in% diag$metric)
})

test_that("accuracy requires reliability values on [0, 1]", {
  fit <- hsquared:::hs_new_fit(
    spec = list(method = "REML", family = list(family = "gaussian")),
    payload = list(y = 1:2),
    result = list(
      reliability = data.frame(id = c("a", "b"), value = c(0.5, 1.2))
    )
  )

  expect_error(
    accuracy(fit),
    "between 0 and 1",
    fixed = TRUE
  )

  malformed <- hsquared:::hs_new_fit(
    spec = list(method = "REML", family = list(family = "gaussian")),
    payload = list(y = 1:2),
    result = list(reliability = data.frame(id = c("a", "b"), se = c(0.1, 0.2)))
  )

  expect_error(
    accuracy(malformed),
    "with a `value` column",
    fixed = TRUE
  )
})

test_that("hsquared_fit residuals require response values", {
  fit <- hsquared:::hs_new_fit(
    spec = list(method = "REML", family = list(family = "gaussian")),
    payload = list(),
    result = list(predictions = data.frame(.fitted = c(1, 2)))
  )

  expect_equal(fitted(fit), c(1, 2))
  expect_error(
    residuals(fit),
    "does not contain response values",
    fixed = TRUE
  )
})

test_that("hsquared_fit residuals check fitted length", {
  fit <- hsquared:::hs_new_fit(
    spec = list(method = "REML", family = list(family = "gaussian")),
    payload = list(y = c(1, 2, 3)),
    result = list(predictions = data.frame(.fitted = c(1, 2)))
  )

  expect_error(
    residuals(fit),
    "same length",
    fixed = TRUE
  )
})

test_that("hs_fit_result names the fit target and does not say planned v0.1 contract", {
  fit <- hsquared:::hs_new_fit(
    spec = list(
      method = "REML",
      family = list(family = "gaussian"),
      target = "two_effect"
    ),
    payload = list(y = 1:2),
    result = list(converged = TRUE)
  )
  msg <- conditionMessage(
    tryCatch(variance_components(fit), error = function(e) e)
  )
  expect_match(msg, 'target = "two_effect"', fixed = TRUE)
  expect_match(msg, "does not contain variance components", fixed = TRUE)
  expect_false(grepl("planned v0.1 contract", msg, fixed = TRUE))
})

test_that("confint and vcov name live experimental interval extractors when present", {
  fit <- hsquared:::hs_new_fit(
    spec = list(
      method = "REML",
      family = list(family = "gaussian"),
      target = "ai_reml"
    ),
    payload = list(y = 1:3),
    result = list(
      heritability = data.frame(term = "animal", estimate = 0.4),
      heritability_interval = data.frame(
        estimate = 0.4,
        lower = 0.1,
        upper = 0.7,
        level = 0.95,
        method = "delta"
      ),
      variance_component_se = data.frame(component = "animal", se = 0.1),
      converged = TRUE
    )
  )
  expect_error(
    stats::confint(fit),
    "Use the experimental `heritability_interval()`",
    fixed = TRUE
  )
  expect_error(
    stats::confint(fit),
    "variance_component_standard_errors()",
    fixed = TRUE
  )
  expect_error(
    stats::confint(fit),
    "not coverage-calibrated",
    fixed = TRUE
  )
  expect_error(
    stats::vcov(fit),
    "heritability_interval()",
    fixed = TRUE
  )
  expect_false(grepl(
    "planned v0.1 contract",
    conditionMessage(tryCatch(stats::confint(fit), error = function(e) e)),
    fixed = TRUE
  ))
})

test_that("confint does not offer heritability_interval() on genomic fits", {
  fit <- hsquared:::hs_new_fit(
    spec = list(
      method = "REML",
      family = list(family = "gaussian"),
      target = "genomic"
    ),
    payload = list(y = 1:3),
    result = list(converged = TRUE)
  )
  msg <- conditionMessage(
    tryCatch(stats::confint(fit), error = function(e) e)
  )
  expect_match(msg, "not available for genomic", fixed = TRUE)
  expect_false(grepl(
    "Use the experimental `heritability_interval()`",
    msg,
    fixed = TRUE
  ))
})

test_that("print.hsquared_fit shows target, formula, and heritability peek", {
  fit <- hsquared:::hs_new_fit(
    call = quote(hsquared(y ~ animal(1 | id, pedigree = ped), data = dat)),
    spec = list(
      method = "REML",
      family = list(family = "gaussian"),
      target = "ai_reml"
    ),
    payload = list(y = 1:3),
    result = list(
      variance_components = data.frame(
        component = c("animal", "residual"),
        estimate = c(0.4, 0.6)
      ),
      heritability = data.frame(term = "animal", estimate = 0.4),
      converged = TRUE
    )
  )
  out <- paste(utils::capture.output(print(fit)), collapse = "\n")
  expect_match(out, "target: ai_reml", fixed = TRUE)
  expect_match(out, "formula: y ~ animal", fixed = TRUE)
  expect_match(out, "heritability: animal=0.4", fixed = TRUE)
})

test_that("heritability() and print() warn when the fit did not converge", {
  failed <- hsquared:::hs_new_fit(
    spec = list(
      method = "REML",
      family = list(family = "gaussian"),
      target = "ai_reml"
    ),
    payload = list(y = 1:4),
    result = list(
      variance_components = data.frame(
        component = c("animal", "residual"),
        estimate = c(1e-16, 1.2)
      ),
      heritability = data.frame(term = "animal", estimate = 9.48e-17),
      diagnostics = list(
        optimizer_status = "not_converged",
        iterations = 27L
      ),
      converged = FALSE
    )
  )

  expect_warning(
    h2 <- heritability(failed),
    "This `hsquared_fit` object did not converge. The heritability number is not an estimate; do not report it.",
    fixed = TRUE
  )
  expect_warning(
    heritability(failed),
    "failed-fit artefact, not evidence that heritability is zero",
    fixed = TRUE
  )
  expect_equal(h2$estimate, 9.48e-17)

  expect_warning(
    printed <- paste(utils::capture.output(print(failed)), collapse = "\n"),
    "This `hsquared_fit` object did not converge. The heritability number is not an estimate; do not report it.",
    fixed = TRUE
  )
  expect_match(printed, "converged: FALSE", fixed = TRUE)
  expect_match(
    printed,
    "heritability: not reportable (fit did not converge)",
    fixed = TRUE
  )
  expect_false(grepl("9.48", printed, fixed = TRUE))
})

test_that("heritability() and print() stay quiet on a converged interior fit", {
  ok <- hsquared:::hs_new_fit(
    spec = list(
      method = "REML",
      family = list(family = "gaussian"),
      target = "ai_reml"
    ),
    payload = list(y = 1:4),
    result = list(
      variance_components = data.frame(
        component = c("animal", "residual"),
        estimate = c(0.4, 0.6)
      ),
      heritability = data.frame(term = "animal", estimate = 0.4),
      diagnostics = list(optimizer_status = "converged"),
      converged = TRUE
    )
  )
  expect_silent(heritability(ok))
  expect_silent(invisible(utils::capture.output(print(ok))))
})
