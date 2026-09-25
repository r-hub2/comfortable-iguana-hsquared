# Opt-in, experimental random-regression (reaction-norm) bridge. The parser /
# payload / normalizer / extractor tests run without Julia; the live fit test is
# skip-guarded on the local HSquared.jl bridge. Provisional grammar (HSquared.jl#61).

# Tiny long-format fixture: 4 base animals with repeated records across age.
hs_rr_fixture <- function() {
  ped <- data.frame(
    id = c("s1", "d1", "s2", "d2", "a", "b", "c", "d"),
    sire = c(NA, NA, NA, NA, "s1", "s1", "s2", "s2"),
    dam = c(NA, NA, NA, NA, "d1", "d1", "d2", "d2"),
    stringsAsFactors = FALSE
  )
  animals <- c("a", "b", "c", "d")
  ages <- c(1, 3, 5)
  long <- expand.grid(
    id = animals,
    age = ages,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  long$sex <- ifelse(long$id %in% c("a", "c"), "m", "f")
  set.seed(42)
  base <- c(a = 2.0, b = 2.4, c = 1.8, d = 2.2)
  slope <- c(a = 0.3, b = 0.1, c = 0.25, d = 0.05)
  long$weight <- base[long$id] +
    slope[long$id] * long$age +
    stats::rnorm(nrow(long), sd = 0.1)
  list(ped = ped, data = long)
}

test_that("rr() parser builds a random-regression spec and payload shape", {
  fx <- hs_rr_fixture()
  ped <- fx$ped
  dat <- fx$data

  spec <- hsquared:::hs_build_model_spec(
    weight ~ sex + animal(rr(age, order = 2) | id, pedigree = ped),
    data = fx$data,
    family = stats::gaussian(),
    REML = TRUE
  )

  expect_equal(spec$random$animal$design, "random_regression")
  rr <- spec$random$animal$random_regression
  expect_equal(rr$covariate, "age")
  expect_equal(rr$order, 2L)
  expect_equal(rr$lower, 1)
  expect_equal(rr$upper, 5)
  expect_equal(rr$values, fx$data$age)
  expect_match(
    spec$bridge$target,
    "fit_random_regression_reml",
    fixed = TRUE
  )

  payload <- hsquared:::hs_build_bridge_payload(spec)
  expect_false(is.null(payload$y))
  expect_null(payload$Y)
  expect_equal(payload$metadata$response_type, "random_regression")
  expect_equal(payload$random_regression$order, 2L)
  expect_equal(payload$metadata$random_regression$covariate, "age")
  expect_equal(payload$metadata$random_regression$lower, 1)
  expect_equal(payload$metadata$random_regression$upper, 5)
  # Long-format record incidence: one row per record, one column per animal.
  expect_s4_class(payload$Z, "dgCMatrix")
  expect_equal(dim(payload$Z), c(nrow(fx$data), nrow(fx$ped)))
})

test_that("rr() default order is 2 (intercept + slope)", {
  fx <- hs_rr_fixture()
  ped <- fx$ped
  dat <- fx$data
  spec <- hsquared:::hs_build_model_spec(
    weight ~ animal(rr(age) | id, pedigree = ped),
    data = fx$data,
    family = stats::gaussian(),
    REML = TRUE
  )
  expect_equal(spec$random$animal$random_regression$order, 2L)
})

test_that("rr() rejects unsupported syntax with named pointers", {
  fx <- hs_rr_fixture()
  ped <- fx$ped
  dat <- fx$data

  # Unknown covariate column.
  expect_error(
    hsquared:::hs_build_model_spec(
      weight ~ animal(rr(missing_col, order = 2) | id, pedigree = ped),
      data = fx$data,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "was not found in `data`",
    fixed = TRUE
  )

  # Planned rr() argument.
  expect_error(
    hsquared:::hs_build_model_spec(
      weight ~ animal(rr(age, order = 2, type = "spline") | id, pedigree = ped),
      data = fx$data,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "planned, not implemented",
    fixed = TRUE
  )

  # Bad order.
  expect_error(
    hsquared:::hs_build_model_spec(
      weight ~ animal(rr(age, order = 0) | id, pedigree = ped),
      data = fx$data,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "must be a single positive integer",
    fixed = TRUE
  )

  # A genuinely non-intercept, non-rr() left-hand side stays rejected.
  expect_error(
    hsquared:::hs_build_model_spec(
      weight ~ animal(age | id, pedigree = ped),
      data = fx$data,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "random-intercept syntax",
    fixed = TRUE
  )
})

test_that("random_regression target is explicitly opt-in", {
  fx <- hs_rr_fixture()
  ped <- fx$ped
  dat <- fx$data

  # Default engine = "fit" rejects rr() as experimental/opt-in.
  expect_error(
    hsquared(
      weight ~ animal(rr(age, order = 2) | id, pedigree = ped),
      data = fx$data
    ),
    "experimental",
    fixed = TRUE
  )

  # engine = "julia" with the wrong target rejects rr().
  expect_error(
    hsquared(
      weight ~ animal(rr(age, order = 2) | id, pedigree = ped),
      data = fx$data,
      control = hs_control(engine = "julia")
    ),
    "requires the opt-in `target = \"random_regression\"`",
    fixed = TRUE
  )

  # target = "random_regression" without an rr() term errors.
  expect_error(
    hsquared(
      weight ~ animal(1 | id, pedigree = ped),
      data = fx$data,
      control = hs_control(
        engine = "julia",
        engine_control = list(target = "random_regression")
      )
    ),
    "requires an",
    fixed = TRUE
  )

  expect_equal(
    hsquared:::hs_validate_julia_target("random_regression"),
    "random_regression"
  )
})

test_that("random-regression result normalizer exposes K_g, coefficients, trajectories", {
  payload <- structure(
    list(
      y = rep(c(1, 2, 3), 4),
      X = matrix(1, nrow = 12L, ncol = 1L),
      ids = c("a", "b", "c", "d"),
      family = "gaussian",
      random_regression = list(
        covariate = "age",
        order = 2L,
        lower = 1,
        upper = 5
      ),
      metadata = list(
        fixed_colnames = "(Intercept)",
        random_regression = list(
          covariate = "age",
          order = 2L,
          lower = 1,
          upper = 5
        )
      )
    ),
    class = c("hs_bridge_payload", "list")
  )
  raw <- list(
    K_g = matrix(c(0.8, 0.1, 0.1, 0.3), 2L),
    sigma_e2 = 0.5,
    beta = 2.1,
    coef_ids = payload$ids,
    coef_values = matrix(seq(0.1, 0.8, length.out = 8L), nrow = 4L),
    loglik = -10.2,
    converged = TRUE,
    iterations = 33L,
    ncoef = 2L
  )

  result <- hsquared:::hs_normalize_random_regression_result(raw, payload)
  fit <- hsquared:::hs_new_fit(
    spec = list(
      method = "REML",
      family = list(family = "gaussian"),
      target = "random_regression"
    ),
    payload = payload,
    result = result
  )

  # K_g is k x k and labelled.
  K_g <- rr_covariance(fit)
  expect_equal(dim(K_g), c(2L, 2L))
  expect_equal(rownames(K_g), c("legendre0", "legendre1"))

  # Predicted coefficients: q x k in long format.
  rc <- random_coefficients(fit)
  expect_equal(nrow(rc), 8L)
  expect_equal(sort(unique(rc$coefficient)), c("legendre0", "legendre1"))

  # Genetic-variance trajectory matches phi(t)' K_g phi(t) at the endpoints.
  vg <- rr_genetic_variance(fit, at = c(1, 5))
  expect_equal(nrow(vg), 2L)
  phi_lo <- hsquared:::hs_legendre_basis(-1, 2L)
  phi_hi <- hsquared:::hs_legendre_basis(1, 2L)
  expect_equal(
    vg$value,
    c(
      drop(phi_lo %*% raw$K_g %*% phi_lo),
      drop(phi_hi %*% raw$K_g %*% phi_hi)
    ),
    tolerance = 1e-10
  )

  # Heritability trajectory in [0, 1].
  h2 <- rr_heritability(fit)
  expect_true(all(h2$value >= 0 & h2$value <= 1))
  expect_equal(
    h2$value[1],
    vg$value[1] / (vg$value[1] + raw$sigma_e2),
    tolerance = 1e-10
  )

  # Correlation surface: symmetric with unit diagonal.
  corr <- rr_correlation(fit, at = c(1, 3, 5))
  expect_equal(dim(corr), c(3L, 3L))
  expect_equal(unname(diag(corr)), rep(1, 3L), tolerance = 1e-10)
  expect_equal(unname(corr), unname(t(corr)), tolerance = 1e-10)

  # Eigen-functions: rotation-invariant decomposition of K_g.
  ef <- rr_eigenfunctions(fit, at = c(1, 3, 5))
  expect_named(
    ef,
    c(
      "covariate",
      "eigenvalues",
      "variance_explained",
      "eigen_coefficients",
      "eigenfunctions"
    )
  )
  expect_equal(length(ef$eigenvalues), 2L)
  expect_equal(
    ef$eigenvalues,
    sort(eigen(raw$K_g, symmetric = TRUE)$values, decreasing = TRUE),
    tolerance = 1e-10
  )
  expect_equal(sum(ef$variance_explained), 1, tolerance = 1e-10)
  expect_equal(dim(ef$eigen_coefficients), c(2L, 2L))
  expect_equal(nrow(ef$eigenfunctions), 3L * 2L)
  expect_equal(sort(unique(ef$eigenfunctions$axis)), c(1L, 2L))

  # reaction-norm autoplot returns a faceted ggplot (variance + heritability).
  p_rn <- autoplot(fit, "reaction_norm")
  expect_s3_class(p_rn, "ggplot")
  # The reaction-norm trajectories are functionals of K_g, so the figure is in
  # the plotting-standard's rotation-invariant set (24-plotting-standard.md §3
  # binding rule): rotation_status MUST be "rotation_invariant".
  m_rn <- attr(p_rn, "hsquared_meta")
  expect_equal(m_rn$type, "reaction_norm")
  expect_equal(m_rn$rotation_status, "rotation_invariant")
  expect_equal(m_rn$interval_status, "descriptive")
  # Recompute path: the panels equal the trajectory extractors (not just a ggplot).
  gv_rn <- rr_genetic_variance(fit)
  h2_rn <- rr_heritability(fit)
  expect_equal(
    p_rn$data$value[p_rn$data$panel == "genetic variance"],
    gv_rn$value
  )
  expect_equal(
    p_rn$data$value[p_rn$data$panel == "heritability"],
    h2_rn$value
  )

  # Payload auto-detect: attach an engine `rr_genetic_variance_plot_data` payload
  # with a DISTINGUISHABLE covariate grid so consumption vs recompute is provable.
  fit_pd <- fit
  fit_pd$result$rr_genetic_variance_plot_data <- list(
    covariate = c(99, 98, 97),
    value = c(1, 1, 1),
    heritability = c(0.5, 0.5, 0.5)
  )
  # default at = NULL -> the payload is consumed (its covariate grid is plotted).
  p_consume <- autoplot(fit_pd, "reaction_norm")
  expect_equal(
    sort(unique(p_consume$data$covariate)),
    c(97, 98, 99)
  )
  # a custom `at` must BYPASS the payload and recompute on the user's grid
  # (in-range: the fitted covariate range here is [1, 5], see #213).
  p_bypass <- autoplot(fit_pd, "reaction_norm", at = c(2, 3, 4))
  expect_equal(
    sort(unique(p_bypass$data$covariate)),
    c(2, 3, 4)
  )
  # rename precedence: with both fields, `value` wins over `genetic_variance`.
  fit_both <- fit
  fit_both$result$rr_genetic_variance_plot_data <- list(
    covariate = c(1, 2),
    value = c(5, 5),
    genetic_variance = c(9, 9),
    heritability = c(0.5, 0.5)
  )
  gvb <- autoplot(fit_both, "reaction_norm")$data
  expect_equal(gvb$value[gvb$panel == "genetic variance"], c(5, 5))
  # partial payload (no heritability) -> all-or-nothing guard recomputes.
  fit_partial <- fit
  fit_partial$result$rr_genetic_variance_plot_data <- list(
    covariate = c(99, 98),
    value = c(1, 1)
  )
  expect_false(any(
    autoplot(fit_partial, "reaction_norm")$data$covariate %in%
      c(98, 99)
  ))

  # rr_eigenfunctions autoplot (recompute path): faceted psi_j(t) curves whose
  # values equal the rr_eigenfunctions() extractor; rotation-invariant meta.
  p_ef <- autoplot(fit, "rr_eigenfunctions")
  expect_s3_class(p_ef, "ggplot")
  ef_rn <- rr_eigenfunctions(fit)
  expect_equal(
    sort(p_ef$data$value),
    sort(ef_rn$eigenfunctions$value),
    tolerance = 1e-10
  )
  expect_equal(
    attr(p_ef, "hsquared_meta")$rotation_status,
    "rotation_invariant"
  )

  # rr_surface autoplot (recompute path): the covariance surface S = phi K_g phi'
  # over the grid; correlation = TRUE has a unit diagonal.
  p_surf <- autoplot(fit, "rr_surface")
  expect_s3_class(p_surf, "ggplot")
  k_g_rn <- rr_covariance(fit)
  pts_rn <- hsquared:::hs_rr_eval_points(fit, NULL, 25L)
  phi_rn <- hsquared:::hs_legendre_design(pts_rn$t, pts_rn$order)
  surf_rn <- phi_rn %*% k_g_rn %*% t(phi_rn)
  expect_equal(
    sort(p_surf$data$value),
    sort(as.numeric(surf_rn)),
    tolerance = 1e-10
  )
  expect_equal(
    attr(p_surf, "hsquared_meta")$rotation_status,
    "rotation_invariant"
  )
  p_corr <- autoplot(fit, "rr_surface", correlation = TRUE)
  diag_corr <- p_corr$data$value[
    p_corr$data$covariate_i == p_corr$data$covariate_j
  ]
  expect_true(all(abs(diag_corr - 1) < 1e-8))

  # Generic fit S3 surfaces work on a random-regression fit.
  expect_equal(stats::nobs(fit), 12L)
  expect_s3_class(stats::logLik(fit), "logLik")
  expect_equal(nrow(fixef(fit)), 1L)
  expect_equal(
    fit_diagnostics(fit)$value[
      fit_diagnostics(fit)$metric == "variance_components_source"
    ],
    "estimated_random_regression_reml"
  )
})

test_that("random-regression extractors reject non-RR fits", {
  payload <- list(y = 1:3, metadata = list(fixed_colnames = "(Intercept)"))
  fit <- hsquared:::hs_new_fit(
    spec = list(
      method = "REML",
      family = list(family = "gaussian"),
      target = "ai_reml"
    ),
    payload = payload,
    result = list(converged = TRUE)
  )
  expect_error(rr_covariance(fit), "random-regression", fixed = TRUE)
  expect_error(rr_genetic_variance(fit), "random-regression", fixed = TRUE)
  expect_error(rr_heritability(fit), "random-regression", fixed = TRUE)
  expect_error(rr_correlation(fit), "random-regression", fixed = TRUE)
  expect_error(random_coefficients(fit), "random-regression", fixed = TRUE)
  expect_error(rr_covariance(42), "random-regression", fixed = TRUE)
})

test_that("R Legendre basis matches the engine convention", {
  # phi_0 = sqrt(1/2); phi_1(t) = sqrt(3/2) t; phi_2(t) = sqrt(5/2)(3t^2-1)/2.
  expect_equal(
    hsquared:::hs_legendre_basis(0, 3L),
    c(
      sqrt(1 / 2),
      0,
      sqrt(5 / 2) * (-1 / 2)
    ),
    tolerance = 1e-12
  )
  expect_equal(
    hsquared:::hs_legendre_basis(1, 3L),
    c(
      sqrt(1 / 2),
      sqrt(3 / 2),
      sqrt(5 / 2)
    ),
    tolerance = 1e-12
  )
  # standardize_covariate maps [lower, upper] -> [-1, 1].
  expect_equal(
    hsquared:::hs_standardize_covariate(c(1, 3, 5), 1, 5),
    c(-1, 0, 1)
  )
})

test_that("hsquared can use the opt-in experimental random-regression bridge", {
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for live random-regression bridge smoke."
  )

  fx <- hs_rr_fixture()
  ped <- fx$ped
  dat <- fx$data

  fit <- hsquared(
    weight ~ sex + animal(rr(age, order = 2) | id, pedigree = ped),
    data = fx$data,
    family = stats::gaussian(),
    REML = TRUE,
    control = hs_control(
      engine = "julia",
      engine_control = list(target = "random_regression", iterations = 400L)
    )
  )

  expect_s3_class(fit, "hsquared_fit")
  expect_equal(fit$spec$target, "random_regression")

  # K_g is k x k and positive definite.
  K_g <- rr_covariance(fit)
  expect_equal(dim(K_g), c(2L, 2L))
  expect_true(all(eigen(K_g, symmetric = TRUE, only.values = TRUE)$values > 0))

  # Residual variance positive.
  expect_true(fit$result$residual_variance > 0)

  # Predicted coefficients: q x k.
  rc <- random_coefficients(fit)
  expect_equal(nrow(rc), nrow(fx$ped) * 2L)

  # Heritability trajectory in [0, 1].
  h2 <- rr_heritability(fit)
  expect_true(all(h2$value >= 0 & h2$value <= 1))

  expect_false(is.null(fit$result$rr_genetic_variance_plot_data))
  expect_false(is.null(fit$result$rr_eigenfunctions_plot_data))
  expect_false(is.null(fit$result$rr_covariance_surface_plot_data))
  gv_payload <- fit$result$rr_genetic_variance_plot_data
  gv_fallback <- rr_genetic_variance(fit)
  expect_equal(range(gv_payload$covariate), c(1, 5))
  expect_equal(gv_payload$value, gv_fallback$value, tolerance = 1e-8)
  expect_equal(gv_payload$heritability, h2$value, tolerance = 1e-8)
  expect_equal(
    dim(fit$result$rr_eigenfunctions_plot_data$eigenfunctions),
    c(25L, 2L)
  )
  expect_equal(
    dim(fit$result$rr_covariance_surface_plot_data$surface),
    c(25L, 25L)
  )
  expect_equal(
    length(unique(autoplot(fit, "reaction_norm", n = 7L)$data$covariate)),
    7L
  )

  # rr_eigenfunctions matches the engine `rr_eigenfunctions` element-wise.
  ef <- rr_eigenfunctions(fit, at = NULL, n = 9L)
  pts <- hsquared:::hs_rr_eval_points(fit, NULL, 9L)
  JuliaCall::julia_assign("Kg_ef", as.matrix(K_g))
  JuliaCall::julia_assign("ts_ef", pts$t)
  JuliaCall::julia_command("out_ef = HSquared.rr_eigenfunctions(Kg_ef, ts_ef);")
  ev_J <- JuliaCall::julia_eval("collect(Float64, out_ef.eigenvalues)")
  ve_J <- JuliaCall::julia_eval("collect(Float64, out_ef.variance_explained)")
  psi_J <- JuliaCall::julia_eval("Matrix{Float64}(out_ef.eigenfunctions)")
  expect_equal(ef$eigenvalues, ev_J, tolerance = 1e-8)
  expect_equal(ef$variance_explained, ve_J, tolerance = 1e-8)
  psi_R <- matrix(ef$eigenfunctions$value, nrow = length(pts$at))
  expect_equal(psi_R, psi_J, tolerance = 1e-8)

  expect_true(isTRUE(fit$result$converged))
})

# A minimal normalized random-regression fit, built through the real result
# normalizer so the fence is exercised against the payload shape the bridge
# actually produces rather than a hand-shaped list.
hs_rr_mock_fit <- function() {
  payload <- structure(
    list(
      y = rep(c(1, 2, 3), 4),
      X = matrix(1, nrow = 12L, ncol = 1L),
      ids = c("a", "b", "c", "d"),
      family = "gaussian",
      random_regression = list(
        covariate = "age",
        order = 2L,
        lower = 1,
        upper = 5
      ),
      metadata = list(
        fixed_colnames = "(Intercept)",
        random_regression = list(
          covariate = "age",
          order = 2L,
          lower = 1,
          upper = 5
        )
      )
    ),
    class = c("hs_bridge_payload", "list")
  )
  raw <- list(
    K_g = matrix(c(0.8, 0.1, 0.1, 0.3), 2L),
    sigma_e2 = 0.5,
    beta = 2.1,
    coef_ids = payload$ids,
    coef_values = matrix(seq(0.1, 0.8, length.out = 8L), nrow = 4L),
    loglik = -10.2,
    converged = TRUE,
    iterations = 33L,
    ncoef = 2L
  )
  hsquared:::hs_new_fit(
    spec = list(
      method = "REML",
      family = list(family = "gaussian"),
      target = "random_regression"
    ),
    payload = payload,
    result = hsquared:::hs_normalize_random_regression_result(raw, payload)
  )
}

# A21/C1+C2 regression guard. The `rr()` argument-naming predicate is `order =`,
# and the reaction-norm heritability accessor is `rr_heritability()`. Two public
# surfaces once carried `rr(t, k = 2)` and credited `heritability()` with the
# h2(t) curve; both spellings are refused by the package, so they were false
# statements about a route these same surfaces call covered.
test_that("rr() argument naming and the RR heritability accessor are pinned", {
  fx <- hs_rr_fixture()
  ped <- fx$ped

  # Accepted spelling: `order =` reaches the opt-in fence, i.e. the grammar
  # itself is accepted and only the opt-in gate stops the call.
  expect_error(
    hsquared::hsquared(
      weight ~ animal(rr(age, order = 2) | id, pedigree = ped),
      data = fx$data,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "Closest working call",
    fixed = TRUE
  )
  err <- tryCatch(
    hsquared::hsquared(
      weight ~ animal(rr(age, order = 2) | id, pedigree = ped),
      data = fx$data,
      family = stats::gaussian(),
      REML = TRUE
    ),
    error = function(e) conditionMessage(e)
  )
  expect_true(grepl("`rr(...)` is not on the default path", err, fixed = TRUE))
  expect_true(grepl("target = \"random_regression\"", err, fixed = TRUE))
  expect_true(grepl("k = 2", err, fixed = TRUE))

  # Rejected spelling: `k =` is refused, and the error names the implemented
  # spelling rather than only the unsupported one.
  expect_error(
    hsquared::hsquared(
      weight ~ animal(rr(age, k = 2) | id, pedigree = ped),
      data = fx$data,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "rr(covariate, order = k)",
    fixed = TRUE
  )

  # `heritability()` is not the RR accessor: it errors and points at
  # `rr_heritability()` instead of falling through to the generic
  # "planned v0.1 contract" miss, which reads as "not implemented".
  fit <- hs_rr_mock_fit()
  expect_error(heritability(fit), "rr_heritability()", fixed = TRUE)
  expect_error(heritability(fit), "random-regression", fixed = TRUE)
  expect_false(grepl(
    "planned v0.1 contract",
    conditionMessage(
      tryCatch(heritability(fit), error = function(e) e)
    ),
    fixed = TRUE
  ))

  # And the implemented accessor still returns the curve.
  expect_gt(nrow(rr_heritability(fit)), 1L)
})

# Claim-surface contract: `docs/` and `vignettes/articles/` are .Rbuildignore'd,
# so these are source-tree-only checks (skipped under `R CMD check`).
test_that("public claim surfaces use the implemented rr() grammar and accessor", {
  surfaces <- c(
    file.path("docs", "design", "06-public-claims-register.md"),
    file.path("docs", "design", "capability-status.md"),
    file.path("vignettes", "articles", "model-status.Rmd")
  )
  for (rel in surfaces) {
    path <- file.path(testthat::test_path("..", ".."), rel)
    skip_if_not(file.exists(path), paste(rel, "is not in the build tarball"))
    lines <- readLines(path, warn = FALSE)

    # No `rr(..., k = ...)` — the parser refuses that argument name.
    expect_identical(
      grep("rr\\([^)]*\\bk\\s*=", lines),
      integer(0),
      label = paste(rel, "uses rr(covariate, order = k)")
    )

    # Any surface crediting an RR heritability accessor must name the real one.
    rr_rows <- grep(
      "random.regression|reaction.norm",
      lines,
      ignore.case = TRUE
    )
    if (length(rr_rows) && any(grepl("h.\\(t\\)", lines[rr_rows]))) {
      expect_true(
        any(grepl("rr_heritability", lines, fixed = TRUE)),
        label = paste(rel, "names rr_heritability() as the h2(t) accessor")
      )
    }
  }
})
