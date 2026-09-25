# Live R<->engine parity guard for the plotting plot-data contract
# (docs/design/24-plotting-standard.md §5/§7; twin HSquared.jl#93). Skip-guarded
# on the local HSquared.jl bridge. As preparers land engine-side, add a parity
# case here so the R `autoplot` recompute path and the engine `*_plot_data`
# preparer cannot drift apart. (A Julia-free consumer-parity test -- payload vs
# recompute tidy frames -- lives in test-autoplot.R and always runs.)

test_that("genetic_correlation_plot_data engine preparer matches cov2cor + h2 [live]", {
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for the live plot-data parity test."
  )
  hsquared:::hs_julia_setup(hsquared:::hs_default_julia_project())

  g_cov <- matrix(
    c(
      2.0,
      0.6,
      0.2,
      0.6,
      1.0,
      0.1,
      0.2,
      0.1,
      0.5
    ),
    nrow = 3,
    byrow = TRUE
  )
  h2 <- c(0.05, 0.40, 0.30)
  JuliaCall::julia_assign("G_par", g_cov)
  JuliaCall::julia_assign("h_par", h2)
  JuliaCall::julia_command(
    "pd_par = HSquared.genetic_correlation_plot_data(G_par; heritabilities = h_par);"
  )
  corr_engine <- JuliaCall::julia_eval(
    "Matrix{Float64}(pd_par.genetic_correlations)"
  )
  h2_engine <- JuliaCall::julia_eval("collect(Float64, pd_par.heritabilities)")
  rot_engine <- JuliaCall::julia_eval("pd_par.rotation_invariant")

  # The engine's `genetic_correlations` is D^-1 G D^-1 == stats::cov2cor(G). The R
  # autoplot consumer trusts this contract for BOTH paths it can take: the
  # bridge-attached payload and the recompute fallback (the fit-stored
  # `genetic_correlation`, the same engine function evaluated at fit time).
  expect_equal(corr_engine, unname(stats::cov2cor(g_cov)), tolerance = 1e-12)
  expect_equal(h2_engine, h2, tolerance = 1e-12)
  expect_true(isTRUE(rot_engine))
})

test_that("autoplot consumes a live-marshalled genetic_correlation_plot_data payload [live]", {
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for the live plot-data parity test."
  )
  hsquared:::hs_julia_setup(hsquared:::hs_default_julia_project())

  g_cov <- matrix(
    c(
      2.0,
      0.6,
      0.2,
      0.6,
      1.0,
      0.1,
      0.2,
      0.1,
      0.5
    ),
    nrow = 3,
    byrow = TRUE
  )
  h2 <- c(0.05, 0.40, 0.30)
  JuliaCall::julia_assign("G_par2", g_cov)
  JuliaCall::julia_assign("h_par2", h2)
  JuliaCall::julia_command(paste(
    "pd2 = HSquared.genetic_correlation_plot_data(G_par2;",
    "traits = [\"t1\", \"t2\", \"t3\"], heritabilities = h_par2);"
  ))
  # Marshal the whole NamedTuple into R exactly as a bridge attachment would, to
  # check JuliaCall hands back a real R matrix (not a flattened vector) and that
  # the consumer unpack path works end-to-end on live engine output.
  pd_list <- JuliaCall::julia_eval(paste(
    "(traits = collect(String, pd2.traits),",
    "genetic_correlations = Matrix{Float64}(pd2.genetic_correlations),",
    "heritabilities = collect(Float64, pd2.heritabilities),",
    "rotation_invariant = pd2.rotation_invariant)"
  ))
  expect_true(is.matrix(pd_list$genetic_correlations))
  expect_equal(dim(pd_list$genetic_correlations), c(3L, 3L))

  fit <- structure(
    list(result = list(genetic_correlation_plot_data = pd_list)),
    class = "hsquared_fit"
  )
  p <- autoplot(fit, "g_matrix")
  cc <- unname(stats::cov2cor(g_cov))
  v_12 <- p$data$value[
    as.character(p$data$row) == "t1" & as.character(p$data$col) == "t2"
  ]
  expect_equal(v_12, cc[1, 2], tolerance = 1e-10)
  # t1 (h2 = 0.05 < 0.1) flags its off-diagonal cells.
  off <- p$data[as.character(p$data$row) != as.character(p$data$col), ]
  expect_true(any(off$low))
})

test_that("genetic_pca_plot_data engine preparer matches eigen_G scree [live]", {
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for the live plot-data parity test."
  )
  hsquared:::hs_julia_setup(hsquared:::hs_default_julia_project())

  g_cov <- matrix(c(2.0, 0.6, 0.6, 0.5), 2, 2)
  JuliaCall::julia_assign("G_pca", g_cov)
  JuliaCall::julia_command("pdg = HSquared.genetic_pca_plot_data(G_pca);")
  ev_engine <- JuliaCall::julia_eval("collect(Float64, pdg.eigenvalues)")
  ve_engine <- JuliaCall::julia_eval("collect(Float64, pdg.variance_explained)")
  rot_engine <- JuliaCall::julia_eval("pdg.rotation_invariant")
  not_loadings <- JuliaCall::julia_eval("pdg.is_eigenstructure_not_loadings")

  # The engine eigenvalues are the rotation-invariant eigenstructure; R's eigen_G
  # recompute (via the fit-stored G) sorts the same eigenvalues descending.
  r_ev <- sort(eigen(g_cov, symmetric = TRUE)$values, decreasing = TRUE)
  expect_equal(sort(ev_engine, decreasing = TRUE), r_ev, tolerance = 1e-10)
  expect_equal(sum(ve_engine), 1, tolerance = 1e-10)
  expect_true(isTRUE(rot_engine))
  expect_true(isTRUE(not_loadings))

  # Consume a marshalled payload end-to-end through autoplot.
  pdg <- JuliaCall::julia_eval(paste(
    "(eigenvalues = collect(Float64, pdg.eigenvalues),",
    "variance_explained = collect(Float64, pdg.variance_explained),",
    "axis_labels = collect(String, pdg.axis_labels),",
    "rotation_invariant = pdg.rotation_invariant)"
  ))
  fit <- structure(
    list(result = list(genetic_pca_plot_data = pdg)),
    class = "hsquared_fit"
  )
  p <- autoplot(fit, "g_geometry")
  expect_equal(p$data$eigenvalue, ev_engine, tolerance = 1e-10)
})

test_that("rr_genetic_variance_plot_data engine preparer matches R hs_rr_variance_values [live]", {
  # The RR parity guard the twin asked to co-own (#93 Q6): R's recompute
  # (hs_rr_variance_values) and the engine preparer must agree on v_g(t) on a
  # seeded K_g / standardized covariate grid.
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for the live plot-data parity test."
  )
  hsquared:::hs_julia_setup(hsquared:::hs_default_julia_project())

  k_g <- matrix(c(0.4, 0.05, 0.05, 0.1), 2, 2)
  ts <- seq(-1, 1, length.out = 5) # standardized covariate in [-1, 1]
  JuliaCall::julia_assign("Kg_rr", k_g)
  JuliaCall::julia_assign("ts_rr", ts)
  JuliaCall::julia_command(
    "pdrr = HSquared.rr_genetic_variance_plot_data(Kg_rr, ts_rr);"
  )
  gv_engine <- JuliaCall::julia_eval("collect(Float64, pdrr.genetic_variance)")
  cov_engine <- JuliaCall::julia_eval("collect(Float64, pdrr.covariate)")

  gv_r <- hsquared:::hs_rr_variance_values(k_g, ts, nrow(k_g))
  expect_equal(gv_engine, gv_r, tolerance = 1e-10)
  expect_equal(cov_engine, ts, tolerance = 1e-12)

  # Consume a marshalled payload end-to-end through autoplot (mirrors the other
  # live tests; needs heritability, so supply a residual to the preparer).
  JuliaCall::julia_command(
    "pdrr2 = HSquared.rr_genetic_variance_plot_data(Kg_rr, ts_rr; residual = 1.0);"
  )
  pdrr <- JuliaCall::julia_eval(paste(
    "(covariate = collect(Float64, pdrr2.covariate),",
    "genetic_variance = collect(Float64, pdrr2.genetic_variance),",
    "heritability = collect(Float64, pdrr2.heritability))"
  ))
  fit <- structure(
    list(result = list(rr_genetic_variance_plot_data = pdrr)),
    class = "hsquared_fit"
  )
  p <- autoplot(fit, "reaction_norm")
  gv_panel <- p$data$value[p$data$panel == "genetic variance"]
  expect_equal(gv_panel, gv_engine, tolerance = 1e-10)

  # rr_eigenfunctions_plot_data: marshal the m x k eigenfunctions matrix and feed
  # the eigenfunction figure; the consumed values match the engine matrix.
  JuliaCall::julia_command(
    "pdee = HSquared.rr_eigenfunctions_plot_data(Kg_rr, ts_rr);"
  )
  ef_engine <- JuliaCall::julia_eval("Matrix{Float64}(pdee.eigenfunctions)")
  pdee <- JuliaCall::julia_eval(paste(
    "(covariate = collect(Float64, pdee.covariate),",
    "eigenfunctions = Matrix{Float64}(pdee.eigenfunctions),",
    "variance_explained = collect(Float64, pdee.variance_explained),",
    "rotation_invariant = pdee.rotation_invariant)"
  ))
  fit_ef <- structure(
    list(result = list(rr_eigenfunctions_plot_data = pdee)),
    class = "hsquared_fit"
  )
  p_ef <- autoplot(fit_ef, "rr_eigenfunctions")
  expect_equal(
    sort(p_ef$data$value),
    sort(as.numeric(ef_engine)),
    tolerance = 1e-10
  )

  # rr_covariance_surface_plot_data: marshal the m x m surface and feed the
  # surface heatmap; the consumed values match the engine surface.
  JuliaCall::julia_command(
    "pdsf = HSquared.rr_covariance_surface_plot_data(Kg_rr, ts_rr);"
  )
  surf_engine <- JuliaCall::julia_eval("Matrix{Float64}(pdsf.surface)")
  pdsf <- JuliaCall::julia_eval(paste(
    "(covariate = collect(Float64, pdsf.covariate),",
    "surface = Matrix{Float64}(pdsf.surface),",
    "is_correlation = pdsf.is_correlation)"
  ))
  fit_sf <- structure(
    list(result = list(rr_covariance_surface_plot_data = pdsf)),
    class = "hsquared_fit"
  )
  p_sf <- autoplot(fit_sf, "rr_surface")
  expect_equal(
    sort(p_sf$data$value),
    sort(as.numeric(surf_engine)),
    tolerance = 1e-10
  )
})

test_that("variance_components_plot_data engine preparer feeds the variance forest [live]", {
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for the live plot-data parity test."
  )

  # A clean ~100-animal pedigree with simulated h2 = 0.4 so fit_ai_reml converges
  # to interior variance components (toy 3-5 animal datasets pin sigma_a2 -> 0).
  ped <- hs_sim_pedigree(n_founder = 20, n_per_gen = 40, n_gen = 2, seed = 7)
  dat <- hs_sim_genedrop_phenotypes(
    ped,
    sigma_a2 = 0.4,
    sigma_e2 = 0.6,
    seed = 7
  )
  fit <- hsquared(
    y ~ animal(1 | id, pedigree = ped),
    data = dat,
    family = stats::gaussian(),
    REML = TRUE
  )
  vcpd <- fit$result$variance_components_plot_data
  expect_false(is.null(vcpd))
  # Engine field names match the consumer's contract, and the estimates equal the
  # fit's variance components / heritability. This is the actual bridge-attached
  # payload, not a manual post-fit preparer call.
  expect_equal(vcpd$term, c("sigma_a2", "sigma_e2", "h2"))
  expect_equal(
    vcpd$estimate[1:2],
    variance_components(fit)$estimate,
    tolerance = 1e-8
  )

  p <- autoplot(fit, "variance")
  expect_s3_class(p, "ggplot")
  expect_true(all(
    c("sigma_a2", "sigma_e2", "h2") %in% as.character(p$data$term)
  ))
  expect_equal(
    p$data$estimate[as.character(p$data$term) == "sigma_a2"],
    vcpd$estimate[1],
    tolerance = 1e-10
  )

  # The preparer ships `lo`/`hi` as NaN where the interval is unavailable; confirm
  # a Julia Float64 NaN survives the JuliaCall bridge as R NaN (positions 1/3),
  # while a finite value passes through (position 2). The consumer turns the NaN
  # into NA -> no whisker, and keeps the finite one.
  mixed <- JuliaCall::julia_eval("collect(Float64, [NaN, 1.0, NaN])")
  expect_true(is.nan(mixed[1]))
  expect_false(is.nan(mixed[2]))
  fit_nan <- structure(
    list(
      result = list(
        variance_components_plot_data = list(
          term = c("sigma_a2", "sigma_e2", "h2"),
          estimate = c(0.5, 0.5, 0.4),
          lo = mixed,
          hi = mixed,
          panel = c(
            "variance components",
            "variance components",
            "heritability"
          ),
          interval_status = "none"
        )
      )
    ),
    class = "hsquared_fit"
  )
  lo_out <- autoplot(fit_nan, "variance")$data$lo
  expect_true(all(is.na(lo_out[c(1, 3)]))) # engine NaN -> R NA
  expect_equal(lo_out[2], 1.0) # finite value survives
})

test_that("breeding_values_plot_data engine preparer matches the R recompute [live]", {
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for live plot-data parity."
  )

  ped <- hs_sim_pedigree(n_founder = 12, n_per_gen = 24, n_gen = 2, seed = 41)
  dat <- hs_sim_genedrop_phenotypes(
    ped,
    sigma_a2 = 0.5,
    sigma_e2 = 0.5,
    seed = 41
  )
  fit <- hsquared(
    y ~ animal(1 | id, pedigree = ped),
    data = dat,
    family = stats::gaussian(),
    REML = TRUE
  )

  bvpd <- fit$result$breeding_values_plot_data
  expect_false(is.null(bvpd))
  expect_equal(bvpd$pev_scale, "validation")
  expect_false("trait" %in% names(bvpd))

  # The recompute fallback agrees with the attached engine payload.
  bv <- breeding_values(fit)
  pev <- prediction_error_variance(fit)
  o <- match(bvpd$id, bv$id)
  expect_false(anyNA(o))
  expect_equal(bv$value[o], bvpd$value, tolerance = 1e-8)
  po <- match(bvpd$id, pev$id)
  expect_equal(pev$value[po], bvpd$pev, tolerance = 1e-8)

  # consume: the bridge-attached payload draws the same EBVs
  p <- autoplot(fit, "breeding_values")
  expect_s3_class(p, "ggplot")
  expect_setequal(round(p$data$value, 8), round(bvpd$value, 8))
})
