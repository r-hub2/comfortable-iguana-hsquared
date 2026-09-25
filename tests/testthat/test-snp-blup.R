test_that("snp_blup is an accepted opt-in engine target", {
  expect_identical(hsquared:::hs_validate_julia_target("snp_blup"), "snp_blup")
})

test_that("snp_blup requires a marker-matrix genomic term", {
  ped <- data.frame(
    id = c("a", "b", "c"),
    sire = c(NA, NA, "a"),
    dam = c(NA, NA, "b")
  )
  dat <- data.frame(y = c(1, 2, 3), id = c("a", "b", "c"))
  expect_error(
    hsquared(
      y ~ animal(1 | id, pedigree = ped),
      data = dat,
      family = stats::gaussian(),
      control = hs_control(
        engine = "julia",
        engine_control = list(
          target = "snp_blup",
          variance_components = c(sigma_g2 = 1, sigma_e2 = 1)
        )
      )
    ),
    "requires a `genomic(1 | id, markers = M)` term",
    fixed = TRUE
  )
})

test_that("snp_blup rejects a supplied-Ginv genomic term (needs raw markers)", {
  ids <- paste0("g", 1:3)
  Ginv <- diag(3)
  dimnames(Ginv) <- list(ids, ids)
  dat <- data.frame(y = c(1, 2, 3), id = ids)
  expect_error(
    hsquared(
      y ~ genomic(1 | id, Ginv = Ginv),
      data = dat,
      family = stats::gaussian(),
      control = hs_control(
        engine = "julia",
        engine_control = list(
          target = "snp_blup",
          variance_components = c(sigma_g2 = 1, sigma_e2 = 1)
        )
      )
    ),
    "requires a `genomic(1 | id, markers = M)` term",
    fixed = TRUE
  )
})

test_that("snp_blup without supplied variances no longer requires them (routes to REML)", {
  # This checks the routing OFF-bridge (the old supplied-variance gate is gone, so
  # the unsupplied path reaches the bridge-availability check, not a "variance
  # required" error). The successful on-bridge REML fit is covered by the dedicated
  # live test below.
  testthat::skip_if(
    hsquared:::hs_julia_bridge_available(),
    "on-bridge REML routing is covered by the live REML SNP-BLUP test"
  )
  ids <- paste0("g", 1:5)
  set.seed(4)
  M <- matrix(stats::rbinom(5 * 20, 2, 0.3), 5, 20)
  rownames(M) <- ids
  dat <- data.frame(y = c(1, 2, 3, 4, 5), id = ids)
  err <- tryCatch(
    hsquared(
      y ~ genomic(1 | id, markers = M),
      data = dat,
      family = stats::gaussian(),
      control = hs_control(
        engine = "julia",
        engine_control = list(target = "snp_blup")
      )
    ),
    error = function(e) conditionMessage(e)
  )
  # not the old supplied-variance gate; the install-guidance error instead
  expect_false(grepl("variance_components", err, fixed = TRUE))
  expect_match(err, "Julia bridge|JuliaCall|HSquared.jl")
})

test_that("snp_blup variance components must be named sigma_g2 / sigma_e2", {
  expect_error(
    hsquared:::hs_validate_snp_blup_variances(c(sigma_a2 = 1, sigma_e2 = 1)),
    "sigma_g2",
    fixed = TRUE
  )
  expect_equal(
    hsquared:::hs_validate_snp_blup_variances(c(sigma_g2 = 2, sigma_e2 = 3)),
    c(sigma_g2 = 2, sigma_e2 = 3)
  )
})

test_that("snp_blup normalizer reports descriptive marker variance shares", {
  markers <- matrix(
    c(
      0,
      1,
      2,
      0,
      0,
      2,
      1,
      1,
      1
    ),
    nrow = 3,
    ncol = 3
  )
  payload <- list(
    ids = paste0("g", 1:3),
    markers = markers,
    marker_names = c("m1", "m2", "m3"),
    y = c(1, 2, 3),
    metadata = list(fixed_colnames = "(Intercept)")
  )
  raw <- list(
    beta = 1,
    gebv = c(0.1, -0.2, 0.3),
    marker_effects = c(1, 2, 3),
    p = colMeans(markers) / 2,
    fitted = c(1.1, 1.8, 3.3),
    nobs = 3
  )

  result <- hsquared:::hs_normalize_julia_snp_blup_result(
    raw,
    payload,
    c(sigma_g2 = 1, sigma_e2 = 2)
  )
  mve <- result$marker_variance_explained
  centered <- sweep(markers, 2, colMeans(markers), check.margin = FALSE)
  expected_contribution <- colMeans(centered^2) * raw$marker_effects^2

  expect_equal(mve$marker, c("m1", "m2", "m3"))
  expect_equal(mve$effect, raw$marker_effects)
  expect_equal(result$marker_allele_frequencies, raw$p)
  expect_equal(mve$contribution, expected_contribution)
  expect_equal(sum(mve$proportion, na.rm = TRUE), 1)
  expect_equal(mve$contribution[3], 0)
  expect_equal(mve$proportion[3], 0)
})

test_that("hsquared fits opt-in SNP-BLUP marker effects from a marker matrix", {
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for live SNP-BLUP."
  )

  set.seed(7)
  na <- 12
  ids <- paste0("g", seq_len(na))
  m <- 60
  M <- matrix(stats::rbinom(na * m, 2, 0.3), na, m)
  rownames(M) <- ids
  colnames(M) <- paste0("snp", seq_len(m))

  n <- 36
  rec <- rep(ids, length.out = n)
  dat <- data.frame(y = 5 + stats::rnorm(n, 0, 1), id = rec)

  fit <- hsquared(
    y ~ genomic(1 | id, markers = M),
    data = dat,
    family = stats::gaussian(),
    control = hs_control(
      engine = "julia",
      engine_control = list(
        target = "snp_blup",
        variance_components = c(sigma_g2 = 1.0, sigma_e2 = 2.0)
      )
    )
  )

  # Marker effects: one per marker, finite, labelled by the marker columns.
  me <- marker_effects(fit)
  expect_equal(nrow(me), m)
  expect_true(all(is.finite(me$effect)))
  expect_equal(me$marker[1], "snp1")

  # Descriptive fitted-marker shares: not scan p-values or QTL evidence.
  mve <- marker_variance_explained(fit)
  expect_equal(nrow(mve), m)
  expect_equal(mve$marker[1], "snp1")
  expect_true(all(is.finite(mve$contribution)))
  expect_true(all(mve$contribution >= 0))
  expect_equal(sum(mve$proportion, na.rm = TRUE), 1, tolerance = 1e-8)

  # Per-individual genomic breeding values: one per genotyped individual.
  bv <- breeding_values(fit)
  expect_equal(nrow(bv), na)
  expect_true(all(is.finite(bv$value)))

  # Supplied-variance solve: h2 is the supplied ratio, provenance "supplied".
  expect_equal(
    heritability(fit)$estimate,
    1.0 / (1.0 + 2.0),
    tolerance = 1e-8
  )
  expect_equal(
    fit_diagnostics(fit)$value[
      fit_diagnostics(fit)$metric == "variance_components_source"
    ],
    "supplied"
  )
})

test_that("hsquared fits REML-estimated SNP-BLUP when variances are unsupplied", {
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for live SNP-BLUP."
  )

  set.seed(13)
  na <- 25
  ids <- paste0("g", seq_len(na))
  m <- 50
  M <- matrix(stats::rbinom(na * m, 2, 0.3), na, m)
  rownames(M) <- ids
  colnames(M) <- paste0("snp", seq_len(m))
  # a genuine marker-driven genetic signal so REML estimates an interior sigma_g2
  Wc <- scale(M, center = TRUE, scale = FALSE)
  beta_true <- stats::rnorm(m, 0, 0.3)
  gv <- stats::setNames(as.numeric(Wc %*% beta_true), ids)
  n <- 75
  rec <- rep(ids, length.out = n)
  dat <- data.frame(y = 4 + gv[rec] + stats::rnorm(n, 0, 1), id = rec)

  # no supplied variance_components -> REML estimation
  fit <- hsquared(
    y ~ genomic(1 | id, markers = M),
    data = dat,
    family = stats::gaussian(),
    control = hs_control(
      engine = "julia",
      engine_control = list(target = "snp_blup")
    )
  )
  expect_s3_class(fit, "hsquared_fit")
  expect_equal(fit$spec$method, "SNP-BLUP-REML")

  vc <- variance_components(fit)
  expect_setequal(vc$component, c("genomic", "residual"))
  expect_true(all(is.finite(vc$estimate) & vc$estimate > 0))
  # variances are ESTIMATED, not the supplied default
  expect_equal(
    fit_diagnostics(fit)$value[
      fit_diagnostics(fit)$metric == "variance_components_source"
    ],
    "estimated_snp_blup_reml"
  )
  expect_equal(nrow(marker_effects(fit)), m)
  expect_equal(nrow(breeding_values(fit)), na)
  expect_true(is.finite(as.numeric(stats::logLik(fit))))
  # the fit genuinely converged (not just a finite loglik) ...
  expect_true(isTRUE(fit$result$converged))
  expect_equal(
    fit_diagnostics(fit)$value[
      fit_diagnostics(fit)$metric == "optimizer_status"
    ],
    "converged"
  )
  # ... and AIC/BIC are defined (df was set on the REML path)
  expect_true(is.finite(stats::AIC(fit)))
  # the variances are ESTIMATED, not the (1, 1) optimizer start
  expect_false(isTRUE(all.equal(
    sort(vc$estimate),
    c(1, 1),
    tolerance = 1e-3
  )))

  # parity: the estimated genomic variance matches a direct engine
  # fit_snp_blup_reml on the same record-level markers (the bridge left hsq_*).
  direct_sg2 <- JuliaCall::julia_eval(
    "HSquared.fit_snp_blup_reml(hsq_y, hsq_X, hsq_markers_rec).sigma_g2"
  )
  expect_equal(
    vc$estimate[vc$component == "genomic"],
    direct_sg2,
    tolerance = 1e-6
  )
})

test_that("GBLUP and REML SNP-BLUP give equivalent GEBVs on the same markers", {
  # The textbook GBLUP <-> SNP-BLUP equivalence (the twin V2-SNPBLUP pinned
  # property): fitting markers as a genomic relationship (GREML) or as marker
  # effects (SNP-BLUP), each REML-estimating its own variances, yields the same
  # per-individual genomic breeding values. The small residual is the ridge the
  # genomic path applies to G (genomic_relationship_inverse), not a discrepancy.
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for live SNP-BLUP."
  )

  set.seed(31)
  na <- 30
  ids <- paste0("g", seq_len(na))
  m <- 80
  M <- matrix(stats::rbinom(na * m, 2, 0.3), na, m)
  rownames(M) <- ids
  colnames(M) <- paste0("snp", seq_len(m))
  Wc <- scale(M, center = TRUE, scale = FALSE)
  gv <- stats::setNames(as.numeric(Wc %*% stats::rnorm(m, 0, 0.25)), ids)
  n <- 90
  rec <- rep(ids, length.out = n)
  dat <- data.frame(y = 3 + gv[rec] + stats::rnorm(n, 0, 1), id = rec)

  fit_gblup <- hsquared(
    y ~ genomic(1 | id, markers = M),
    data = dat,
    family = stats::gaussian(),
    control = hs_control(
      engine = "julia",
      engine_control = list(target = "genomic")
    )
  )
  fit_snp <- hsquared(
    y ~ genomic(1 | id, markers = M),
    data = dat,
    family = stats::gaussian(),
    control = hs_control(
      engine = "julia",
      engine_control = list(target = "snp_blup")
    )
  )
  bg <- breeding_values(fit_gblup)
  bs <- breeding_values(fit_snp)
  bg <- bg[order(bg$id), ]
  bs <- bs[order(bs$id), ]
  expect_equal(bg$id, bs$id)
  # GEBVs are equivalent: near-perfect correlation, small relative difference.
  expect_gt(stats::cor(bg$value, bs$value), 0.999)
  expect_lt(max(abs(bg$value - bs$value)) / stats::sd(bg$value), 0.02)
})
