# Opt-in, experimental direct-maternal correlated model (Phase 4):
# 2x2 G_dm with direct additive genetic + maternal additive genetic effects,
# fitted via target = "direct_maternal".
#
# Tests cover:
#   (A) target routing: "direct_maternal" accepted / "maternal_genetic" without
#       correct target gives a helpful error
#   (B) payload-v2 emitter: hs_build_bridge_payload() produces two pedigree
#       blocks (animal + maternal) for a maternal_genetic formula — the fit
#       function reassembles these into a correlated Julia block
#   (C) accessor fences: heritability(), direct_heritability(), genetic_correlation(),
#       direct_variance(), partner_variance(), direct_maternal_covariance()
#   (D) live parity test (skip-guarded: only runs when JuliaCall + HSquared.jl
#       are available)

# ---- shared fixtures --------------------------------------------------------

make_dm_ped <- function() {
  data.frame(
    id   = c("g1", "g2", "g3", "g4", "g5", "g6"),
    sire = c(NA,   NA,   "g1", "g1", "g2", "g2"),
    dam  = c(NA,   NA,   "g2", "g2", "g1", "g1"),
    stringsAsFactors = FALSE
  )
}

make_dm_dat <- function() {
  data.frame(
    y   = c(1.2, 0.8, 2.1, 1.9, 1.5, 2.3, 1.7, 0.9),
    id  = c("g3", "g4", "g5", "g6", "g3", "g4", "g5", "g6"),
    dam = c("g2", "g2", "g1", "g1", "g2", "g2", "g1", "g1"),
    stringsAsFactors = FALSE
  )
}

# ---- (A) target routing ----------------------------------------------------- #

test_that("target 'direct_maternal' is accepted by hs_validate_julia_target", {
  expect_equal(
    hsquared:::hs_validate_julia_target("direct_maternal"),
    "direct_maternal"
  )
})

test_that("hs_effect_targets includes 'direct_maternal' for maternal_genetic", {
  targets <- hsquared:::hs_effect_targets("maternal_genetic")
  expect_true("direct_maternal" %in% targets)
  expect_true("two_effect" %in% targets)
})

test_that("hs_second_effect_target still suggests 'two_effect' as default", {
  # The default suggestion in the error message should still be two_effect, not
  # direct_maternal (which is the opt-in correlated path).
  expect_equal(
    hsquared:::hs_second_effect_target("maternal_genetic"),
    "two_effect"
  )
})

test_that("maternal_genetic with wrong target gives allowed-targets error", {
  ped <- make_dm_ped()
  dat <- make_dm_dat()
  expect_error(
    hsquared(
      y ~ animal(1 | id, pedigree = ped) + maternal_genetic(1 | dam),
      data = dat,
      family = stats::gaussian(),
      control = hs_control(
        engine = "julia",
        engine_control = list(target = "ai_reml")
      )
    ),
    # Should mention two_effect as the first allowed target
    "two_effect",
    fixed = TRUE
  )
})

test_that("target 'direct_maternal' without maternal_genetic term errors", {
  ped <- data.frame(
    id = c("a", "b", "c"), sire = c(NA, NA, "a"), dam = c(NA, NA, "b"),
    stringsAsFactors = FALSE
  )
  dat <- data.frame(y = c(1, 2), id = c("b", "c"), stringsAsFactors = FALSE)
  expect_error(
    hsquared(
      y ~ animal(1 | id, pedigree = ped),
      data = dat,
      family = stats::gaussian(),
      control = hs_control(
        engine = "julia",
        engine_control = list(target = "direct_maternal")
      )
    ),
    "maternal_genetic",
    fixed = TRUE
  )
})

# ---- (B) payload emitter ---------------------------------------------------- #

test_that("payload has two pedigree blocks for a maternal_genetic formula", {
  ped <- make_dm_ped()
  dat <- make_dm_dat()
  spec <- hsquared:::hs_build_model_spec(
    y ~ animal(1 | id, pedigree = ped) + maternal_genetic(1 | dam),
    data = dat,
    family = stats::gaussian(),
    REML = TRUE
  )
  payload <- hsquared:::hs_build_bridge_payload(spec)
  re <- payload$random_effects

  expect_length(re, 2L)
  expect_equal(re[[1L]]$name, "animal")
  expect_equal(re[[1L]]$type, "pedigree")
  expect_equal(re[[2L]]$name, "maternal")
  expect_equal(re[[2L]]$type, "pedigree")
})

test_that("block1 Z is Zd (record->animal) for maternal_genetic formula", {
  ped <- make_dm_ped()
  dat <- make_dm_dat()
  spec <- hsquared:::hs_build_model_spec(
    y ~ animal(1 | id, pedigree = ped) + maternal_genetic(1 | dam),
    data = dat,
    family = stats::gaussian(),
    REML = TRUE
  )
  payload <- hsquared:::hs_build_bridge_payload(spec)
  # Z (block1) is record->animal: dimensions n_records x n_pedigree_animals
  Zd <- as.matrix(payload$random_effects[[1L]]$Z)
  expect_equal(nrow(Zd), nrow(dat))
  expect_equal(ncol(Zd), length(payload$ids))
  # Each row sums to 1 (each record maps to exactly one animal)
  expect_equal(rowSums(Zd), rep(1, nrow(dat)))
})

test_that("block2 Z is Zm (record->dam) for maternal_genetic formula", {
  ped <- make_dm_ped()
  dat <- make_dm_dat()
  spec <- hsquared:::hs_build_model_spec(
    y ~ animal(1 | id, pedigree = ped) + maternal_genetic(1 | dam),
    data = dat,
    family = stats::gaussian(),
    REML = TRUE
  )
  payload <- hsquared:::hs_build_bridge_payload(spec)
  Zm <- as.matrix(payload$random_effects[[2L]]$Z)
  # Zm is record->dam: n_records x n_pedigree_animals (dam levels are pedigree ids)
  expect_equal(nrow(Zm), nrow(dat))
  # Each row sums to 1 (each record has exactly one dam)
  expect_equal(rowSums(Zm), rep(1, nrow(dat)))
})

test_that("block2 pedigree carries the same rows as block1 for maternal_genetic", {
  ped <- make_dm_ped()
  dat <- make_dm_dat()
  spec <- hsquared:::hs_build_model_spec(
    y ~ animal(1 | id, pedigree = ped) + maternal_genetic(1 | dam),
    data = dat,
    family = stats::gaussian(),
    REML = TRUE
  )
  payload <- hsquared:::hs_build_bridge_payload(spec)
  re <- payload$random_effects
  # Both blocks should share the pedigree structure (same founders / ordering)
  expect_equal(re[[1L]]$pedigree$id, re[[2L]]$pedigree$id)
  expect_equal(re[[1L]]$pedigree$sire, re[[2L]]$pedigree$sire)
  expect_equal(re[[1L]]$pedigree$dam,  re[[2L]]$pedigree$dam)
})

test_that("relmat_status is 'build_in_julia' for both blocks", {
  ped <- make_dm_ped()
  dat <- make_dm_dat()
  spec <- hsquared:::hs_build_model_spec(
    y ~ animal(1 | id, pedigree = ped) + maternal_genetic(1 | dam),
    data = dat,
    family = stats::gaussian(),
    REML = TRUE
  )
  payload <- hsquared:::hs_build_bridge_payload(spec)
  re <- payload$random_effects
  expect_equal(re[[1L]]$relmat_status, "build_in_julia")
  expect_equal(re[[2L]]$relmat_status, "build_in_julia")
})

# ---- (C) accessor fences ---------------------------------------------------- #

# Fabricate a minimal direct_maternal hsquared_fit to test extractor behaviour
# without a live Julia engine.
#
# Components are the source of truth. Default `r_am` is the defining identity
# `σ_dm / sqrt(σ²_ad · σ²_am)` so derived-estimand assertions can hold
# (A21 GAP / C4). Callers may still pass `r_am` as a free knob for plumbing
# or boundary-warning tests.
make_dm_fit <- function(r_am = NULL, converged = TRUE) {
  sigma_ad <- 0.30
  sigma_am <- 0.15
  sigma_dm <- -0.10
  if (is.null(r_am)) {
    r_am <- sigma_dm / sqrt(sigma_ad * sigma_am)
  }
  result <- list(
    variance_components = data.frame(
      component = c("direct", "maternal", "covariance", "residual"),
      estimate  = c(sigma_ad, sigma_am, sigma_dm, 0.55),
      stringsAsFactors = FALSE
    ),
    heritability = data.frame(
      term     = "direct",
      # sigma_P = sigma_ad + sigma_am + sigma_dm + sigma_e2
      #         = 0.30 + 0.15 + (-0.10) + 0.55 = 0.90  (Willham 1972)
      # Old denominator (0.30+0.15+0.55 = 1.00) was wrong; updated.
      estimate = 0.30 / (0.30 + 0.15 + (-0.10) + 0.55),
      stringsAsFactors = FALSE
    ),
    genetic_correlation = data.frame(
      term_1   = "direct",
      term_2   = "maternal",
      estimate = r_am,
      stringsAsFactors = FALSE
    ),
    direct_variance  = 0.30,
    partner_variance = 0.15,
    covariance       = -0.10,
    breeding_values  = data.frame(id = "g1", value = 0.1, stringsAsFactors = FALSE),
    random_effects   = list(
      animal   = data.frame(id = "g1", value = 0.1, stringsAsFactors = FALSE),
      maternal = data.frame(id = "g1", value = 0.0, stringsAsFactors = FALSE)
    ),
    maternal_effects  = data.frame(id = "g1", value = 0.0, stringsAsFactors = FALSE),
    fixed_effects     = c(`(Intercept)` = 1.8),
    loglik            = -10.5,
    nobs              = 8L,
    converged         = converged,
    diagnostics       = list(target = "direct_maternal")
  )
  structure(
    list(
      spec    = list(target = "direct_maternal", method = "REML",
                     family = list(family = "gaussian", link = "identity")),
      payload = NULL,
      result  = result,
      engine  = "HSquared.jl"
    ),
    class = "hsquared_fit"
  )
}

test_that("hs_fit_is_direct_maternal detects direct_maternal fit", {
  fit <- make_dm_fit()
  expect_true(hsquared:::hs_fit_is_direct_maternal(fit))
})

test_that("hs_fit_is_direct_maternal returns FALSE for a non-dm fit", {
  # Fabricate a simple v0.1 fit (no direct_variance field)
  fit_plain <- structure(
    list(spec = list(target = NULL), result = list(heritability = 0.3)),
    class = "hsquared_fit"
  )
  expect_false(hsquared:::hs_fit_is_direct_maternal(fit_plain))
})

test_that("direct_heritability() returns labelled data frame with fence", {
  fit <- make_dm_fit()
  out <- direct_heritability(fit)
  expect_s3_class(out, "data.frame")
  expect_equal(out$term, "direct")
  expect_true(is.numeric(out$estimate))
  expect_true(!is.null(attr(out, "interpretation")))
  # The fence must mention direct and the phenotypic variance sigma^2_P
  expect_match(attr(out, "interpretation"), "sigma^2_P", fixed = TRUE)
})

test_that("heritability() on direct_maternal fit warns and returns labelled triple", {
  fit <- make_dm_fit()
  expect_warning(
    out <- heritability(fit),
    "labelled",
    fixed = TRUE
  )
  expect_s3_class(out, "data.frame")
  # Must contain all four components
  expect_equal(
    out$component,
    c("h2_direct", "m2_maternal", "h2_total_willham", "r_am")
  )
  # h2_direct numerics — sigma_P = 0.30+0.15+(-0.10)+0.55 = 0.90
  expect_equal(out$estimate[out$component == "h2_direct"],
               0.30 / 0.90, tolerance = 1e-10)
  expect_equal(out$estimate[out$component == "m2_maternal"],
               0.15 / 0.90, tolerance = 1e-10)
  expect_equal(out$estimate[out$component == "h2_total_willham"],
               (0.30 + 1.5 * (-0.10) + 0.5 * 0.15) / 0.90, tolerance = 1e-10)
  expect_true(!is.null(attr(out, "interpretation")))
  expect_match(attr(out, "interpretation"), "Willham", fixed = TRUE)
})

test_that("total_heritability() returns Willham h2_T data frame", {
  fit <- make_dm_fit()
  out <- total_heritability(fit)
  expect_s3_class(out, "data.frame")
  expect_equal(out$term, "total_willham")
  # sigma_P = 0.30 + 0.15 + (-0.10) + 0.55 = 0.90
  # h2_T = (0.30 + 1.5*(-0.10) + 0.5*0.15) / 0.90
  #       = (0.30 - 0.15 + 0.075) / 0.90 = 0.225 / 0.90 = 0.25
  expect_equal(out$estimate,
               (0.30 + 1.5 * (-0.10) + 0.5 * 0.15) / 0.90,
               tolerance = 1e-10)
  expect_true(!is.null(attr(out, "interpretation")))
  expect_match(attr(out, "interpretation"), "Willham", fixed = TRUE)
  expect_match(attr(out, "interpretation"), "MASS SELECTION", fixed = TRUE)
})

test_that("total_heritability() can be lower than direct_heritability() when r_am < 0", {
  fit <- make_dm_fit(r_am = -0.4)
  h2_T <- total_heritability(fit)$estimate
  h2_d <- direct_heritability(fit)$estimate
  # With sigma_dm = -0.10, h2_T = 0.25 < h2_d = 0.30/0.90 ~ 0.333
  expect_lt(h2_T, h2_d)
})

test_that("total_heritability() errors on non-dm fit", {
  fit_plain <- structure(
    list(spec = list(target = NULL), result = list()),
    class = "hsquared_fit"
  )
  expect_error(total_heritability(fit_plain), "direct_maternal", fixed = TRUE)
})

test_that("total_heritability() errors on non-hsquared_fit object", {
  expect_error(total_heritability(list(x = 1)), "direct_maternal", fixed = TRUE)
})

test_that("genetic_correlation() on direct_maternal fit returns r_am data frame", {
  fit <- make_dm_fit(r_am = -0.4)
  out <- genetic_correlation(fit)
  expect_s3_class(out, "data.frame")
  expect_equal(out$term_1, "direct")
  expect_equal(out$term_2, "maternal")
  expect_equal(out$estimate, -0.4)
})

# A21 GAP: r_am is a covered derived estimand. The gate requires an identity
# test against the defining function of the components, not a plumbing
# round-trip. The default fixture now derives r_am from those components.
test_that("r_am identity: genetic_correlation equals cov / sqrt(var_d * var_m)", {
  fit <- make_dm_fit()
  sigma_ad <- 0.30
  sigma_am <- 0.15
  sigma_dm <- -0.10
  expected <- sigma_dm / sqrt(sigma_ad * sigma_am)
  expect_equal(genetic_correlation(fit)$estimate, expected, tolerance = 1e-12)
  h <- suppressWarnings(heritability(fit))
  expect_equal(h$estimate[h$component == "r_am"], expected, tolerance = 1e-12)
})

# A21 C6: direct_heritability() reads result$heritability; heritability()
# recomputes h2_direct from variance_components. They must stay identical.
test_that("direct_heritability() matches heritability() h2_direct", {
  fit <- make_dm_fit()
  h <- suppressWarnings(heritability(fit))
  expect_equal(
    direct_heritability(fit)$estimate,
    h$estimate[h$component == "h2_direct"],
    tolerance = 1e-12
  )
})

# A21 C6 fence prose. docs/ is .Rbuildignore'd, so this is source-tree-only
# (skipped under R CMD check), matching the RR grammar-surface guard.
test_that("capability-status maternal fences distinguish the two m2 denominators", {
  path <- file.path(
    testthat::test_path("..", ".."),
    "docs",
    "design",
    "capability-status.md"
  )
  skip_if_not(file.exists(path), "capability-status.md is not in the build tarball")
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  expect_true(
    grepl("maternal_proportion()", text, fixed = TRUE),
    info = "two-effect fence names maternal_proportion()"
  )
  expect_true(
    grepl("no covariance", text, fixed = TRUE),
    info = "two-effect fence names the no-covariance m2 denominator"
  )
  expect_true(
    grepl("omits the covariance", text, fixed = TRUE),
    info = "direct-maternal fence names the other m2 denominator"
  )
  expect_true(
    grepl("m2_maternal", text, fixed = TRUE),
    info = "direct-maternal fence names m2_maternal"
  )
})

test_that("genetic_correlation() warns when |r_am| >= 0.99", {
  fit <- make_dm_fit(r_am = 0.999)
  expect_warning(
    genetic_correlation(fit),
    "boundary",
    fixed = TRUE
  )
})

test_that("genetic_correlation() does NOT warn for a negative r_am within bounds", {
  fit <- make_dm_fit(r_am = -0.4)
  expect_no_warning(genetic_correlation(fit))
})

test_that("direct_variance() returns sigma_ad as a single numeric", {
  fit <- make_dm_fit()
  expect_equal(direct_variance(fit), 0.30)
})

test_that("partner_variance() returns sigma_am as a single numeric", {
  fit <- make_dm_fit()
  expect_equal(partner_variance(fit), 0.15)
})

test_that("direct_maternal_covariance() returns sigma_dm (may be negative)", {
  fit <- make_dm_fit()
  expect_equal(direct_maternal_covariance(fit), -0.10)
})

test_that("direct_heritability() errors on non-dm fit", {
  fit_plain <- structure(
    list(
      spec   = list(target = NULL),
      result = list(heritability = 0.3)
    ),
    class = "hsquared_fit"
  )
  expect_error(
    direct_heritability(fit_plain),
    "direct_maternal",
    fixed = TRUE
  )
})

test_that("direct_heritability() errors on non-hsquared_fit object", {
  expect_error(
    direct_heritability(list(x = 1)),
    "direct_maternal",
    fixed = TRUE
  )
})

test_that("direct_variance() errors on non-dm fit", {
  fit_plain <- structure(
    list(spec = list(target = NULL), result = list()),
    class = "hsquared_fit"
  )
  expect_error(direct_variance(fit_plain), "direct_maternal", fixed = TRUE)
})

test_that("partner_variance() errors on non-dm fit", {
  fit_plain <- structure(
    list(spec = list(target = NULL), result = list()),
    class = "hsquared_fit"
  )
  expect_error(partner_variance(fit_plain), "direct_maternal", fixed = TRUE)
})

test_that("direct_maternal_covariance() errors on non-dm fit", {
  fit_plain <- structure(
    list(spec = list(target = NULL), result = list()),
    class = "hsquared_fit"
  )
  expect_error(
    direct_maternal_covariance(fit_plain),
    "direct_maternal",
    fixed = TRUE
  )
})

test_that("maternal_effects() is accessible from a direct_maternal fit", {
  fit <- make_dm_fit()
  out <- maternal_effects(fit)
  expect_s3_class(out, "data.frame")
})

# ---- (D) live parity test (skip-guarded) ------------------------------------ #

# This test exercises the full R -> Julia bridge for the direct-maternal model.
# It requires:
#   1. Julia installed
#   2. JuliaCall R package installed
#   3. A local HSquared.jl checkout (via HSQUARED_JULIA_PROJECT env var or
#      a detectable local path)
#   4. HSquared.jl version that includes fit_direct_maternal_reml
# The test is skipped unless both JuliaCall and the bridge are available.

skip_if_no_julia_bridge <- function() {
  if (!requireNamespace("JuliaCall", quietly = TRUE)) {
    skip("JuliaCall not installed — skipping live parity test")
  }
  project <- hsquared:::hs_default_julia_project()
  if (!hsquared:::hs_julia_bridge_available(project)) {
    skip("HSquared.jl Julia project not found — skipping live parity test")
  }
}

test_that("live R<->engine parity: direct_maternal fit returns converged result", {
  hs_skip_live_julia()
  skip_if_no_julia_bridge()

  ped <- make_dm_ped()
  dat <- make_dm_dat()

  fit_dm <- hsquared(
    y ~ animal(1 | id, pedigree = ped) + maternal_genetic(1 | dam),
    data   = dat,
    family = stats::gaussian(),
    REML   = TRUE,
    control = hs_control(
      engine = "julia",
      engine_control = list(target = "direct_maternal")
    )
  )

  # Shape checks: the fit must be an hsquared_fit with the right structure
  expect_s3_class(fit_dm, "hsquared_fit")
  expect_equal(fit_dm$spec$target, "direct_maternal")

  vc <- variance_components(fit_dm)
  expect_s3_class(vc, "data.frame")
  expect_setequal(
    vc$component,
    c("direct", "maternal", "covariance", "residual")
  )
  # All variances must be finite; covariance may be negative
  expect_true(all(is.finite(vc$estimate)))
  # direct and maternal variances must be non-negative
  expect_gte(vc$estimate[vc$component == "direct"], 0)
  expect_gte(vc$estimate[vc$component == "maternal"], 0)

  # Genetic correlation in [-1, 1], and the A21 GAP identity against
  # the returned components (not a range-only check).
  r_am <- suppressWarnings(genetic_correlation(fit_dm)$estimate)
  expect_gte(r_am, -1)
  expect_lte(r_am, 1)
  sigma_ad_live <- vc$estimate[vc$component == "direct"]
  sigma_am_live <- vc$estimate[vc$component == "maternal"]
  sigma_dm_live <- vc$estimate[vc$component == "covariance"]
  if (sigma_ad_live > 0 && sigma_am_live > 0) {
    expect_equal(
      r_am,
      sigma_dm_live / sqrt(sigma_ad_live * sigma_am_live),
      tolerance = 1e-8
    )
  }

  # direct_heritability in [0, 1]
  h2d <- direct_heritability(fit_dm)$estimate
  expect_gte(h2d, 0)
  expect_lte(h2d, 1)

  # Scalar extractors return finite numerics
  expect_true(is.finite(direct_variance(fit_dm)))
  expect_true(is.finite(partner_variance(fit_dm)))
  expect_true(is.finite(direct_maternal_covariance(fit_dm)))

  # BLUPs: direct animal effects
  bv <- breeding_values(fit_dm)
  expect_s3_class(bv, "data.frame")
  expect_true("id" %in% names(bv))

  # Maternal EBVs
  mat_eff <- maternal_effects(fit_dm)
  expect_s3_class(mat_eff, "data.frame")

  # Parity: direct Julia call vs R bridge
  # Build the engine inputs directly and call fit_direct_maternal_reml,
  # then compare sigma_ad to within 1e-6.
  spec <- hsquared:::hs_build_model_spec(
    y ~ animal(1 | id, pedigree = ped) + maternal_genetic(1 | dam),
    data   = dat,
    family = stats::gaussian(),
    REML   = TRUE
  )
  payload <- hsquared:::hs_build_bridge_payload(spec)
  Zd_mat <- as.matrix(payload$random_effects[[1L]]$Z)
  Zm_mat <- as.matrix(payload$random_effects[[2L]]$Z)

  hsquared:::hs_julia_setup(hsquared:::hs_default_julia_project())
  JuliaCall::julia_assign("hsq_parity_y", payload$y)
  JuliaCall::julia_assign("hsq_parity_X", payload$X)
  JuliaCall::julia_assign("hsq_parity_Zd", Zd_mat)
  JuliaCall::julia_assign("hsq_parity_Zm", Zm_mat)
  JuliaCall::julia_assign(
    "hsq_parity_ped_id",
    as.character(payload$pedigree$id)
  )
  JuliaCall::julia_assign(
    "hsq_parity_ped_sire",
    hsquared:::hs_parent_for_julia(payload$pedigree$sire)
  )
  JuliaCall::julia_assign(
    "hsq_parity_ped_dam",
    hsquared:::hs_parent_for_julia(payload$pedigree$dam)
  )
  JuliaCall::julia_command(paste(
    "let ped = HSquared.normalize_pedigree(",
    "hsq_parity_ped_id, hsq_parity_ped_sire, hsq_parity_ped_dam);",
    "global hsq_parity_Ainv = Matrix(HSquared.pedigree_inverse(ped)); end"
  ))
  JuliaCall::julia_command(paste(
    "global hsq_parity_ref = HSquared.fit_direct_maternal_reml(",
    "hsq_parity_y, hsq_parity_X, hsq_parity_Zd, hsq_parity_Zm,",
    "hsq_parity_Ainv);"
  ))
  ref_sigma_ad <- JuliaCall::julia_eval(
    "Float64(hsq_parity_ref.variance_components.sigma_ad)"
  )
  ref_sigma_am <- JuliaCall::julia_eval(
    "Float64(hsq_parity_ref.variance_components.sigma_am)"
  )

  # Bridge result must match the direct Julia call to within 1e-6
  expect_equal(direct_variance(fit_dm), ref_sigma_ad, tolerance = 1e-6)
  expect_equal(partner_variance(fit_dm), ref_sigma_am, tolerance = 1e-6)
})
