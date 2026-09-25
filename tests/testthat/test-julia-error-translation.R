# Translate raw Julia errors at the R-Julia bridge into structured
# `hsquared_error` conditions (hsquared#214, #217): the dense-cells cap
# (`max_dense_cells`, forwarded through `engine_control`) and the
# `single_step()` H^-1 construction path's `ridge`/`blend_weight` boundary
# both used to surface as an untranslated Julia stack trace. The unit test
# (a) needs no Julia; (b)/(c)/(e) are live-gated on the local HSquared.jl
# bridge (see repository AGENTS.md).

test_that("hs_julia_fit() translates an error raised inside expr", {
  err <- tryCatch(
    hsquared:::hs_julia_fit(stop("boom"), hint = "try this instead"),
    error = function(e) e
  )
  expect_equal(
    class(err),
    c("hsquared_julia_error", "hsquared_error", "error", "condition")
  )
  expect_match(conditionMessage(err), "boom", fixed = TRUE)
  expect_match(conditionMessage(err), "try this instead", fixed = TRUE)
})

test_that("hs_julia_fit() translates an error with no hint supplied", {
  err <- tryCatch(
    hsquared:::hs_julia_fit(stop("boom")),
    error = function(e) e
  )
  expect_equal(
    class(err),
    c("hsquared_julia_error", "hsquared_error", "error", "condition")
  )
  expect_match(conditionMessage(err), "boom", fixed = TRUE)
})

test_that("hs_julia_fit() passes a non-error result through unchanged", {
  expect_equal(hsquared:::hs_julia_fit(1 + 1), 2)
  expect_equal(hsquared:::hs_julia_fit("value", hint = "unused"), "value")
})

test_that("hs_control() rejects a non-positive max_dense_cells without Julia", {
  expect_error(
    hs_control(engine_control = list(max_dense_cells = -1)),
    "max_dense_cells",
    fixed = TRUE
  )
})

test_that("hs_control() accepts a positive integer max_dense_cells", {
  ctl <- hs_control(engine_control = list(max_dense_cells = 50))
  expect_equal(ctl$engine_control$max_dense_cells, 50)
})

test_that("the lever-free scale hint names no control (non-guarded routes)", {
  err <- tryCatch(
    hsquared:::hs_julia_fit(stop("boom"), hint = hsquared:::hs_dense_scale_hint),
    error = function(e) e
  )
  expect_false(grepl("max_dense_cells", conditionMessage(err), fixed = TRUE))
})

test_that("the dense route hint names max_dense_cells (guarded routes)", {
  err <- tryCatch(
    hsquared:::hs_julia_fit(stop("boom"), hint = hsquared:::hs_dense_route_hint),
    error = function(e) e
  )
  expect_match(conditionMessage(err), "max_dense_cells", fixed = TRUE)
})

test_that("each hs_fit_julia_*_payload() call site uses the hint matching what it forwards", {
  # Only the two routes that actually forward `max_dense_cells` to a guarded
  # engine entry point (`fit_variance_components`, `fit_repeatability_reml`)
  # may name it; naming it elsewhere points the user at a control that
  # target does not accept (hsquared#212's defect class, R3-3).
  body_text_of <- function(fn) {
    paste(deparse(body(getFromNamespace(fn, "hsquared"))), collapse = "\n")
  }

  # Each guarded route names ITS OWN dense hint. The repeatability route has a
  # named sparse escape (`scale_method = "auto"`), so it carries a specific
  # hint rather than the generic "use a sparse route" wording; the invariant
  # that matters is that whatever hint it uses still names `max_dense_cells`,
  # which is asserted on the constants directly below.
  guarded <- c(
    hs_fit_julia_payload = "hs_dense_route_hint",
    hs_fit_julia_repeatability_payload = "hs_repeatability_dense_route_hint"
  )
  for (fn in names(guarded)) {
    txt <- body_text_of(fn)
    expect_match(txt, guarded[[fn]], fixed = TRUE, info = fn)
    expect_false(grepl("hs_dense_scale_hint", txt, fixed = TRUE), info = fn)
  }

  # The actual user-facing contract: every dense hint names the control that
  # lifts the ceiling, and the repeatability one also names the sparse route.
  expect_match(
    hsquared:::hs_dense_route_hint, "max_dense_cells",
    fixed = TRUE
  )
  expect_match(
    hsquared:::hs_repeatability_dense_route_hint, "max_dense_cells",
    fixed = TRUE
  )
  expect_match(
    hsquared:::hs_repeatability_dense_route_hint, "scale_method",
    fixed = TRUE
  )

  unguarded_dense <- c(
    "hs_fit_julia_henderson_mme_payload",
    "hs_fit_julia_metafounder_payload",
    "hs_fit_julia_sparse_reml_payload",
    "hs_fit_julia_ai_reml_payload",
    "hs_fit_julia_nongaussian_payload",
    "hs_fit_julia_two_effect_payload",
    "hs_fit_julia_direct_maternal_payload",
    "hs_fit_julia_n_effect_payload",
    "hs_fit_julia_multivariate_payload",
    "hs_fit_julia_random_regression_payload",
    "hs_fit_julia_genomic_payload",
    "hs_fit_julia_snp_blup_payload",
    "hs_fit_julia_snp_blup_reml_payload"
  )
  expect_length(unguarded_dense, 13L)
  for (fn in unguarded_dense) {
    txt <- body_text_of(fn)
    expect_match(txt, "hs_dense_scale_hint", fixed = TRUE, info = fn)
    expect_false(grepl("hs_dense_route_hint", txt, fixed = TRUE), info = fn)
  }

  single_step <- c(
    "hs_fit_julia_single_step_construct_payload",
    "hs_fit_julia_metafounder_single_step_payload"
  )
  for (fn in single_step) {
    txt <- body_text_of(fn)
    expect_match(txt, "hs_single_step_ridge_hint", fixed = TRUE, info = fn)
  }
})

# --- live tests (skip-guarded on the local HSquared.jl bridge) --------------

test_that("single_step_construct at ridge = 0 on an ill-conditioned genotyped block raises a translated hsquared_error [live]", {
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for this live test."
  )
  hsquared:::hs_julia_setup(hsquared:::hs_default_julia_project())

  ped <- data.frame(
    id = c("c1", "s", "d", "c2", "c3"),
    sire = c("s", NA, NA, "s", "s"),
    dam = c("d", NA, NA, "d", "d"),
    stringsAsFactors = FALSE
  )
  dat <- data.frame(
    y = c(1.1, 2.2, 3.0, 2.5, 1.8),
    id = ped$id,
    stringsAsFactors = FALSE
  )
  # c1 and c2 carry an identical genotype row: only 2 informative markers span
  # 3 genotyped animals, so the sample genomic relationship block is
  # rank-deficient (not positive definite) at the documented default
  # ridge = 0 (issue #214 comments, wave 3).
  m <- matrix(
    c(0, 1, 0, 1, 2, 0),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("c1", "c2", "c3"), c("m1", "m2"))
  )

  err <- tryCatch(
    hsquared(
      y ~ single_step(1 | id, pedigree = ped, markers = m),
      data = dat,
      control = hs_control(
        engine = "julia",
        engine_control = list(target = "single_step_construct")
      )
    ),
    error = function(e) e
  )
  expect_s3_class(err, "hsquared_error")
  expect_match(conditionMessage(err), "ridge|blend_weight", perl = TRUE)
})

test_that("target = repeatability honours max_dense_cells [live]", {
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for this live test."
  )
  hsquared:::hs_julia_setup(hsquared:::hs_default_julia_project())

  ped <- data.frame(
    id = c("a", "b", "c", "d", "e"),
    sire = c(NA, NA, NA, "a", "a"),
    dam = c(NA, NA, NA, "b", "c")
  )
  set.seed(11)
  ids <- rep(c("a", "b", "c", "d", "e"), each = 3)
  pe <- stats::setNames(stats::rnorm(5, 0, 1.0), c("a", "b", "c", "d", "e"))
  dat <- data.frame(
    y = 2 + pe[ids] + stats::rnorm(15, 0, 0.7),
    id = ids
  )

  err <- tryCatch(
    hsquared(
      y ~ animal(1 | id, pedigree = ped) + permanent(1 | id),
      data = dat,
      control = hs_control(
        engine = "julia",
        engine_control = list(target = "repeatability", max_dense_cells = 50)
      )
    ),
    error = function(e) e
  )
  expect_s3_class(err, "hsquared_error")
  expect_match(conditionMessage(err), "max_dense_cells = 50", fixed = TRUE)

  fit <- hsquared(
    y ~ animal(1 | id, pedigree = ped) + permanent(1 | id),
    data = dat,
    control = hs_control(
      engine = "julia",
      engine_control = list(target = "repeatability", max_dense_cells = 1e7)
    )
  )
  expect_s3_class(fit, "hsquared_fit")
  expect_equal(fit$spec$target, "repeatability")
})

test_that("the default animal() route honours max_dense_cells [live]", {
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for this live test."
  )
  hsquared:::hs_julia_setup(hsquared:::hs_default_julia_project())

  ped <- hs_sim_pedigree(n_founder = 10, n_per_gen = 5, n_gen = 2, seed = 1)
  dat <- hs_sim_genedrop_phenotypes(ped, sigma_a2 = 0.4, sigma_e2 = 0.6, seed = 1)
  expect_equal(nrow(ped), 20L)

  err <- tryCatch(
    hsquared(
      y ~ animal(1 | id, pedigree = ped),
      data = dat,
      control = hs_control(
        engine = "julia",
        engine_control = list(max_dense_cells = 50)
      )
    ),
    error = function(e) e
  )
  expect_s3_class(err, "hsquared_error")
  expect_match(conditionMessage(err), "max_dense_cells = 50", fixed = TRUE)
})
