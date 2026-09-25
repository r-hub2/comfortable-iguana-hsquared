# hsquared#212: `engine_control$initial`/`$iterations` were silently
# discarded (no error, no warning) on several Julia payload builders:
#   * hs_fit_julia_n_effect_payload()             (target = "multi_effect")
#   * hs_fit_julia_direct_maternal_payload()       (target = "direct_maternal")
#   * hs_fit_julia_single_step_construct_payload() (`iterations` only)
#   * hs_fit_julia_metafounder_single_step_payload() (`iterations` only)
#   * the multi_effect_ratio_interval() refit (a separate code path from the
#     main multi_effect fit above)
#
# (a) Julia-gated [live]: for each affected target, a default-control fit vs.
#     an extreme-`initial` + `iterations = 1` probe fit must produce a
#     DIFFERENT variance-component estimate, or the probe must fail to
#     converge -- proving `initial`/`iterations` now reach the Julia call.
# (b) NOT Julia-gated: supplying an `engine_control` key a target does not
#     honour errors with class `hsquared_unsupported_syntax`, naming the key
#     and the target (`hs_engine_control_forwarding()`).

hs212_max_rel_diff <- function(vc_base, vc_probe) {
  max(abs(vc_probe - vc_base) / pmax(abs(vc_base), 1e-8))
}

# ---- (a) live fixtures ------------------------------------------------------

hs212_dm_ped <- function() {
  data.frame(
    id   = c("g1", "g2", "g3", "g4", "g5", "g6"),
    sire = c(NA,   NA,   "g1", "g1", "g2", "g2"),
    dam  = c(NA,   NA,   "g2", "g2", "g1", "g1"),
    stringsAsFactors = FALSE
  )
}

hs212_dm_dat <- function() {
  data.frame(
    y   = c(1.2, 0.8, 2.1, 1.9, 1.5, 2.3, 1.7, 0.9),
    id  = c("g3", "g4", "g5", "g6", "g3", "g4", "g5", "g6"),
    dam = c("g2", "g2", "g1", "g1", "g2", "g2", "g1", "g1"),
    stringsAsFactors = FALSE
  )
}

hs212_me_ped <- function() {
  data.frame(
    id   = c("a", "b", "c", "d", "e", "f", "g", "h"),
    sire = c(NA,  NA,  NA,  NA,  "a", "a", "c", "c"),
    dam  = c(NA,  NA,  NA,  NA,  "b", "b", "d", "d"),
    stringsAsFactors = FALSE
  )
}

hs212_me_dat <- function(ped, seed = 212) {
  set.seed(seed)
  nest <- c("nst1", "nst1", "nst2", "nst2", "nst3", "nst3", "nst4", "nst4")
  year <- c("y1", "y2", "y1", "y2", "y1", "y2", "y1", "y2")
  nest_e <- stats::setNames(
    stats::rnorm(4, 0, 0.6),
    c("nst1", "nst2", "nst3", "nst4")
  )
  year_e <- stats::setNames(stats::rnorm(2, 0, 0.5), c("y1", "y2"))
  data.frame(
    y = 3 + nest_e[nest] + year_e[year] + stats::rnorm(8, 0, 0.7),
    id = ped$id,
    nest = nest,
    year = year,
    stringsAsFactors = FALSE
  )
}

hs212_ss_ped <- function(seed = 212) {
  hs_sim_pedigree(n_founder = 10, n_per_gen = 10, n_gen = 1, seed = seed)
}

hs212_ss_markers <- function(ped, seed = 212) {
  geno <- utils::tail(ped$id, 10L)
  set.seed(seed)
  matrix(
    stats::rbinom(length(geno) * 30L, 2L, 0.3),
    nrow = length(geno),
    dimnames = list(geno, paste0("snp", seq_len(30L)))
  )
}

hs212_mf_group <- function(ped) {
  out <- rep("", nrow(ped))
  out[is.na(ped$sire) | is.na(ped$dam)] <- "base"
  names(out) <- as.character(ped$id)
  out
}

hs212_skip_unless_bridge <- function() {
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for a live fit."
  )
}

# ---- (a) live: initial/iterations must move the fit ------------------------

test_that("direct_maternal: extreme initial + iterations=1 moves the fit [live, #212]", {
  hs212_skip_unless_bridge()

  ped <- hs212_dm_ped()
  dat <- hs212_dm_dat()

  baseline <- hsquared(
    y ~ animal(1 | id, pedigree = ped) + maternal_genetic(1 | dam),
    data = dat,
    control = hs_control(
      engine = "julia",
      engine_control = list(target = "direct_maternal")
    )
  )
  probe <- hsquared(
    y ~ animal(1 | id, pedigree = ped) + maternal_genetic(1 | dam),
    data = dat,
    control = hs_control(
      engine = "julia",
      engine_control = list(
        target = "direct_maternal",
        initial = list(G_dm = diag(2) * 500, sigma_e2 = 500),
        iterations = 1L
      )
    )
  )

  vc_base <- variance_components(baseline)$estimate
  vc_probe <- variance_components(probe)$estimate
  max_rel_diff <- hs212_max_rel_diff(vc_base, vc_probe)
  expect_true(
    max_rel_diff > 0.01 || !isTRUE(probe$result$converged),
    label = sprintf("max_rel_diff = %.3e", max_rel_diff)
  )
})

test_that("multi_effect: extreme initial + iterations=1 moves the fit [live, #212]", {
  hs212_skip_unless_bridge()

  ped <- hs212_me_ped()
  dat <- hs212_me_dat(ped)

  baseline <- hsquared(
    y ~ animal(1 | id, pedigree = ped) + (1 | nest) + (1 | year),
    data = dat,
    control = hs_control(
      engine = "julia",
      engine_control = list(target = "multi_effect")
    )
  )
  probe <- hsquared(
    y ~ animal(1 | id, pedigree = ped) + (1 | nest) + (1 | year),
    data = dat,
    control = hs_control(
      engine = "julia",
      engine_control = list(
        target = "multi_effect",
        # K = 3 blocks (animal, nest, year) + residual -> length 4.
        initial = c(500, 500, 500, 500),
        iterations = 1L
      )
    )
  )

  vc_base <- variance_components(baseline)$estimate
  vc_probe <- variance_components(probe)$estimate
  max_rel_diff <- hs212_max_rel_diff(vc_base, vc_probe)
  expect_true(
    max_rel_diff > 0.01 || !isTRUE(probe$result$converged),
    label = sprintf("max_rel_diff = %.3e", max_rel_diff)
  )
})

test_that("single_step_construct: extreme initial + iterations=1 moves the fit [live, #212]", {
  hs212_skip_unless_bridge()

  ped <- hs212_ss_ped(seed = 212)
  dat <- hs_sim_genedrop_phenotypes(ped, sigma_a2 = 0.4, sigma_e2 = 0.6, seed = 212)
  markers <- hs212_ss_markers(ped, seed = 212)

  baseline <- hsquared(
    y ~ single_step(1 | id, pedigree = ped, markers = markers, ridge = 0.05),
    data = dat,
    control = hs_control(
      engine = "julia",
      engine_control = list(target = "single_step_construct")
    )
  )
  probe <- hsquared(
    y ~ single_step(1 | id, pedigree = ped, markers = markers, ridge = 0.05),
    data = dat,
    control = hs_control(
      engine = "julia",
      engine_control = list(
        target = "single_step_construct",
        initial = c(sigma_a2 = 1e6, sigma_e2 = 1e6),
        iterations = 1L
      )
    )
  )

  vc_base <- variance_components(baseline)$estimate
  vc_probe <- variance_components(probe)$estimate
  max_rel_diff <- hs212_max_rel_diff(vc_base, vc_probe)
  expect_true(
    max_rel_diff > 0.01 || !isTRUE(probe$result$converged),
    label = sprintf("max_rel_diff = %.3e", max_rel_diff)
  )
})

test_that("metafounder_single_step: extreme initial + iterations=1 moves the fit [live, #212]", {
  hs212_skip_unless_bridge()

  ped <- hs212_ss_ped(seed = 213)
  dat <- hs_sim_genedrop_phenotypes(ped, sigma_a2 = 0.4, sigma_e2 = 0.6, seed = 213)
  markers <- hs212_ss_markers(ped, seed = 213)
  group <- hs212_mf_group(ped)
  Gamma <- matrix(0, nrow = 1, dimnames = list("base", "base"))

  baseline <- hsquared(
    y ~ single_step(
      1 | id,
      pedigree = ped,
      markers = markers,
      group = group,
      Gamma = Gamma,
      ridge = 0.05
    ),
    data = dat,
    control = hs_control(
      engine = "julia",
      engine_control = list(target = "metafounder_single_step")
    )
  )
  probe <- hsquared(
    y ~ single_step(
      1 | id,
      pedigree = ped,
      markers = markers,
      group = group,
      Gamma = Gamma,
      ridge = 0.05
    ),
    data = dat,
    control = hs_control(
      engine = "julia",
      engine_control = list(
        target = "metafounder_single_step",
        initial = c(sigma_a2 = 1e6, sigma_e2 = 1e6),
        iterations = 1L
      )
    )
  )

  vc_base <- variance_components(baseline)$estimate
  vc_probe <- variance_components(probe)$estimate
  max_rel_diff <- hs212_max_rel_diff(vc_base, vc_probe)
  expect_true(
    max_rel_diff > 0.01 || !isTRUE(probe$result$converged),
    label = sprintf("max_rel_diff = %.3e", max_rel_diff)
  )
})

# ---- (b) not Julia-gated: unsupported engine_control keys error ------------
#
# `hs_engine_control_forwarding()` runs before any target-specific formula
# check, so a minimal animal-only formula suffices regardless of `target`.

hs212_ped <- function() {
  data.frame(
    id   = c("a", "b", "c", "d"),
    sire = c(NA,  NA,  "a", "a"),
    dam  = c(NA,  NA,  "b", "b"),
    stringsAsFactors = FALSE
  )
}

hs212_dat <- function() {
  data.frame(
    y  = c(1, 2, 3, 4),
    id = c("a", "b", "c", "d"),
    stringsAsFactors = FALSE
  )
}

test_that("direct_maternal does not honour em_warmup (unsupported key errors)", {
  ped <- hs212_ped()
  dat <- hs212_dat()

  cnd <- expect_error(
    hsquared(
      y ~ animal(1 | id, pedigree = ped),
      data = dat,
      control = hs_control(
        engine = "julia",
        engine_control = list(target = "direct_maternal", em_warmup = 3L)
      )
    ),
    class = "hsquared_unsupported_syntax"
  )
  expect_match(conditionMessage(cnd), "em_warmup", fixed = TRUE)
  expect_match(conditionMessage(cnd), "direct_maternal", fixed = TRUE)
})

test_that("multi_effect does not honour variance_components (unsupported key errors)", {
  ped <- hs212_ped()
  dat <- hs212_dat()

  cnd <- expect_error(
    hsquared(
      y ~ animal(1 | id, pedigree = ped),
      data = dat,
      control = hs_control(
        engine = "julia",
        engine_control = list(
          target = "multi_effect",
          variance_components = c(sigma_a2 = 1, sigma_e2 = 1)
        )
      )
    ),
    class = "hsquared_unsupported_syntax"
  )
  expect_match(conditionMessage(cnd), "variance_components", fixed = TRUE)
  expect_match(conditionMessage(cnd), "multi_effect", fixed = TRUE)
})

test_that("single_step_construct does not honour em_warmup (unsupported key errors)", {
  ped <- hs212_ped()
  dat <- hs212_dat()

  cnd <- expect_error(
    hsquared(
      y ~ animal(1 | id, pedigree = ped),
      data = dat,
      control = hs_control(
        engine = "julia",
        engine_control = list(target = "single_step_construct", em_warmup = 1L)
      )
    ),
    class = "hsquared_unsupported_syntax"
  )
  expect_match(conditionMessage(cnd), "em_warmup", fixed = TRUE)
  expect_match(conditionMessage(cnd), "single_step_construct", fixed = TRUE)
})

test_that("metafounder_single_step does not honour marginal (unsupported key errors)", {
  ped <- hs212_ped()
  dat <- hs212_dat()

  cnd <- expect_error(
    hsquared(
      y ~ animal(1 | id, pedigree = ped),
      data = dat,
      control = hs_control(
        engine = "julia",
        engine_control = list(target = "metafounder_single_step", marginal = "laplace")
      )
    ),
    class = "hsquared_unsupported_syntax"
  )
  expect_match(conditionMessage(cnd), "marginal", fixed = TRUE)
  expect_match(conditionMessage(cnd), "metafounder_single_step", fixed = TRUE)
})

test_that("henderson_mme does not honour initial (unsupported key errors)", {
  ped <- hs212_ped()
  dat <- hs212_dat()

  cnd <- expect_error(
    hsquared(
      y ~ animal(1 | id, pedigree = ped),
      data = dat,
      control = hs_control(
        engine = "julia",
        engine_control = list(
          target = "henderson_mme",
          variance_components = c(sigma_a2 = 1, sigma_e2 = 1),
          initial = c(sigma_a2 = 1, sigma_e2 = 1)
        )
      )
    ),
    class = "hsquared_unsupported_syntax"
  )
  expect_match(conditionMessage(cnd), "initial", fixed = TRUE)
  expect_match(conditionMessage(cnd), "henderson_mme", fixed = TRUE)
})

test_that("random_regression does not honour initial (unsupported key errors)", {
  ped <- hs212_ped()
  dat <- hs212_dat()

  cnd <- expect_error(
    hsquared(
      y ~ animal(1 | id, pedigree = ped),
      data = dat,
      control = hs_control(
        engine = "julia",
        engine_control = list(target = "random_regression", initial = c(sigma_a2 = 1, sigma_e2 = 1))
      )
    ),
    class = "hsquared_unsupported_syntax"
  )
  expect_match(conditionMessage(cnd), "initial", fixed = TRUE)
  expect_match(conditionMessage(cnd), "random_regression", fixed = TRUE)
})

test_that("the default fit_animal_model target does not honour iterations", {
  ped <- hs212_ped()
  dat <- hs212_dat()

  cnd <- expect_error(
    hsquared(
      y ~ animal(1 | id, pedigree = ped),
      data = dat,
      control = hs_control(
        engine = "julia",
        engine_control = list(iterations = 50L)
      )
    ),
    class = "hsquared_unsupported_syntax"
  )
  expect_match(conditionMessage(cnd), "iterations", fixed = TRUE)
  expect_match(conditionMessage(cnd), "fit_animal_model", fixed = TRUE)
})

# ---- (b) unit coverage of hs_engine_control_forwarding() itself ------------

test_that("hs_engine_control_forwarding returns the honoured subset without error", {
  control <- hs_control(
    engine = "julia",
    engine_control = list(
      target = "ai_reml",
      initial = c(sigma_a2 = 1, sigma_e2 = 1),
      iterations = 50L,
      em_warmup = 2L
    )
  )
  out <- hsquared:::hs_engine_control_forwarding(control, "ai_reml")
  expect_equal(sort(names(out)), sort(c("em_warmup", "initial", "iterations")))
})

test_that("hs_engine_control_forwarding is a no-op for an empty engine_control", {
  control <- hs_control(engine = "julia")
  expect_equal(
    hsquared:::hs_engine_control_forwarding(control, "direct_maternal"),
    list()
  )
})

test_that("max_dense_cells passes the forwarding gate on its two guarded targets", {
  expect_silent(hsquared:::hs_engine_control_forwarding(
    hs_control(engine = "julia",
               engine_control = list(target = "repeatability", max_dense_cells = 50)),
    "repeatability"
  ))
  expect_silent(hsquared:::hs_engine_control_forwarding(
    hs_control(engine = "julia", engine_control = list(max_dense_cells = 50)),
    "fit_animal_model"
  ))
  expect_error(
    hsquared:::hs_engine_control_forwarding(
      hs_control(engine = "julia",
                 engine_control = list(target = "random_regression", max_dense_cells = 50)),
      "random_regression"
    ),
    class = "hsquared_unsupported_syntax"
  )
})

# R1-1: the PR body claims the two single-step sites are "byte for byte"
# reproductions of the pre-fix call string, which Rose found false (the
# validated `iterations` is now forwarded unconditionally). Behaviour is
# unchanged today only because R's dispatch default (100L) equals the
# payload functions' own formal default (100L) -- a coincidence, not a
# guarantee. There is no existing seam in this file that captures or mocks
# the `julia_command` string, so this pins the two halves of that
# coincidence directly (both non-live): the dispatch-time default that
# `hs_engine_control_value()` resolves to when `iterations` is unsupplied,
# and the payload functions' own formal default. If either drifts from the
# other, this test catches it before a live Julia run would.
test_that("single_step_construct/metafounder_single_step default dispatch pins iterations = 100L (R1-1)", {
  default_control <- hs_control(engine = "julia")
  expect_identical(
    hsquared:::hs_engine_control_value(default_control, "iterations", 100L),
    100L
  )
  expect_identical(
    formals(hsquared:::hs_fit_julia_single_step_construct_payload)$iterations,
    100L
  )
  expect_identical(
    formals(hsquared:::hs_fit_julia_metafounder_single_step_payload)$iterations,
    100L
  )
})
