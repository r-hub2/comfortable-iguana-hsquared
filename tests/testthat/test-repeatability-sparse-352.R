# Scalable animal + permanent-environment model (HSquared.jl#352).
#
# The repeatability target fitted only through the DENSE estimator, whose engine
# guard is nobs^2 + nanimals^2 <= max_dense_cells. Any real repeated-measures
# pedigree blows that ceiling long before the model is the problem, so the
# standard repeated-measures animal model was unreachable in practice and V_A
# absorbed V_PE. `scale_method = "auto"` routes the SAME model through the
# engine's sparse K-effect AI-REML.
#
# The old `length(blocks) < 3L` floor in the multi-effect payload was an R-side
# accident: every engine entry point asserts only K >= 1.

# ---- pure R: no Julia required -------------------------------------------

test_that("engine_control accepts scale_method on the repeatability target", {
  expect_silent(
    hs_control(
      engine = "julia",
      engine_control = list(target = "repeatability", scale_method = "auto")
    )
  )
  expect_silent(
    hs_control(
      engine = "julia",
      engine_control = list(target = "repeatability", scale_method = "dense")
    )
  )
})

test_that("an unknown engine_control key is still refused on repeatability", {
  # Guards the hsquared#212 defect class: naming a lever the target does not
  # accept. Adding scale_method must not open the allowlist generally.
  expect_error(
    hsquared(
      y ~ animal(1 | id, pedigree = data.frame(
        id = c("a", "b", "c"), sire = c(NA, NA, "a"), dam = c(NA, NA, "b")
      )) + permanent(1 | id),
      data = data.frame(y = c(1, 2, 3, 1.5, 2.5), id = c("a", "b", "c", "a", "b")),
      control = hs_control(
        engine = "julia",
        engine_control = list(target = "repeatability", nonsense_key = 1)
      )
    ),
    "nonsense_key"
  )
})

test_that("the multi-effect payload floor is two blocks, not three", {
  # The floor is what made `animal + ONE i.i.d. effect` unreachable. A
  # one-block payload must still be refused; a two-block payload must get PAST
  # the block check (it then stops at the bridge-availability gate, which is
  # mocked to FALSE so this stays a pure-R test).
  fake <- structure(
    list(random_effects = list(list(name = "animal"))),
    class = "hs_bridge_payload"
  )
  expect_error(
    hsquared:::hs_fit_julia_n_effect_payload(fake),
    "needs at least two"
  )

  two <- structure(
    list(random_effects = list(list(name = "animal"), list(name = "permanent"))),
    class = "hs_bridge_payload"
  )
  testthat::local_mocked_bindings(
    hs_julia_bridge_available = function(...) FALSE,
    .package = "hsquared"
  )
  expect_error(
    hsquared:::hs_fit_julia_n_effect_payload(two),
    "requires Julia"
  )
})

test_that("permanent() is routed to the repeatability target, not multi_effect", {
  # The PE model keeps ONE surface. The router names the target and prints the
  # closest working call, which is better than multi_effect silently accepting
  # the same model under generic block labels.
  ped <- data.frame(
    id = c("a", "b", "c"), sire = c(NA, NA, "a"), dam = c(NA, NA, "b")
  )
  dat <- data.frame(y = c(1, 2, 3, 1.5, 2.5), id = c("a", "b", "c", "a", "b"))
  expect_error(
    hsquared(
      y ~ animal(1 | id, pedigree = ped) + permanent(1 | id),
      data = dat,
      control = hs_control(
        engine = "julia",
        engine_control = list(target = "multi_effect")
      )
    ),
    'needs `target = "repeatability"`'
  )
})

# ---- live Julia ----------------------------------------------------------

# Repeated-records dataset with gene-dropped breeding values. Drawing `a`
# i.i.d. instead of dropping it down the pedigree produces data with NO
# A-structured variance, which the fit then correctly assigns to the PE block --
# a degenerate fixture that looks like an engine bug. Hence hs_sim_genedrop_bv().
hs_rep_fixture <- function(n_per_gen, n_gen, reps, seed = 11,
                           sigma_a2 = 1.0, sigma_pe2 = 0.6, sigma_e2 = 0.8) {
  ped <- hs_sim_pedigree(n_founder = 40, n_per_gen = n_per_gen, n_gen = n_gen,
                         seed = seed)
  a <- hs_sim_genedrop_bv(ped, sigma_a2 = sigma_a2, seed = seed)
  set.seed(seed + 1)
  pe <- stats::setNames(
    stats::rnorm(nrow(ped), 0, sqrt(sigma_pe2)), ped$id
  )
  dat <- do.call(rbind, lapply(seq_len(reps), function(r) {
    data.frame(
      id = ped$id,
      y = 2 + a[ped$id] + pe[ped$id] +
        stats::rnorm(nrow(ped), 0, sqrt(sigma_e2)),
      stringsAsFactors = FALSE
    )
  }))
  list(ped = ped, dat = dat, n_animal = nrow(ped), nobs = nrow(dat))
}

hs_rep_fit <- function(fx, scale_method) {
  hsquared(
    y ~ animal(1 | id, pedigree = fx$ped) + permanent(1 | id),
    data = fx$dat,
    family = stats::gaussian(),
    control = hs_control(
      engine = "julia",
      engine_control = list(
        target = "repeatability", scale_method = scale_method
      )
    )
  )
}

test_that("dense and sparse routes agree below the dense ceiling", {
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required."
  )
  fx <- hs_rep_fixture(n_per_gen = 60, n_gen = 2, reps = 3)
  # Precondition: this fixture is UNDER the ceiling, so both routes can run.
  expect_lt(fx$nobs^2 + fx$n_animal^2, 1e6)

  fd <- hs_rep_fit(fx, "dense")
  fa <- hs_rep_fit(fx, "auto")

  vd <- variance_components(fd)
  va <- variance_components(fa)
  expect_equal(vd$component, c("animal", "permanent", "residual"))
  expect_equal(va$component, c("animal", "permanent", "residual"))
  # Two independent estimators on identical data: agreement here is the whole
  # justification for offering the sparse route as the same model.
  expect_equal(va$estimate, vd$estimate, tolerance = 1e-3)
  expect_equal(
    heritability(fa)$estimate, heritability(fd)$estimate,
    tolerance = 1e-3
  )
  expect_equal(
    repeatability(fa)$estimate, repeatability(fd)$estimate,
    tolerance = 1e-3
  )
})

test_that("the two routes' loglik differ by exactly the REML constant", {
  # NOT cosmetic: logLik() is not comparable across scale_method, so an AIC or
  # LRT spanning the two routes is wrong by this constant. Pinned so the offset
  # cannot drift silently, and so reconciling it engine-side fails loudly here
  # rather than passing unnoticed.
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required."
  )
  fx <- hs_rep_fixture(n_per_gen = 60, n_gen = 2, reps = 3)
  fd <- hs_rep_fit(fx, "dense")
  fa <- hs_rep_fit(fx, "auto")

  nfixed <- length(fd$result$fixed_effects)
  expected_gap <- (fx$nobs - nfixed) / 2 * log(2 * pi)
  expect_equal(
    fd$result$loglik - fa$result$loglik, expected_gap,
    tolerance = 1e-4
  )
})

test_that("the sparse route fits a model the dense route refuses", {
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required."
  )
  fx <- hs_rep_fixture(n_per_gen = 120, n_gen = 3, reps = 4)
  # Precondition: this fixture is OVER the ceiling.
  expect_gt(fx$nobs^2 + fx$n_animal^2, 1e6)

  # The refusal must name the lever that fixes it, not just "use a sparse route".
  expect_error(hs_rep_fit(fx, "dense"), "scale_method")

  fit <- hs_rep_fit(fx, "auto")
  expect_s3_class(fit, "hsquared_fit")
  expect_true(fit$result$converged)
  vc <- variance_components(fit)
  expect_equal(vc$component, c("animal", "permanent", "residual"))
  expect_true(all(vc$estimate > 0))
  # The point of the slice: V_A and V_PE are SEPARATED, not summed into V_A.
  # A single animal effect on this data would report roughly their sum.
  expect_gt(vc$estimate[vc$component == "permanent"], 0.1)
  expect_gt(heritability(fit)$estimate, 0)
  expect_lt(heritability(fit)$estimate, 1)
  # Per-individual effects come back id-labelled on both blocks.
  expect_equal(nrow(permanent_effects(fit)), fx$n_animal)
  expect_equal(nrow(breeding_values(fit)), fx$n_animal)
  expect_setequal(permanent_effects(fit)$id, fx$ped$id)
})

test_that("provenance names the estimator that actually ran", {
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required."
  )
  fx <- hs_rep_fixture(n_per_gen = 60, n_gen = 2, reps = 3)
  expect_equal(
    hs_rep_fit(fx, "dense")$result$diagnostics$variance_components,
    "estimated_repeatability_reml"
  )
  fa <- hs_rep_fit(fx, "auto")
  expect_equal(
    fa$result$diagnostics$variance_components,
    "estimated_repeatability_sparse_multi_effect_aireml"
  )
  expect_equal(fa$result$diagnostics$scale_method, "auto")
  # The sparse route DOES carry a repeatability interval now, formed from the
  # already-fitted components via multi_effect_sum_ratio_interval() rather than
  # by refitting densely.
  expect_s3_class(fa$result$repeatability_interval, "data.frame")
})

test_that("the sparse route's repeatability interval matches the dense one", {
  # Two different interval code paths -- the dense one REFITS from the raw
  # matrices, the sparse one differentiates the sparse loglik at the fitted
  # components. Below the ceiling both can run and must agree.
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required."
  )
  fx <- hs_rep_fixture(n_per_gen = 60, n_gen = 2, reps = 3)
  rid <- repeatability_interval(hs_rep_fit(fx, "dense"))
  ria <- repeatability_interval(hs_rep_fit(fx, "auto"))
  expect_equal(ria$estimate, rid$estimate, tolerance = 1e-3)
  expect_equal(ria$se, rid$se, tolerance = 1e-3)
  expect_equal(ria$lower, rid$lower, tolerance = 1e-3)
  expect_equal(ria$upper, rid$upper, tolerance = 1e-3)
  expect_true(ria$lower < ria$estimate && ria$estimate < ria$upper)
  expect_true(ria$lower > 0 && ria$upper < 1)
})

test_that("the sparse interval survives the above-ceiling fit", {
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required."
  )
  fx <- hs_rep_fixture(n_per_gen = 120, n_gen = 3, reps = 4)
  ri <- repeatability_interval(hs_rep_fit(fx, "auto"))
  expect_true(is.finite(ri$estimate) && is.finite(ri$se))
  expect_true(ri$lower < ri$estimate && ri$estimate < ri$upper)
})

test_that("animal + ONE bare (1 | group) effect now fits via multi_effect", {
  # The payoff of lowering the block floor from three to two: K = 2 through the
  # generic multi-effect target. Under the old floor this stopped at
  # "needs at least three random-effect blocks" no matter what the engine could
  # do, which is the half of HSquared.jl#352 that is not the PE model.
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required."
  )
  ped <- hs_sim_pedigree(n_founder = 40, n_per_gen = 60, n_gen = 2, seed = 5)
  a <- hs_sim_genedrop_bv(ped, sigma_a2 = 1.0, seed = 5)
  set.seed(99)
  grp <- factor(sample(paste0("g", 1:8), nrow(ped), replace = TRUE))
  geff <- stats::setNames(stats::rnorm(8, 0, sqrt(0.5)), levels(grp))
  dat <- data.frame(
    id = ped$id,
    grp = grp,
    y = 2 + a[ped$id] + geff[as.character(grp)] +
      stats::rnorm(nrow(ped), 0, sqrt(0.8)),
    stringsAsFactors = FALSE
  )

  fit <- hsquared(
    y ~ animal(1 | id, pedigree = ped) + (1 | grp),
    data = dat,
    family = stats::gaussian(),
    control = hs_control(
      engine = "julia",
      engine_control = list(target = "multi_effect")
    )
  )
  expect_s3_class(fit, "hsquared_fit")
  vc <- variance_components(fit)
  expect_equal(nrow(vc), 3L) # animal + grp + residual
  expect_true(all(vc$estimate >= 0))
})

test_that("scale_method = auto says it ignores initial and iterations", {
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required."
  )
  fx <- hs_rep_fixture(n_per_gen = 60, n_gen = 2, reps = 3)
  expect_warning(
    hsquared(
      y ~ animal(1 | id, pedigree = fx$ped) + permanent(1 | id),
      data = fx$dat,
      control = hs_control(
        engine = "julia",
        engine_control = list(
          target = "repeatability",
          scale_method = "auto",
          initial = c(sigma_a2 = 2, sigma_pe2 = 2, sigma_e2 = 2)
        )
      )
    ),
    "ignores `initial`"
  )
})

test_that("the sparse route returns variance-component and h2 standard errors", {
  # Before HSquared.jl#352's engine half, `variance_component_standard_errors`
  # had exactly one method (`fit::AnimalModelFit`), so EVERY multi-effect fit
  # had no SEs at all. Fitting the PE model correctly but losing the SEs would
  # have traded one reported problem for the other.
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required."
  )
  fx <- hs_rep_fixture(n_per_gen = 120, n_gen = 3, reps = 4)
  fit <- hs_rep_fit(fx, "auto")

  se <- variance_component_standard_errors(fit)
  expect_s3_class(se, "data.frame")
  # Same components, in the same order as variance_components(), so the two
  # tables bind side by side.
  expect_equal(se$component, c("animal", "permanent", "residual"))
  expect_equal(se$component, variance_components(fit)$component)
  expect_true(all(is.finite(se$se)))
  expect_true(all(se$se > 0))

  # Public shape is data.frame(term, se) matching heritability() (hsquared#236).
  # Internal storage stays a bare numeric; the extractor normalizes.
  h2se <- heritability_standard_error(fit)
  expect_s3_class(h2se, "data.frame")
  expect_equal(names(h2se), c("term", "se"))
  expect_equal(h2se$term, heritability(fit)$term)
  expect_true(is.finite(h2se$se) && h2se$se > 0)
  # A sane SE keeps the +/- 1 SE band inside (0, 1).
  h2 <- heritability(fit)$estimate
  expect_gt(h2 - h2se$se, 0)
  expect_lt(h2 + h2se$se, 1)

  # The permanent block's ratio is a variance-explained proportion, NOT a
  # heritability, so it is carried under its own name rather than as an h2 row.
  pse <- fit$result$permanent_proportion_se
  expect_true(is.numeric(pse) && length(pse) == 1L && is.finite(pse) && pse > 0)
})

test_that("the dense route still has no standard errors, and says so", {
  # Unchanged behaviour on the covered route: an honest refusal naming the
  # missing field, NOT a silent NA.
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required."
  )
  fx <- hs_rep_fixture(n_per_gen = 60, n_gen = 2, reps = 3)
  expect_error(
    variance_component_standard_errors(hs_rep_fit(fx, "dense")),
    "does not contain"
  )
})
