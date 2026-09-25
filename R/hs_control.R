#' Create hsquared control options
#'
#' `hs_control()` records execution and storage controls for hsquared model
#' calls. The default `engine = "fit"` fits the validated v0.1 Gaussian
#' animal model through the `HSquared.jl` engine (average-information REML);
#' `engine = "validate"` parses and validates the contract without fitting,
#' returning the validated spec invisibly; and `engine = "julia"` exposes
#' advanced engine controls.
#'
#' @param engine Execution engine. `"fit"` (default) fits the v0.1 Gaussian
#'   animal model via the `HSquared.jl` engine; this requires a local Julia,
#'   the `JuliaCall` package, and an `HSquared.jl` checkout. `"validate"`
#'   validates the parser and bridge payload without fitting, then returns the
#'   validated spec invisibly (after a confirming message). `"julia"` exposes
#'   the advanced opt-in bridge with explicit `target` control.
#' @param backend Planned compute backend. One of `"auto"`, `"cpu"`,
#'   `"threads"`, `"cuda"`, `"amdgpu"`, `"metal"`, or `"oneapi"`.
#' @param accelerator Planned accelerator preference. One of `"auto"`,
#'   `"none"`, `"gpu"`, `"cuda"`, `"amdgpu"`, `"metal"`, or `"oneapi"`.
#' @param precision Planned numeric precision. One of `"float64"` or
#'   `"float32"`.
#' @param save Planned fitted-object storage mode. One of `"minimal"`,
#'   `"full"`, or `"tiny"`.
#' @param engine_control A named list for engine-specific controls. The current
#'   experimental Julia bridge recognizes `julia_project`, `initial`,
#'   `iterations`, `em_warmup`, `target`, `variance_components`, `marginal`,
#'   and `max_dense_cells`.
#'   `julia_project` is honoured by every `target`; supplying a key a given
#'   `target` does not honour errors (hsquared#212) rather than being
#'   silently ignored. Per-target honoured keys (besides `julia_project`):
#'   * `fit_animal_model` (default): `initial`, `max_dense_cells`.
#'   * `henderson_mme`, `metafounder`, `snp_blup`: `variance_components`.
#'   * `ai_reml`: `initial`, `iterations`, `em_warmup`.
#'   * `sparse_reml`, `two_effect`, `direct_maternal`,
#'     `genomic`, `single_step`, `single_step_construct`,
#'     `metafounder_single_step`, `relmat`, `precision`: `initial`,
#'     `iterations`.
#'   * `repeatability`: `initial`, `iterations`, `max_dense_cells`,
#'     `scale_method`.
#'   * `multi_effect`: `initial`, `iterations`, `scale_method` -- `initial`/
#'     `iterations` are honoured on the `scale_method = "dense"` (default)
#'     route only; the opt-in `scale_method = "auto"` route does not yet
#'     forward them (HSquared.jl#343, a known remaining gap).
#'   * `multivariate`: `initial`, `iterations`, `genetic_structure`, `rank`.
#'   * `random_regression`: `iterations` (no `initial`).
#'   * `nongaussian`: `marginal`, `iterations`, `initial` (a list with
#'     `sigma_a2`), `restart_check`.
#'   `max_dense_cells` bounds `nobs^2 + nanimals^2` on
#'   the engine's dense-validation fitters (hsquared#214, #217): the default Julia
#'   target `target = "fit_animal_model"` (via `HSquared.fit_animal_model()` /
#'   `fit_variance_components()`) and `target = "repeatability"`. Both require
#'   `engine = "julia"`. It has **no effect** under the default `engine = "fit"`
#'   path, which routes to the sparse-capable `HSquared.fit_ai_reml()` and enforces
#'   no dense-cell cap at all. It must be a
#'   single positive integer; the default, `1e6`, mirrors the engine's own
#'   `DEFAULT_MAX_DENSE_CELLS` unchanged. Raise it to fit a larger dense
#'   problem at the cost of memory and time, or switch to a sparse route
#'   (`target = "ai_reml"`/`"sparse_reml"`) instead of raising it indefinitely.
#'   Exceeding the cap now raises an `hsquared_error` naming the observed cell
#'   count and the effective cap, rather than a raw Julia trace.
#'   On `target = "repeatability"`, `max_dense_cells` reaches the point fit only.
#'   The opt-in repeatability-coefficient interval is computed by a separate engine
#'   entry point that exposes no such control and runs at the engine default
#'   (`1e6`), so **raising the cap above `1e6` returns the point fit with the
#'   interval silently absent** (`NULL`) rather than with an error. Tracked on the
#'   engine side; not fixed here. This applies to `scale_method = "dense"` only:
#'   the sparse route forms its interval from the fitted components and never
#'   reaches that entry point, so it is unaffected by `max_dense_cells`.
#'
#'   `target` selects which Julia estimator the `engine = "julia"` bridge runs;
#'   it has no effect under the default `engine = "fit"` path. The supported
#'   targets are `"fit_animal_model"`, `"ai_reml"`, `"sparse_reml"`,
#'   `"henderson_mme"`, `"repeatability"`, `"two_effect"`, `"multi_effect"`,
#'   `"direct_maternal"`, `"random_regression"`, `"genomic"`,
#'   `"single_step"`, `"single_step_construct"`, `"metafounder"`,
#'   `"metafounder_single_step"`, `"snp_blup"`, `"relmat"`, `"precision"`,
#'   `"multivariate"`, and `"nongaussian"`, described below. Covered opt-in
#'   routes (validation scale; not the default path) include `"two_effect"`
#'   (`common_env()`), `"direct_maternal"`, `"multi_effect"`, and
#'   `"random_regression"` at `k = 2`. `marginal` applies only to
#'   `target = "nongaussian"`.
#'   With `engine = "julia"` and no `target`, the bridge defaults to
#'   `target = "fit_animal_model"`: it surfaces the Julia-owned
#'   `HSquared.fit_animal_model()` dense NelderMead optimizer, honouring the
#'   `REML` flag. This is **not** the same estimator as the default
#'   `engine = "fit"` path, which runs the validated average-information REML
#'   estimator (`HSquared.fit_ai_reml()`, the same one reached by
#'   `target = "ai_reml"`). The two paths target the same Gaussian animal model
#'   but differ at the optimizer level, so their estimates can differ slightly.
#'   Reach the validated estimator with `engine = "fit"` (the ordinary path) or
#'   with `engine = "julia"` and `target = "ai_reml"`.
#'   `target = "henderson_mme"` is a supplied-variance validation path and
#'   requires `variance_components` with named `sigma_a2` and `sigma_e2` values.
#'   `target = "metafounder"` is an experimental supplied-variance validation
#'   path for `metafounder(1 | id, pedigree = ped, group = mf_group, Gamma =
#'   Gamma)`. It builds the Julia-owned animal-only `A^Gamma` relationship;
#'   `Gamma` and the variance components are supplied, not estimated.
#'   `target = "sparse_reml"` is an experimental, opt-in validation path that
#'   surfaces the Julia-owned `HSquared.fit_sparse_reml()` REML-only sparse
#'   optimizer; it accepts `initial` (named `sigma_a2`/`sigma_e2`) and
#'   `iterations`. It is not the default, not production fitting, and not a
#'   variance-component estimation claim for the public R interface.
#'   `target = "ai_reml"` exposes the same average-information REML estimator
#'   (`HSquared.fit_ai_reml()`) that the default `engine = "fit"` path uses,
#'   with explicit `initial` and `iterations` control. This is the validated
#'   v0.1 estimator for the univariate Gaussian animal model; the `engine = "fit"`
#'   default is the ordinary way to reach it. `em_warmup` (a non-negative integer,
#'   default `0`) optionally runs that many EM-REML warm-start steps before the
#'   average-information iterations; it applies to this estimator on both the
#'   default `engine = "fit"` path and `target = "ai_reml"`. `em_warmup = 0` (the
#'   default) is byte-identical to the previous behaviour. A small positive value
#'   (e.g. `3`) can rescue convergence from poor starting variances; on an
#'   identified fit it converges to the same optimum.
#'   `target = "repeatability"` is an experimental, opt-in path for the
#'   repeatability (permanent-environment) model. It requires
#'   `animal(1 | id, pedigree = ped) + permanent(1 | id)` in the formula and
#'   surfaces the Julia-owned `HSquared.fit_repeatability_reml()` REML-only
#'   optimizer (three-component `initial` with `sigma_a2`/`sigma_pe2`/`sigma_e2`).
#'   It is REML only, not the default, and the additive and permanent-environment
#'   variances are identifiable only with repeated records per individual.
#'   `scale_method` selects the estimator. `"dense"` (the default) is the
#'   covered validation-scale `fit_repeatability_reml()` above, which the engine
#'   refuses once `nobs^2 + nanimals^2` exceeds `max_dense_cells` -- a ceiling
#'   any field-scale repeated-measures pedigree passes easily.
#'   `scale_method = "auto"` fits the SAME model through the engine's sparse
#'   K-effect AI-REML (`fit_multi_effect(method = :auto)`), expressing it as the
#'   animal block plus an identity-relationship permanent-environment block, and
#'   so has no dense ceiling. Below the ceiling the two agree to REML tolerance.
#'   Two differences are deliberate on the sparse route: `initial` and
#'   `iterations` are NOT honoured (the engine picks its own; a supplied value
#'   warns), and `loglik` carries the REML normalising constant
#'   `-(n - p)/2 * log(2 * pi)` that the dense route omits -- so
#'   log-likelihoods, AIC, or likelihood-ratio tests must NOT be compared
#'   across `scale_method` on the raw `loglik` (HSquared.jl #365). Fit
#'   `diagnostics` expose `loglik_convention`, `loglik_full_constant_offset`
#'   (add to `loglik` for the full-constant scale), and
#'   `loglik_comparable_across_routes` when the linked engine is new enough;
#'   the engine helper `comparable_loglik(fit)` does the same conversion.
#'   `repeatability_interval()` IS available on both
#'   routes and agrees between them; the sparse one is formed from the fitted
#'   components (`multi_effect_sum_ratio_interval()`) rather than by refitting
#'   densely, and both are asymptotic logit-scale delta intervals that are NOT
#'   coverage-calibrated.
#'   `target = "two_effect"` is an experimental, opt-in path for two-effect
#'   models. It requires `animal(1 | id, pedigree = ped)` plus a second random
#'   effect -- `common_env(1 | group)` (an IID common-environment effect) or
#'   `maternal_genetic(1 | dam)` (a maternal genetic effect carrying the pedigree
#'   relationship) -- and surfaces the Julia-owned
#'   `HSquared.fit_two_effect_reml()` REML-only optimizer (three-component
#'   `initial` with `sigma_a2`/`sigma_c2`/`sigma_e2`). It is REML only and not
#'   the default. The `common_env()` leg is covered at validation scale; the
#'   `maternal_genetic()` leg on this same target stays experimental. The
#'   covered correlated sibling is `target = "direct_maternal"`.
#'   `target = "multi_effect"` is an opt-in path for
#'   `animal(1 | id, pedigree = ped)` plus one or more bare `(1 | group)` IID
#'   intercepts. It is covered at validation scale (the independent
#'   generalization of two-effect; not the default path). Random slopes and
#'   correlated `(x || group)` terms remain rejected. The animal-block ratio is
#'   narrow-sense h2; other blocks are variance-explained proportions, not
#'   heritabilities. `initial` is a plain numeric vector of length K + 1 (one
#'   value per block, in formula order, plus the residual); it and
#'   `iterations` are honoured on the `scale_method = "dense"` route (see
#'   `engine_control` above for the `"auto"` gap).
#'   `target = "direct_maternal"` is an opt-in path for
#'   `animal(1 | id, pedigree = ped) + maternal_genetic(1 | dam)`. It estimates
#'   the correlated 2x2 direct-maternal genetic covariance and is covered at
#'   validation scale (not the default path). `heritability()` returns the
#'   labelled Willham triple (direct h2_d, maternal m2, total h2_T, r_am), never
#'   a bare scalar. `initial` is a list with `G_dm` (a 2x2 matrix) and
#'   `sigma_e2` (a positive scalar); it and `iterations` are honoured.
#'   `target = "random_regression"` is an opt-in reaction-norm path for
#'   `animal(rr(covariate, order = k) | id, pedigree = ped)`. It is covered at
#'   `k = 2` (linear reaction norm). `rr_heritability()` returns h2(t) as a
#'   curve; `heritability()` errors on this result and names
#'   `rr_heritability()`. `k >= 3` stays experimental.
#'   `target = "genomic"` is the explicit experimental narrow Gaussian REML
#'   genomic route. It accepts
#'   `genomic(1 | id, markers = M)`, using sample allele frequencies,
#'   unweighted VanRaden method 1, and ridge `0.01`, or
#'   `genomic(1 | id, Ginv = Ginv)`, whose construction method, allele-frequency
#'   source, ridge, and scale denominator remain unknown. `heritability()`
#'   returns the genomic variance-component ratio on the declared relationship
#'   scale, `sigma_g2 / (sigma_g2 + sigma_e2)`. Because relationship scaling
#'   enters through the kernel, this is not generally the fraction of average
#'   marginal phenotypic variance and is not pedigree-, founder-base-,
#'   population-, or universal narrow-sense heritability. Genomic
#'   `heritability_interval()` and
#'   `heritability_standard_error()` are unavailable. The R capability remains
#'   covered at validation scale, experimental, and opt-in; the explicit route
#'   does not move `public_covered_count` from 7. A fail-closed boundary candidate matched
#'   an independent oracle on all 240 sealed holdouts, correcting 30
#'   classifications with no losses, but one cell had a 5.99x p95 runtime ratio
#'   against the frozen 3x cap. The seeds are spent, the nine-cell campaign did
#'   not run, and the route must remain opt-in pending performance work and a new
#'   untouched holdout.
#'   `target = "single_step"` remains an experimental, opt-in path for
#'   `single_step(1 | id, Hinv = Hinv)` (a supplied single-step relationship
#'   inverse). Both targets surface `HSquared.fit_ai_reml()` and are REML-only;
#'   single-step is not a default route.
#'   `target = "single_step_construct"` fits
#'   `single_step(1 | id, pedigree = ped, markers = M)` after the engine builds
#'   `H^-1` from the pedigree and genotyped-subset markers. `target =
#'   "metafounder_single_step"` fits
#'   `single_step(1 | id, pedigree = ped, markers = M, group = mf_group, Gamma =
#'   Gamma)` through the Julia-owned supplied-`Gamma` `H^Gamma` path. Both are
#'   experimental, opt-in, dense/validation-scale, REML-only, and not
#'   comparator-validated; `Gamma` is supplied, not estimated. Both honour
#'   `initial` (named `sigma_a2`/`sigma_e2`) and `iterations`.
#'   `target = "snp_blup"` is an experimental, opt-in path for the SNP-BLUP /
#'   RR-BLUP marker-effect model. It requires `genomic(1 | id, markers = M)` (a
#'   raw marker matrix) and estimates per-marker effects (`marker_effects()`) and
#'   per-individual genomic breeding values. `variance_components` is **optional**:
#'   supply named `sigma_g2` (genomic) and `sigma_e2` (residual) for a
#'   supplied-variance solve (`HSquared.fit_snp_blup()`), or omit them to have
#'   `hsquared()` **estimate** `sigma_g2`/`sigma_e2` by REML from the markers
#'   (`HSquared.fit_snp_blup_reml()`). Not the default.
#'   `target = "relmat"` is an experimental, opt-in path for
#'   `relmat(1 | id, K = K)` (a supplied dense relationship matrix; the parser
#'   marshals the inverse). `target = "precision"` is the same experimental path
#'   for `precision(1 | id, Q = Q)` (a supplied precision/inverse). Neither is
#'   covered or the default; the supplied matrix is provenance, not an estimate.
#'   The `animal(1 | id, pedigree = ped)` route rejects selfing (rows with the
#'   same known sire and dam) in v0.1, with no argument that reaches the
#'   engine's `allow_selfing` flag; `relmat(1 | id, K = A)` with a hand-built
#'   or `AGHmatrix`-built relationship matrix `A` is the current workaround
#'   for a selfing or hermaphroditic pedigree.
#'   `target = "multivariate"` names the experimental multivariate Gaussian
#'   animal model. Naming it is optional: a `cbind()` Gaussian response with an
#'   `animal()` term auto-routes to this target on the default path, and under
#'   `engine = "julia"` with no `target`. The t = 2 unstructured capability is
#'   covered at validation scale, while remaining experimental. It requires a
#'   `cbind(trait1, trait2, ...)` response with `animal(1 | id, pedigree = ped)`,
#'   surfaces the Julia-owned `HSquared.fit_multivariate_reml()` REML-only
#'   optimizer, and returns G/R covariance matrices, genetic and residual
#'   correlations, per-trait heritability, and cross-trait breeding values.
#'   Three or more traits and structured covariance are experimental. The
#'   reserved `genetic_structure` control currently accepts `"unstructured"` and
#'   `"diagonal"` on the R bridge. `"diagonal"` is the rotation-free structured
#'   subset: off-diagonal genetic covariances are fixed at zero. `"lowrank"`
#'   remains planned. `"factor_analytic"` is planned on the R surface (Julia
#'   `V4-FA` is engine-covered at HSquared.jl `60895208` / #300) and is not
#'   activated on the R bridge; not covered. The future `rank` control is also
#'   reserved and currently errors instead of being ignored.
#'
#'   `target = "nongaussian"` is an experimental, opt-in conditional GLMM for
#'   `poisson(log)` or `binomial(logit)` (binary 0/1 or
#'   `cbind(successes, failures)` counts) with one intercept and
#'   `animal(1 | id, pedigree = ped)`. The `marginal` control selects a Laplace
#'   marginal likelihood approximation (`"laplace"`, default) or a variational
#'   ELBO (`"variational"`; aliases `"la"`/`"va"`). It reports the ratified
#'   conditional three-field contract: Poisson latent and count-scale observation
#'   h2; logit latent, liability, and numerically integrated observation-scale h2
#'   for Bernoulli or common-trial Binomial input. Varying trials return literal
#'   `NaN` with `"varying_trials_no_scalar_estimand"`, never a trial-count-averaged
#'   scalar. The ELBO is a lower bound on the marginal
#'   log-likelihood, so variational and Laplace `logLik`/`AIC` are **not**
#'   comparable. This path remains experimental and not coverage-calibrated.
#'   `initial` (hsquared#225) is a list with `sigma_a2`. The engine fits the
#'   single variance component with a **bracketed** Brent search over
#'   `log(sigma_a2)` on `log(initial$sigma_a2) +/- 6` -- there is no start
#'   value, and `initial` sets the **centre of the bracket**, so supplying it
#'   moves the whole search window. Unsupplied, the engine's own hard-coded
#'   `sigma_a2 = 1.0` centres it, giving `[exp(-6), exp(6)]`; a true `sigma_a2`
#'   outside that window cannot be reached without an `initial` on the scale of
#'   the data, which is what makes `initial` the retry lever for a boundary
#'   refusal. `restart_check` (hsquared#225, logical,
#'   default `FALSE`) opts into the engine's two-start restart
#'   (HSquared.jl#327): it refits once from a bumped second start and flags
#'   `boundary = TRUE` when the estimate moves with the start, catching a
#'   boundary a single fit can otherwise miss. The result carries `boundary`
#'   (hsquared#222) next to `converged`: `TRUE` means the fitted `sigma_a2` is
#'   a function of the search start, not the data. In practice a
#'   `boundary = TRUE` fit is refused before it reaches the R result at all --
#'   the Julia payload builder raises, translated into a classed
#'   `hsquared_julia_error` naming `initial`/`restart_check` as the retry
#'   levers -- so `boundary` on a returned fit is `FALSE`. Read it as
#'   `fit$result$boundary`, or via `fit_diagnostics()`'s `search_boundary` /
#'   `search_boundary_condition` rows (hsquared#230); those rows are
#'   unrelated to that function's `at_boundary` rows, which flag a variance
#'   component at or near zero.
#'
#' @return An object of class `"hs_control"`.
#' @export
hs_control <- function(
  engine = c("fit", "validate", "julia"),
  backend = c("auto", "cpu", "threads", "cuda", "amdgpu", "metal", "oneapi"),
  accelerator = c("auto", "none", "gpu", "cuda", "amdgpu", "metal", "oneapi"),
  precision = c("float64", "float32"),
  save = c("minimal", "full", "tiny"),
  engine_control = list()
) {
  engine <- match.arg(engine)
  backend <- match.arg(backend)
  accelerator <- match.arg(accelerator)
  precision <- match.arg(precision)
  save <- match.arg(save)

  if (!is.list(engine_control)) {
    stop("`engine_control` must be a list.", call. = FALSE)
  }
  if (length(engine_control) > 0L) {
    names_ok <- !is.null(names(engine_control)) &&
      all(nzchar(names(engine_control)))
    if (!names_ok) {
      stop("`engine_control` must be a named list.", call. = FALSE)
    }
  }
  if ("max_dense_cells" %in% names(engine_control)) {
    hs_validate_max_dense_cells(engine_control[["max_dense_cells"]])
  }

  structure(
    list(
      engine = engine,
      backend = backend,
      accelerator = accelerator,
      precision = precision,
      save = save,
      engine_control = engine_control
    ),
    class = "hs_control"
  )
}

hs_engine_control_value <- function(control, name, default) {
  if (!name %in% names(control$engine_control)) {
    return(default)
  }
  control$engine_control[[name]]
}
