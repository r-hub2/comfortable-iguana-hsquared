# hsquared (development version)

* **Authorship.** Added co-authors Yefeng Yang (ORCID 0000-0002-8610-4016) and
  Szymon Drobniak (ORCID 0000-0001-8101-6247) across `DESCRIPTION`, `inst/CITATION`,
  and `CITATION.cff`; Shinichi Nakagawa remains maintainer and corresponding author.
  No version bump; experimental **0.9.0** unchanged.

* **hsquared#236: `heritability_standard_error()` returns `data.frame(term, se)`.**
  Same public table shape as [heritability()], so
  `merge(heritability(fit), heritability_standard_error(fit), by = "term")`
  works. Sibling `variance_component_standard_errors()` stays
  `data.frame(component, se)`. Estimands and internal storage are unchanged
  (fits may still carry a bare numeric `heritability_se`); only the extractor
  return shape is aligned. **Breaking** for callers that treated the return as a
  length-1 numeric (use `$se`). No covered flip; `public_covered_count` stays
  **7**; version stays **0.9.0**.

* **Genomic ordinary-call honesty (T3).** Default-path `genomic()` now refuses
  with a pasteable `target = "genomic"` next call and names the held ordinary
  activation (design-44 G5; sealed `BOUNDARY_HOLDOUT_FAIL`). No silent default
  route, no covered flip, `public_covered_count` stays **7**, experimental
  **0.9.0** unchanged. Activation remains an owner science ticket.

* **hsquared#237 T1 close: `cbind` + `permanent` stays an honest fence, with a
  design note.** The engine still has no multi-trait permanent-environment
  fitter (`G0 ⊗ A + P0 ⊗ I`), so the named reject from parity-09 remains the
  live path. This slice banks `docs/design/57-mv-pe-cbind-permanent-237.md`,
  strengthens the multivariate repeated-records ABSORB warning to cite #237 and
  the two closest live routes (univariate `target = "repeatability"`, or
  `cbind` without `permanent()`), and adds capability / validation-debt rows.
  **No promotion** — `public_covered_count` stays **7**; version stays **0.9.0**.

* **R↔Julia `engine = "julia"` parity inventory at experimental 0.9.0 (parity-09).**
  Durable matrix `docs/design/55-r-julia-engine-julia-parity-09.md` classifies every
  Julia-fittable surface as REACHABLE / HONEST-ERROR / SILENT-GAP. SILENT-GAPs SG1–SG3
  closed with named unsupported-syntax tips (`cbind`+`permanent` #237; matfree
  engine-only; FA/lowrank stay planned). Skip-guarded smoke
  `tests/testthat/test-engine-julia-parity-smoke.R` names every REACHABLE target.
  **No promotion** — `public_covered_count` stays **7**; version stays **0.9.0**.
  Twin pointer: HSquared.jl `docs/design/12-bridge-compatibility.md`.

* **The sparse repeatability route asks the engine for its uncertainty once,
  not three times (HSquared.jl#370).** `target = "repeatability"` with
  `scale_method = "auto"` wanted three post-fit quantities — the
  variance-component standard errors, the ratio standard errors, and the
  repeatability interval — and called a separate engine function for each. All
  three are derived from the same asymptotic covariance, and that covariance is
  a finite-difference Hessian of the REML log-likelihood, the most expensive
  post-fit quantity in the engine, so it was being built three times. The route
  now calls `multi_effect_uncertainty()` once and reads all three from it.
  Measured on the great tit development data (11,856 records, 10,937 pedigree
  rows), end to end in R including the bridge: **6.15 s to 1.26 s**, with every
  reported number unchanged to the last printed digit — variance components,
  all three standard errors, the heritability standard error, and the
  repeatability interval and its endpoints. `hsquared` still works against an
  engine checkout that predates `multi_effect_uncertainty()`: it falls back to
  the three separate calls when the function is absent. (Szymon Drobniak.)

* **Fitting repeated records without a `permanent()` term now warns
  (HSquared.jl#352).** With more than one record per individual and only the
  additive effect in the model, the additive variance absorbs the
  permanent-environment variance and both it and the heritability are inflated
  — silently, until now. The warning names the counts that triggered it and the
  exact call that separates the two components, and is suppressible with
  `suppressWarnings()` when the single-effect model is intended.
  `formula_status()`'s `permanent(1 | id)` row now names the repeated-measures
  model and the `scale_method` lever. (Szymon Drobniak.)

* **Fixed the great tit development notebook's standard-error columns.**
  `dev-test/great_tit_animal_model.qmd` hardcoded `SE = NA_real_` for the
  `hsquared` rows in both the variance-component and heritability comparison
  tables and never called `variance_component_standard_errors()` or
  `heritability_standard_error()`. That, not an engine failure, is where the
  "hsquared produces no standard errors" report came from; the engine returns
  them on this route and always did. (Szymon Drobniak.)

* **The animal + permanent-environment model now fits at real pedigree scale
  (HSquared.jl#352).** `target = "repeatability"` gains
  `scale_method`. The default `"dense"` is unchanged; `scale_method = "auto"`
  fits the same `animal(1 | id, pedigree = ped) + permanent(1 | id)` model
  through the engine's sparse K-effect AI-REML, which has no
  `nobs^2 + nanimals^2` dense ceiling. Previously any field-scale
  repeated-measures dataset was refused by that ceiling, so users fell back to
  a single animal effect and `V_A` silently absorbed `V_PE`. Below the ceiling
  the two estimators agree to REML tolerance. On the sparse route `initial` and
  `iterations` are not honoured (a supplied value now warns), no
  `loglik` carries the REML
  normalising constant the dense route omits, so log-likelihoods must not be
  compared across `scale_method`. `repeatability_interval()` works on both
  routes and they agree. The dense route's refusal now names
  `scale_method = "auto"` instead of saying only "use a sparse route".
  (Szymon Drobniak.)

* **The sparse repeatability route returns standard errors
  (HSquared.jl#352).** `variance_component_standard_errors()` and
  `heritability_standard_error()` now resolve on a
  `scale_method = "auto"` repeatability fit, backed by new K-effect engine
  functions (`multi_effect_variance_component_standard_errors()`,
  `multi_effect_ratio_standard_errors()`). Previously the engine defined
  standard errors only for a single-animal-effect fit and for the multivariate
  fit, so every multi-effect model had none. The engine refuses at a flat or
  boundary optimum rather than returning `NaN`, and that refusal surfaces
  through the usual one-warning path rather than as a silent `NA`. The
  permanent block's ratio is reported as a variance-explained proportion, not
  as a heritability. The `"dense"` route still has no standard errors and still
  says so. (Szymon Drobniak.)

* **`animal(...)` plus a single `(1 | group)` effect now fits
  (HSquared.jl#352).** The multi-effect bridge required at least three
  random-effect blocks, an R-side accident — every engine entry point accepts
  `K >= 1`. The floor is now two, so `target = "multi_effect"` reaches the
  two-block model instead of erroring. (Szymon Drobniak.)

* Fixed stale documentation: `R/julia-bridge.R` comments described the Julia
  `prediction_error_variance()`/`reliability()` extractor defaults as
  `:dense`, which HSquared.jl#355 made false (the default is now `:auto`,
  selecting the inverse by storage; `:dense` remains the explicit validation
  oracle). No code or behaviour change. (HSquared.jl#364, Szymon Drobniak.)

* **Post-fit engine failures now surface instead of vanishing
  (HSquared.jl#351).** The bridge computes the experimental standard errors,
  intervals, and plot data inside Julia `try` blocks so a failure there never
  aborts the fit; previously the `catch` swallowed the error, so a fit whose
  `variance_component_standard_errors()` or `heritability_standard_error()`
  threw came back with those fields silently absent and no message, and
  nobody could tell "undefined here" from "the engine threw". Every such
  block (the SE trio, the multivariate covariance SEs, the repeatability /
  two-effect / n-effect ratio intervals, and all plot-data attachments,
  including an `OutOfMemoryError` there) now records the engine's message.
  After the fit, `hsquared()` raises ONE warning naming each engine function
  that failed with the first 200 characters of its message, and attaches the
  full record as `attr(fit, "bridge_errors")` (a named character vector).
  What is computed is unchanged; a clean fit carries no attribute and no
  warning. The underlying throw on the reported real-data fit is a separate
  engine question. No version bump; no capability-status change.
  Credit: Szymon Drobniak (finding).

* **`fit_diagnostics()` reports `search_boundary` for non-Gaussian fits,
  distinct from `at_boundary` (#230 / HSquared.jl#327).** Fits whose result
  carries a `boundary` field (the non-Gaussian bridge route) now gain a
  `search_boundary` row (mirroring `result$boundary`) and a
  `search_boundary_condition` row, which explains the bracket rail when
  `TRUE` and reads `"interior"` when `FALSE`.
  Gaussian animal-model fits (`engine = "fit"`) never carry a `boundary`
  field, so they gain neither row. `search_boundary` is unrelated to the
  existing `at_boundary` rows, which flag an estimated variance component at
  or near zero; `search_boundary` instead reports whether the optimizer's
  single-variance search stopped on the rail of its log-scale bracket. No
  version bump; no capability-status change.

* **The non-Gaussian bridge (`target = "nongaussian"`) forwards
  `engine_control$initial`/`$restart_check` and surfaces `boundary`
  (#222, #225 / HSquared.jl#342).** `initial` (a list with `sigma_a2`)
  recentres the engine's log-scale search bracket for the single variance
  component: `HSquared.fit_laplace_reml()` runs a bracketed Brent search over
  `log(sigma_a2)` on `log(initial$sigma_a2) +/- 6`, so supplying `initial`
  moves the whole search window, not just a starting point. Unsupplied, the
  engine's hard-coded `sigma_a2 = 1.0` centres it, giving `[exp(-6), exp(6)]`
  -- a true `sigma_a2` outside that window cannot be reached without an
  `initial` on the scale of the data; `restart_check` (logical, default
  `FALSE`) opts into the engine's two-start restart. The R result now
  carries `boundary` next to `converged` (read it as `fit$result$boundary`;
  a `fit_diagnostics()` row is tracked as #230, distinct from the existing
  `at_boundary` rows): `TRUE` means the fitted `sigma_a2` is a function of
  the search start, not the data. A fit that lands on its search boundary is refused by
  the Julia payload builder and surfaces as a classed
  `hsquared_julia_error`/`hsquared_error` naming `initial`/`restart_check` as
  the retry levers (previously the advice named these levers but the R bridge
  had no way to supply them); the R-side `hs_ng09_boundary()` guard also
  aborts with an `hsquared_boundary_refused`/`hsquared_error` condition on any
  future route that returns `boundary = TRUE` without refusing first. No
  version bump; no capability-status change.

* **Random-regression basis and evaluation points now error on out-of-range
  covariates instead of clamping (#213 / PR #219).** `hs_legendre_basis()`
  errors beyond the engine's own `1e-10` tolerance, matching the Julia
  `legendre_basis` boundary exactly (values within tolerance are still
  clamped to `[-1, 1]`, as before). `hs_rr_eval_points()` now validates a
  user-supplied `at` against the fitted covariate range and errors with the
  same `hsquared_out_of_range`/`hsquared_error`-classed condition, naming the
  fitted range and the offending value(s), before any random-regression
  trajectory extractor (`rr_genetic_variance()`, `rr_heritability()`,
  `rr_correlation()`, `rr_eigenfunctions()`, the `reaction_norm`/`rr_surface`
  plot types) can silently return the nearest fitted endpoint labelled with
  the requested out-of-range value.

* **Julia bridge errors are now translated into classed R conditions, and a
  new `max_dense_cells` lever raises the dense-validation-size cap on two
  guarded targets (#214, #217 / PR #220).** Every Julia-bridge fit call now
  raises an `hsquared_julia_error`/`hsquared_error`-classed condition
  carrying the original Julia message plus a one-line hint (e.g. pointing at
  `ridge`/`blend_weight` for an ill-conditioned `single_step()` block),
  instead of a raw, untranslated Julia error and stack trace.
  `hs_control(engine_control = list(max_dense_cells = ...))` raises or lowers
  the dense-validation-size cap (default `1e6`, unchanged) on two Julia
  targets: `target = "fit_animal_model"` (the default target under
  `engine = "julia"`) and `target = "repeatability"`, which is now guarded by
  the same cap for the first time. Both require `engine = "julia"`; the lever
  has **no effect** under the default `engine = "fit"` path, which routes to
  the sparse-capable `HSquared.fit_ai_reml()` and enforces no dense-cell cap
  at all. No other target was given the lever: `target = "multi_effect"` and
  `target = "direct_maternal"` do reach engine fitters that enforce the same
  guard, but they reach them through `HSquared.fit_payload_v2()`, which
  accepts no `max_dense_cells` argument, so those routes stay pinned at the
  engine default `1e6` with no R-side control — their bridge errors therefore
  carry a lever-free size hint rather than naming `max_dense_cells`. Repeatability's caveat: raising `max_dense_cells` above
  the default disables the repeatability *interval* silently (it returns
  `NULL`) rather than erroring, because `HSquared.repeatability_interval`
  runs the engine default internally and is not itself configurable by this
  lever — see `?hs_control`.

* **`engine_control$initial`/`iterations` are now forwarded on every target
  that documents them, and an unsupported key now errors (#212 / PR #223).**
  Five sites previously accepted `initial` and/or `iterations` with no
  effect: `target = "multi_effect"`, `target = "direct_maternal"`,
  `single_step_construct()`, `metafounder_single_step()`, and the
  `multi_effect_ratio_interval()` refit. All five now forward both controls
  into their Julia calls. A new per-target forwarding check errors with class
  `hsquared_unsupported_syntax` (naming the key and the target) when a caller
  supplies an `engine_control` key their chosen target does not honour; the
  full per-target table is in `?hs_control`. Known engine-side gap, not fixed
  here: `target = "multi_effect"` with `scale_method = "auto"` still does not
  forward `initial`/`iterations` (HSquared.jl#343).

* **Documentation family: engine contracts carried into the R surface, plus
  vignette fixtures (#208, #209, #210, #211, #215, #216, #218 / PR #221).**
  `?prediction_error_variance`, `?reliability`, and
  `?covariance_standard_errors` now state the engine's own PEV, reliability,
  accuracy, and Fisher-z interval-construction formulas, including that
  reliability/accuracy values are not clipped. `?genomic_markers` and the
  genomic-prediction vignette document the required 0/1/2 (or `[0, 2]`
  dosage) marker coding; a centered `-1/0/1` matrix must be recoded, and the
  out-of-range check now raises a classed `hsquared_error` instead of a bare
  `stop()`. `?heritability_interval` scopes the cited coverage figures to the
  design they were measured on (`q = 120`, interior h², level 0.95, DRAC job
  47925485). `sigma_P` is now written `sigma^2_P` throughout user-visible
  roxygen prose and interpretation strings, with no arithmetic change.
  `?animal`/`?hs_control` document that the `pedigree =` route rejects
  selfing rows in v0.1, with `relmat(1 | id, K = A)` named as the current
  workaround. `?factor_g_extractors` clarifies that `specific_variance()` is
  withheld only because `factor_analytic`/`lowrank` are not yet activated on
  the R-to-Julia bridge, not because `Psi` is unidentified. Six vignette
  articles gained a fixture chunk so their fit chunks resolve every object
  and stop only at the expected "requires the HSquared.jl Julia engine"
  message. `vignettes/articles/twin-boundary.Rmd` gained a "Twin contract
  rule" section pointing at the matching paragraph landed in HSquared.jl
  PR #338.

None of the four entries above changes `DESCRIPTION`, bumps the version, moves
a capability-status row, or changes `public_covered_count` (stays **7**). All
are R-repo-only bridge/validation/documentation fixes from the H² twin
independent test campaign.

# hsquared 0.9.0 (experimental release)

* **Ratified three-field non-Gaussian boundary.** The opt-in Poisson/log and
  Binomial/logit route retains its explicit scale boundaries: observation-scale
  h² is defined for Bernoulli and common-trial Binomial inputs, while varying
  trial counts return the named non-scalar sentinel rather than an averaged
  scalar. This does not increase `public_covered_count`, add a default route,
  or make a non-Gaussian REML claim.

* **Retained interval evidence, without promotion.** H0, H1, and H3 evidence
  is reported only at its declared route and claim level; no package-wide
  nominal-coverage or covered-status claim is created. The release remains
  experimental and not production; CRAN availability is tracked separately.

* **Default genomic route remains deferred.** The covered genomic-GREML
  estimand remains available only through explicit
  `engine = "julia", target = "genomic"`; the ordinary call is again the
  pedigree animal-model route. This restores the 0.9 scope boundary, changes
  no covered count, and does not alter the scientific fixture or evidence.

# hsquared 0.8.0

* **Interval-coverage claim levels for the univariate animal model (2000-rep
  C1 confirm; H0 Layer B).** The pre-registered C1 campaign
  (`docs/design/34-interval-recovery-pre-registration.md`) returned at the
  2000-rep confirm tier (HSquared.jl DRAC job 47925485). Mapped through
  doc-34 §4 at the worst interpretable interior cell: the **h² interval**
  is **directional-conservative** on all three legs — delta over-covers
  (0.969), profile is in-band (0.950), and bootstrap **mildly
  under-covers versus nominal 0.95** at Ĉ **0.924** (still ≥0.90). Do
  **not** describe h² as "never under-covers." The **σ²a profile** and
  **σ²a bootstrap** intervals are directional-conservative (profile
  0.963/0.947/0.956; bootstrap worst 0.918). The **σ²a delta/Wald**
  interval is **experimental-only** — it under-covers (0.897) and is
  shown as a point estimate ± SE only, not as a calibrated or
  conservative interval. `directional-conservative` is **not**
  coverage-calibrated-at-nominal: no leg is promoted to `point`. No
  capability is promoted — `public_covered_count` stays **7**. This is
  an honesty / claim-level change, not a covered flip, not a version
  bump, and not a 0.9 / 1.0 claim. **H0 Layer B ratified** (owner
  `ratify H0 Layer B`, 2026-09-07; Rose CLEAN WITH NITS; no covered flip).

* **Single-step opt-in routes are now distinguished in user-facing status
  guidance.** `single_step(..., Hinv = Hinv)` consumes a supplied relationship
  inverse, whereas `single_step(..., pedigree = ped, markers = M)` (or the
  matching `hs_data()` bundle) constructs H^-1. Both remain experimental and
  opt-in; `public_covered_count` stays **7**.

* **R factor-analytic boundary is now documented and contract-tested.** The
  planned R surface is documented separately from the Julia `V4-FA`
  engine-covered route; `cov = fa(...)` and
  `genetic_structure = "factor_analytic"` remain explicit planned-path
  rejections. This does not expose R loadings or an R-public FA fit, and
  `public_covered_count` stays **7**.

* **Twin-boundary article added for 0.9 prep.** New pkgdown page documents the
  R-public versus Julia-engine split, `public_covered_count = 7`, experimental
  0.8.0 (0.9 not released), and the FA/SS engine-only versus R-public fence.
  No covered flip.

* **Experimental number 0.7.0 → 0.8.0.** Both twins lockstep. Marks the
  0.8 engine pillar pair (Julia FA + single-step covered). R FA stays
  **planned** (#169). R `single_step()` stays **opt-in partial**.
  `public_covered_count` stays **7**. Experimental label **retained**.
  Not production, not CRAN, not 1.0.

* **R FA reject text catch-up (not a public FA).** Parser and
  `engine_control` still reject `cov = fa(K)` and
  `genetic_structure = "factor_analytic"`. Errors and `formula_status()`
  now name Julia `V4-FA` as engine-covered (`60895208` / #300) instead of
  the stale experimental #292 8/10 wording. R FA stays **planned**. Count
  stays **7**. No `cov = fa` parse. No version bump.

* **Single-step default-path errors now name both opt-in routes.** Rejections
  distinguish supplied-`Hinv` fitting (`target = "single_step"`) from
  engine-built `H^-1` construction (`target = "single_step_construct"`).
  Both remain experimental and opt-in; neither is a covered default R formula
  path. `public_covered_count` stays **7**.

# hsquared 0.7.0

* **Single-step H⁻¹ is no longer planned-only.** The Julia engine has
  construction AGREE vs AGHmatrix::Hmatrix and an n=240 recovery PASS
  (HSquared.jl #295). R keeps the existing opt-in
  `target = "single_step_construct"` / supplied-`Hinv` routes and does
  **not** promote `single_step()` to the ordinary default path. Status is
  **partial / experimental**, not covered. Darwin UNSIGNED.
  `public_covered_count` stays **7**. No General / CRAN / 1.0.

* **Factor-analytic G: engine-covered pointer, R row stays planned.** Julia
  `V4-FA` is engine-covered (HSquared.jl #300 / `60895208`; S4 evidence on
  #292). R `formula_status()` reserves `cov = fa()` and errors that the
  formula/bridge are not activated. The R capability row stays **planned**
  (post-G10 pointer). Not covered. `public_covered_count` stays **7**.
  Experimental stays **0.7.0**. No General / CRAN / 1.0.

* **Historical default genomic GREML activation (superseded for 0.9).**
  Owner YES 2026-09-03. Narrow `genomic(1 | id, markers = M)` /
  `genomic(1 | id, Ginv = Q)` now auto-routes on `engine = "fit"` (and on
  `engine = "julia"` with no explicit target) to the same covered
  validation-scale GREML estimand as `target = "genomic"`. Session warn-once
  keeps relationship-scale / experimental fences. Design-44 nine-cell recovery
  campaign is confirmatory / may still run — not claimed PASS. `public_covered_count`
  stayed **7**. Single-step / SNP-BLUP stayed opt-in. The 0.9 candidate
  restores explicit genomic routing. No General / CRAN / 1.0.

* **#7 - R genomic GREML `partial -> covered` (validation-scale, opt-in).**
  Owner ink: overnight approvals #5-#10 (auto-flip #7 when design-41 section 3 + Rose
  CLEAN). Darwin SIGN + Rose tip CLEAN after G1-G6 gap-clear (R #154 / Julia #287).
  `public_covered_count` moves **6 -> 7**. Version **0.7.0** experimental label
  **retained**. Scope: explicit `engine = "julia", target = "genomic"` on
  `genomic(1 | id, markers = M)` / `Ginv = Q`; estimand `genomic_variance_ratio`
  on \(K_\lambda = G + 0.01 I\) (design-51). **A later maintainer decision
  temporarily activated the default route for this same estimand (2026-09-03);
  the 0.9 candidate restores explicit routing.** **NOT** single-step / SNP-BLUP /
  APY / GPU / intervals. **NO-ANCHOR DISCLOSURE:** no clean Mrode genomic-\(h^2\)
  pin. Evidence: design-53 SUPERSEDE; Totoro exact-`G` PASS; Boole RATIFIED;
  canon VanRaden lock; promotion parity (DP-10 CI caveat). No General / CRAN / 1.0.

# hsquared 0.6.0

* **R multivariate Gaussian animal model `partial → covered` (validation-scale).**
  Owner ink: maintainer sign-off (2026-09-02). `public_covered_count` moves
  **5 → 6**. Scope: default-route `cbind(trait1, trait2) ~ animal(1 | id, pedigree = ped)`
  unstructured t=2 G0/R0 under design-41 §3 (items 1–6 + 8; #7 intervals deferred;
  #9 DoD in this flip). Experimental label **retained**. NOT production, NOT CRAN,
  NOT General force-merge. k≥3 and `genetic_structure = "diagonal"` stay experimental.
  Interval calibration stays uncalibrated. **NO-ANCHOR DISCLOSURE:** Mrode 5.1 is
  supplied-G0/R0 only; estimated G0/R0 have no published textbook anchor.
  Evidence pointers: capability-status multivariate row; A29 CLEAN; Wave B merge
  `55076481`; A27 Darwin SIGN; A26 local parity (DP-10 C).

# hsquared 0.5.0


This is an **experimental** numbered release (D-41). It is **not** a CRAN
submission. Experimental banners, lifecycle badges, `.onAttach`, pkgdown
callout, and DESCRIPTION wording stay on. `public_covered_count` stays
**5**. No covered flip. Multivariate stays `partial`. First public
*number* is 0.5.0; first public CRAN tarball is a later owner step after
Julia 0.5.0 is in General.

## New features

* **The `cbind()` multivariate animal model routes on the default call, and the documentation now says so.** A `cbind(trait1, trait2, ...)` Gaussian response with `animal(1 | id, pedigree = ped)` selects the multivariate REML fitter without `engine`/`target` arguments, and `engine = "julia"` with no `target` auto-selects `"multivariate"` for symmetry. That behaviour landed earlier (default routing, against the grammar predicate frozen in `docs/design/38`), but no public surface was updated with it: the multivariate article still said the model was "fitted **only** through `engine = "julia"` and `target = "multivariate"`", and the model-status article said `cbind(...)` responses "fit only through" that path. Both were false. The vignettes, README, `formula_status()`, `validation_status()` notes, the multivariate extractor error, capability-status, and the public-claims register now describe the default route, and the explicit spelling remains accepted so existing scripts keep working. **This is a reachability correction, not a promotion:** the multivariate capability stays `partial`/experimental, its covered-flip gate is unchanged, and `public_covered_count` remains 5. `formula_status()`'s fitting status for the `cbind()` row changes from `"fitted (opt-in multivariate)"` to `"fitted (default route, experimental multivariate)"`.

* **Genomic GREML hardening plus a banked negative activation gate.** The explicit `engine = "julia", target = "genomic"` route now has stricter marker/ID validation, engine-owned construction provenance, a deterministic cross-twin fixture, independent base-R construction, exact marker-versus-supplied-Q identity, and a fresh hash-pinned `blupf90+` point-estimate comparison. Marker input uses sample allele frequencies, unweighted VanRaden method 1, and `K = G + 0.01I`; supplied `Ginv` construction remains unknown. `heritability()` returns `genomic_variance_ratio = sigma_g2 / (sigma_g2 + sigma_e2)`, the genomic variance-component ratio on the declared relationship scale—not generally an average marginal phenotypic-variance fraction or pedigree/founder-base/population heritability. Ratio intervals and standard errors are unavailable. A fail-closed boundary candidate matched an independent base-R oracle on all 240 sealed holdouts, corrected 30 boundary classifications with no losses, and produced no invalid or unresolved fits. The conjunctive activation gate nevertheless failed because `n120_m600_r050` had a 5.99x candidate/default p95 runtime ratio against the frozen 3x cap. The 240 seeds are spent, the nine-cell recovery campaign did not run, and default routing was not added. The holdout used independent HWE/no-LD markers and provides no robustness evidence for LD, population structure, imputation, base-frequency misspecification, real panels, or production genotype data. The R row remains `partial`/experimental, the route is not a production genotype pipeline, and `public_covered_count` remains 5.

* **Large-scale multi-effect fitting: `engine_control` `scale_method = "auto"` (opt-in, experimental).** For the arbitrary-N independent `(1 | g)` multi-effect model (`target = "multi_effect"`), `hs_control(engine = "julia", engine_control = list(target = "multi_effect", scale_method = "auto"))` lets the engine route large problems the dense factorization cannot reach through its **matrix-free Monte-Carlo REML** fit (never forming/factorizing the mixed-model coefficient matrix). The default `scale_method = "dense"` is unchanged and byte-identical. At validation scale `"auto"` uses the sparse-exact estimator, which reduces exactly to the covered dense result (verified live: dense-vs-auto variance components agree to ~3e-5 on a 200-animal fixture); only when the problem exceeds the direct-factorization budget does it fall back to the experimental matrix-free fit, whose estimates carry a Monte-Carlo standard error. **The covered claim is unchanged** (validation-scale, exact path); the large-scale matrix-free path is opt-in **experimental** (no calibrated intervals / no external comparator through it yet).

* **Direct–maternal correlated (2×2 G) animal model — covered at validation scale (opt-in).** `hsquared()` now fits the correlated direct–maternal model through `hs_control(engine = "julia", engine_control = list(target = "direct_maternal"))` on `animal(1 | id, pedigree = ped) + maternal_genetic(1 | dam)`, surfacing the Julia-owned `HSquared.fit_direct_maternal_reml()`. This is the **first correlated random-effect model** (a 2×2 genetic covariance `G_dm` over the direct and maternal breeding values), distinct from the INDEPENDENT two-effect maternal leg (`target = "two_effect"`). It estimates four variance components (`σ²_ad` direct, `σ²_am` maternal, `σ_dm` the direct–maternal covariance, `σ²e` residual) and adds `genetic_correlation()` (the direct–maternal `r_am`), `direct_heritability()`, and `total_heritability()` (Willham total `h²_T`). **Willham fence:** `heritability()` on a `direct_maternal` fit returns the **labelled triple** — direct `h²_d = σ²_ad/σ_P`, maternal variance ratio `m² = σ²_am/σ_P`, Willham total `h²_T = (σ²_ad + 1.5·σ_dm + 0.5·σ²_am)/σ_P`, and `r_am` — **never a bare scalar** (`σ_P = σ²_ad + σ²_am + σ_dm + σ²e`); direct `h²` is NOT the total heritability, and a negative `r_am` is real and expected (`h²_T < h²_d` when `r_am < 0`). **Covered at VALIDATION scale** (mirrors the twin `V4-DIRECT-MATERNAL`, now covered): a PRE-DECLARED 48-seed bias/MCSE recovery gate PASSED (48/48 converged, all four `|bias| ≤ 2·MCSE`) and a `sommer` 4.4.5 `covm()` same-estimand REML comparator AGREE (`σ_dm` 1.1e-2, other entries ≤ 4.2e-3 relative). Experimental, opt-in, REML-only, dense/validation-scale (single relationship `A`, n ≤ ~1000; `|r_am| → 1` rides on `converged`); not the default, not ML, not production. A 2nd independent same-estimand comparator on a different lineage (`blupf90+` AIREMLF90 2×2-G; WOMBAT not installed), broader-DGP / larger-scale recovery, and a Mrode Ch.7 anchor remain owed (covered does not retire them).

* **Common-environment two-effect leg and arbitrary-N independent `(1 | g)` multi-effect model — covered at validation scale (opt-in).** The opt-in common-environment leg (`animal(1 | id) + common_env(1 | group)`, `target = "two_effect"`; additive genetic + IID common environment, A2 = I) and its arbitrary-N generalization (`animal(1 | id, pedigree = ped) + (1 | g1) + (1 | g2) + ...`, `target = "multi_effect"`; one pedigree animal component plus any number of INDEPENDENT i.i.d. effects) are now **covered at validation scale**, mirroring the twin `V3-TWOEFFECT-REML` / `V3-NEFFECT-REML` covered gates: a pre-declared 48-seed bias/MCSE recovery gate PASSED + a same-estimand REML comparator AGREE (`blupf90+` ~1e-5, `sommer` cross-check ~2e-5 for the two-effect leg; `sommer` for the multi-effect generalization, with exact live R–Julia parity). The animal-block ratio is narrow-sense `h²`; the other blocks are variance-explained proportions, not heritabilities; the per-component ratio intervals are asymptotic and not coverage-calibrated. The **maternal-genetic leg** of `target = "two_effect"` (A2 = pedigree) uses the same estimator with exact live parity but **stays experimental** (its own recovery gate + comparator on the maternal-A2 design are owed). INDEPENDENT effects only (not correlated / random-regression / non-Gaussian); opt-in, not the default, not ML, not production.

* **Genome-wide significance calibration for `gwas()` (`genome_wide = TRUE`).** `gwas(fit, markers, method = "single", genome_wide = TRUE, n_permutations =, seed =)` adds a genome-wide-calibrated `genome_wide_p` column, surfacing the Julia-owned `HSquared.genome_wide_marker_scan()`. For each analysis the permutation null is **rebuilt from the analysed phenotype** (`y` permuted conditional on `X`, re-scanned `n_permutations` times) and the genome-wide p is the exact add-one permutation p `(1 + #{null max >= observed})/(n_permutations + 1)` (significant when `genome_wide_p <= 0.05`). This is the **exact per-dataset rule** whose family-wise type-I control is validated at validation and production scale on HSquared.jl (the per-dataset add-one REBUILD gate, type-I 0.0504/0.0542 at alpha = 0.05); the anti-conservative `(1-alpha)` quantile threshold (HSquared.jl #202) is **not** used, and the mildly anti-conservative fixed-null-reuse simulation shortcut is **not** used. The result carries a `calibration` attribute (method `permutation_addone`; `empirical_type1 = NA` because the per-dataset rule's validity is by construction + externally validated, named in `validation_reference`; the chi-square genome-wide threshold is carried on the LOD scale). **SCOPE: fixed-effect / intercept-only only** -- `genome_wide = TRUE` requires `method = "single"`; the relatedness-corrected mixed-model / LOCO genome-wide null is a different, not-yet-validated calibration and is rejected. The nominal `p_value`/`bonferroni_p`/`bh_qvalue` columns remain NOT genome-wide calibrated (`print()` now scopes the caveat to those columns). Verified live element-wise against the engine. Experimental; the covered-claim promotion is twin-gated.

* **ggplot2 visualization layer (`autoplot()`).** New `autoplot()` methods render hsquared results as `ggplot2` figures, in the style of the `brms`/`bayesplot` ecosystem and consistent with the sister packages `drmTMB`/`gllvmTMB`: `autoplot(fit, "variance")` draws a variance-component **and** per-trait `h²` forest with experimental 95% intervals (`± 1.96·SE`, labelled asymptotic/REML); `autoplot(fit, "breeding_values")` draws a sorted EBV caterpillar with `± 1.96·√PEV` bands (faceted by trait for multivariate fits); `autoplot(fit, "g_matrix")` draws a **rotation-invariant** genetic-correlation heatmap of `G` for multivariate fits (raw factor loadings are never plotted — the ratified cross-lane convention; off-diagonal cells involving a low-`h²` trait are flagged as imprecise); `autoplot(fit, "g_geometry")` draws a rotation-invariant genetic-eigenstructure **scree** (eigenvalues + % variance explained; axis directions / loadings are never drawn); `autoplot(fit, "reaction_norm")` draws the genetic-variance and heritability trajectories of a random-regression fit; `autoplot(fit, "rr_eigenfunctions")` draws the rotation-invariant eigenfunctions `ψ_j(t)` of `K_g` (signs arbitrary, span-ambiguous under repeated eigenvalues); `autoplot(fit, "rr_surface")` draws the genetic covariance/correlation surface over the covariate grid; and `autoplot(gwas_result, "manhattan")` (default) / `autoplot(gwas_result, "qq")` draw a Manhattan plot and a QQ plot (the latter with a `lambda_GC` genomic-inflation diagnostic), both carrying the uncalibrated-significance banner. Also new: `hs_recovery_forest()` for known-truth recovery studies (bias `± 2·MCSE`, intervals covering zero = no detectable bias) and an exported `theme_hsquared()`. The figures are **uncertainty-first** and the helpers are modular (each takes a tidy data frame and returns a `ggplot`) so they can be factored into a shared visualization package later. Each figure also carries an `attr(p, "hsquared_meta")` honest-status attribute (type / interval status / rotation status / notes), implementing the §3 R-drawing contract of the cross-lane plotting architecture shared with `HSquared.jl` (`13-plotting-layer.md`). `ggplot2` is now a dependency; the base-R `plot()` method is unchanged.
* The Julia bridge now attaches the available engine `*_plot_data` payloads at fit time for standard animal-model, multivariate, and random-regression fits. `autoplot()` still recomputes when a payload is absent or when a random-regression user supplies a custom grid.
* `formula_status()` now lists the ratified planned missing-data grammar: `missing = miss_control(response = "include")` for future masked-response rows and `mi(x)` with `missing = miss_control(predictor = "model")` for a future modelled missing predictor. This is a grammar contract only; `mi()`, `miss_control()`, `impute_model()`, and `imputed()` are not exported and no missing-data fitting path is active.
* `gwas_table()` and `lod_scores()` now accept an already-computed `hs_gwas` object from `gwas(fit, markers)`, returning the current uncalibrated marker-scan table and marker-level LOD scores. Fit-level GWAS/QTL/eQTL tables, map joins, and calibrated thresholds remain planned (#23).
* Mirrored the Julia-owned genomic GBLUP/SNP-BLUP target fixture from HSquared.jl PR #140 (`008ea4d`) into R tests and validation/comparator ledgers. The Julia-free fixture check recomputes the supplied-frequency VanRaden `G`, `Ginv`, supplied-variance GBLUP MME solution, and SNP-BLUP marker-effect-to-GEBV route agreement. The fixture is available for future external comparator runs, but no AGHmatrix/rrBLUP/sommer/JWAS/BGLR/BLUPF90 comparator evidence is claimed and the genomic/SNP-BLUP rows remain partial.
* Mirrored the Julia-owned marker-scan result payload fixture from HSquared.jl PR #142 (`f9fbbb1`) into R tests. The fixture pins `gwas()` payload-normalization parity without live Julia, but it does not activate calibrated thresholds, map-annotated QTL/eQTL tables, external comparator evidence, or a covered claim.
* Reserved `metafounder_effects()` as an exported, error-only extractor so the future metafounder result surface has a stable name without implying that explicit metafounder solutions are currently returned. Existing `gamma_matrix()` and `metafounder_groups()` remain provenance-only extractors for supplied `Gamma` and group assignments.
* Added experimental live R bridges for supplied-`Gamma` metafounder relationships. `metafounder(1 | id, pedigree = ped, group = mf_group, Gamma = Gamma)` now parses, validates ID-keyed metafounder `group` labels and a finite symmetric positive-semidefinite supplied `Gamma`, and fits an animal-only supplied-variance `A^Gamma` model through `engine_control = list(target = "metafounder", variance_components = c(sigma_a2 = ..., sigma_e2 = ...))` by calling the Julia-owned `metafounder_animal_model()` path. `single_step(1 | id, pedigree = ped, markers = M, group = mf_group, Gamma = Gamma)` continues to fit through `target = "metafounder_single_step"` by calling `fit_metafounder_single_step_reml()`. `gamma_matrix(fit)` and `metafounder_groups(fit)` return the supplied `Gamma` matrix and group assignments as fitted-object provenance for these paths; they are not estimates. Skip-guarded live tests pin the animal-only `Gamma = 0` reduction to the ordinary Henderson MME supplied-variance path, the single-step `Gamma = 0` reduction to ordinary single-step construction, and nonzero-`Gamma` sensitivity probes with stable labels/dimensions. No `Gamma` estimation, returned metafounder-specific effects, BLUPF90 comparator evidence, production-scale claim, or covered support is claimed.
* Reconciled the PEV/reliability standard-field bridge status. Default, sparse, and explicit AI-REML Julia result-payload paths consume engine `prediction_error_variance` and `reliability` fields when present (current engines emit them via `:selinv`) and only use dense extractor calls as a backward-compatible fallback; the supplied-variance Henderson MME target now attaches dense validation-path PEV/reliability unconditionally. The bridge remains partial for multivariate per-trait PEV/reliability, production sparse strategy, and comparator validation (#21/#43).
* **Random-regression (reaction-norm) k=2 model — covered at validation scale (opt-in).** A new opt-in target surfaces the Julia-owned `HSquared.fit_random_regression_reml()`: `hsquared(weight ~ sex + animal(rr(age, order = 2) | id, pedigree = ped), data = long_records, family = gaussian(), REML = TRUE, control = hs_control(engine = "julia", engine_control = list(target = "random_regression")))`. `rr(covariate, order = k)` on the left of the `animal()` bar fits a `k`-coefficient normalized-Legendre polynomial of a within-individual covariate (default `order = 2` = intercept + slope); the covariate is standardized to `[-1, 1]` over its observed range and the bounds are recorded so extractors can re-standardize a user-supplied `at =` on the original scale. New extractors: `rr_covariance()` (the `k x k` coefficient genetic covariance `K_g`), `random_coefficients()` (per-animal predicted Legendre coefficients), and the reaction-norm trajectories `rr_genetic_variance(fit, at =)`, `rr_heritability(fit, at =)`, `rr_correlation(fit, at =)`, and `rr_eigenfunctions(fit, at =)` (the rotation-invariant eigen-decomposition of `K_g` as covariate functions — eigenvalues, proportion of genetic variance explained, sign-canonicalized eigen-coefficients, and eigenfunctions `psi_j(t)`; live-verified `== HSquared.rr_eigenfunctions()` to ~1e-15) — all computed in R from `K_g` and the recorded basis, matching the engine convention. The grammar (`rr(...)` inside `animal()`) was **ratified by the Julia twin** on `HSquared.jl#61` (normalized Legendre on standardized `t ∈ [-1, 1]`; Kirkpatrick/Meyer/Schaeffer). Covered at validation scale at k=2 (mirrors the twin random-regression REML validation row: a pre-declared 48-seed bias/MCSE recovery gate PASSED + a `sommer` 4.4.5 `leg()` same-estimand REML comparator AGREE; exact live R–Julia parity), opt-in, REML-only, single-effect, univariate, dense (k=2 only; k≥3 experimental); mirrors the twin `#54` engine. **Heterogeneous residual variance and a permanent-environment term are still planned** — with the current homogeneous residual and no permanent-environment effect, `rr_heritability()` can OVERSTATE `h^2(t)` for repeated-records designs (test-day, growth curves). Multivariate random regression and combining `rr()` with a second random effect are planned, not implemented (#54).
* **Experimental `gwas(fit, markers)` post-fit marker scan.** Runs a dense, supplied-variance, **relatedness-corrected** mixed-model (GLS) Wald marker scan on a fitted Gaussian animal model, reusing the fit's estimated variance components and pedigree relationship (surfacing the Julia-owned `HSquared.mixed_model_marker_scan()`). Returns a per-marker table (`effect`, `se`, `z`, `chisq`, `p_value`, `bonferroni_p`, `bh_qvalue`, `lod`). **The p-values are NOT genome-wide calibrated** — they are nominal Wald p-values plus Bonferroni/BH over the supplied markers only, with no R threshold activation, no permutation-backed cutoff, and no external comparator. HSquared.jl PR #134 banked a fixed-marker-panel calibration smoke harness, but it does not activate R significance thresholds or provide realistic-LD production calibration; the wrapper applies one whole-pedigree relationship correction. `gwas(fit, markers, method = "single")` additionally surfaces the relatedness-**un**corrected single-marker (OLS) scan (`HSquared.single_marker_scan()`) as a naive contrast — it carries a `scan_method` attribute and `print()`/`autoplot()` flag the absence of any relatedness correction. `gwas(fit, markers, method = "loco", marker_groups = chrom)` surfaces the leave-one-group-out scan (`HSquared.loco_mixed_model_marker_scan()`): a marker is corrected by a genomic (VanRaden) relationship built from all **other** marker groups, so its own signal does not leak into the background relationship. The LOCO relationship is genomic while the reused variance components are pedigree-estimated (a scale mismatch surfaced in `print()`/docs); `print()`/`autoplot()` flag the LOCO correction and `print()` restates the calibration caveat for every method. Verified live to match the engine element-wise (each method, including a per-group LOCO precision-selection check). Experimental, dense/validation-scale; the reserved tabular `gwas_table()`/`qtl_table()`/`eqtl_table()` extractors stay reserved for the planned map-annotated API (#23).
* **Experimental, opt-in single-step H⁻¹ *construction*.** `hsquared()` now accepts `single_step(1 | id, pedigree = ped, markers = M)` and fits it through `engine_control = list(target = "single_step_construct")`: the engine builds the pedigree inverse `A⁻¹` and dense `A` from the pedigree and the genomic relationship `G` from the genotyped-subset markers, assembles the single-step relationship inverse `H⁻¹` (Aguilar et al. 2010), and fits by REML — so you no longer have to precompute `Hinv` yourself (the supplied-`Hinv` form still works). Genotyped animals must be in the pedigree, but phenotyped animals need **not** be genotyped (the point of single-step), and GEBVs are returned for **all** pedigree animals. Construction knobs (`tau`, `omega`, `blend_weight`, `ridge`) are exposed but not comparator-validated. When the data is an `hs_data()` container that bundles a pedigree and genotypes, `single_step(1 | id)` resolves both from the bundle (the `animal(1 | id)` precedent), so neither `pedigree =` nor `markers =` is required; explicit arguments still override the bundle. Experimental, opt-in, REML-only, dense/validation-scale; mirrors the twin `V2-SSHINV` (partial). Verified live: marker-row-order invariance (the genotyped-rows alignment guard), id-labelled GEBVs covering ungenotyped animals, a differs-from-pedigree-model anchor, and the `hs_data()` shorthand fitting identically to the explicit call. Not the default; promotion past `partial` is twin-gated (`docs/design/25`).
* `formula_status()` now reports the `single_step(1 | id)` `hs_data()` bundle shorthand as its own parsed opt-in row, separate from explicit `single_step(1 | id, pedigree = ped, markers = M)` and supplied-`Hinv` single-step forms.
* **G-matrix geometry / evolvability extractors** (Hansen & Houle 2008) for opt-in multivariate fits: `eigen_G()` (the genetic eigenstructure — variance per genetic axis + the genetic principal components, the reserved name now implemented), `g_max()` (the leading genetic axis), `mean_evolvability()`, and the directional `evolvability()`, `variance_along_gradient()`, `respondability()`, `conditional_evolvability()`, and `autonomy()`. These are **rotation-invariant** functionals of the estimated `G` (not of factor loadings), so they are well defined for any multivariate fit and need no loading-rotation convention — the agreed cross-lane convention for what the structured-covariance bridge may surface (`HSquared.jl#42`/`#55`). Computed in R from `genetic_covariance(fit)` and verified to match the engine's `evolvability.jl` definitions by a live parity test; experimental, REML-only, no standard errors, carrying the multivariate fit's `partial` status (#55).
* **Experimental, opt-in non-Gaussian (GLMM) animal model.** `hsquared()` accepts `poisson(log)` or `binomial(logit)` on an intercept-only `animal(1 | id, pedigree = ped)` model through `hs_control(engine = "julia", engine_control = list(target = "nongaussian"))`. The versioned 0.9 conditional three-field result reports Poisson latent and count-scale observation h2; and logit latent, liability, and numerically integrated observation-scale h2 for Bernoulli or common-trial Binomial input. Varying-trial Binomial deliberately reports `h2_observation = NaN` with `h2_observation_undefined_reason = "varying_trials_no_scalar_estimand"`: no trial-count averaging or scalar observation-scale claim is made. `binomial(logit)` accepts binary 0/1 Bernoulli input or `cbind(successes, failures)` counts with scalar or varying positive trial totals. The objective is selected with `engine_control$marginal`: `"laplace"` is a Laplace marginal-likelihood approximation and `"variational"` (aliases `"la"`/`"va"`) is an ELBO; neither is called REML/AI-REML and their `logLik`/`AIC` values are not comparable. Experimental and opt-in only: no coverage calibration, external same-estimand comparator, capability promotion, or covered claim; not the default (#44). The "Fitting models" article now has a section, *Binary and count traits (experimental): which heritability to report*, that walks a varying-clutch fledging example and explains in plain terms that the liability-scale row is the one to report, why the observation-scale row is `NaN` when trial sizes vary, and how to get a data-scale value from `QGglmm` if one is needed.
* **Experimental:** `heritability_interval()` extracts a large-sample confidence interval for `h²` from the default Gaussian animal-model fit, when a local Julia engine returns one. It is a REML-only, asymptotic (logit delta-method) interval — mirroring the engine row `V1-HERIT-CI` (`partial`), so it is not coverage-calibrated and is unreliable at small `n`. It is reported as a point estimate plus bounds, not a validated capability (#11).
* **Experimental:** `variance_component_standard_errors()` and `heritability_standard_error()` return large-sample (delta-method) standard errors from the REML average-information matrix, when the default Gaussian fit's engine provides them. Same `V1-HERIT-CI` (`partial`) caveats: asymptotic, REML-only, omitted near a variance-component boundary (ill-conditioned AI matrix), not coverage-calibrated, not a validated capability.
* **Experimental:** `repeatability_interval()` returns a logit delta-method confidence interval for the repeatability coefficient `t = (Va + Vpe) / Vp` from the opt-in repeatability model, when the engine provides one. It mirrors the engine row `V3-REPEAT-REML` (`partial`): engine-internal self-consistency tested (recovery of `t` + interval bracketing on seeded fixtures) but with no external comparator, no `h²` interval, and no deep-pedigree validation; reported as a point estimate plus bounds only (#12).
* `summary()`/`print()` for `hsquared_fit` now display the experimental heritability confidence interval, variance-component and heritability standard errors, and repeatability interval when a fit carries them, clearly labelled experimental and asymptotic (#28).
* New "Benchmark: hsquared vs sommer and pedigreemm" article documents the v0.1 Gaussian animal-model fit agreeing with `sommer` and the published gryphon anchor within the signed-off band, with reproducing code and the `pedigreemm` one-sided log-likelihood floor (#31).
* Added a base-graphics `plot()` method for `hsquared_fit`: `type = "variance"` plots the variance components (with experimental `+/- 1.96 SE` whiskers when present) and `type = "residuals"` plots residuals against fitted values (#30).
* **Experimental:** `covariance_standard_errors()` returns delta-method standard errors for the multivariate genetic/residual covariance and correlation matrices and per-trait `h²`, for an opt-in **unstructured** multivariate fit when the engine provides them. Mirrors the engine multivariate REML validation row (`partial`): asymptotic, REML-only, unstructured-only, not coverage-calibrated. The strict per-seed recovery gate is still a non-pass (7/12 seeds), but the 12-seed bias/MCSE study (twin `HSquared.jl#78`/`#79`) shows **no detectable bias** (`|bias| ≤ 2·MCSE` for all six covariance parameters) — reported, not yet a validated capability (#26).
* New "A worked animal-model analysis (gryphon)" article walks one univariate animal model end to end — fit, heritability with experimental CI/SE, breeding values and accuracy, and the diagnostic `plot()` — on the gryphon teaching dataset (#29).
* Added a second published external-canon anchor: a CI-runnable, Julia-free test pinning the reference Henderson solver to the published Mrode (2014) Example 3.2 **sire-model** solutions (p.48), extending the canon to a second model class (#32).
* **Experimental:** the opt-in multivariate target now accepts `engine_control = list(genetic_structure = "diagonal")` (a diagonal genetic covariance — per-trait genetic variances with zero genetic covariances; no rotation ambiguity). `"lowrank"`/`"factor_analytic"` remain gated on a validated rotation convention. New `covariance_structure_lrt(constrained, full)` reports the experimental diagonal-vs-unstructured likelihood-ratio test from two multivariate fits on the same data — statistic `2*Δloglik`, with `df` the difference in the two fits' genetic-covariance parameter counts (the off-diagonal genetic covariances, i.e. `t(t-1)/2` for the diagonal-in-unstructured case), an interior null (`boundary = FALSE`), χ². It mirrors the engine multivariate REML validation row (`partial`): asymptotic, REML-only, not a validated test (#47).
* **Validation evidence (multivariate, t = 2).** Two reproducible studies plus one Bayesian agreement probe now back the experimental multivariate target (all `.Rbuildignore`d, not part of the build). `data-raw/multivariate-recovery-study.R` is a 100-replicate **cold-start** known-truth recovery study in which all six G0/R0 elements, the genetic correlation, and both per-trait `h²` fall within bias ± 2·MCSE (100/100 converged; EBV accuracy 0.79/0.74) — **no detectable bias**, agreeing with the twin's 12-seed study (`HSquared.jl#78`/`#79`); being cold-started from the identity, it is not a warm-start artifact. `data-raw/multivariate-comparator-study.R` is a **full-unstructured-residual** `sommer` external comparator that reproduces the engine's serialized `phase4_multitrait_parity` G0/R0/β/`h²`/EBV to ≤ 8e-5 and additionally recovers the off-diagonal residual covariance the in-suite diagonal-residual `sommer` check cannot. `data-raw/multivariate-mcmcglmm-agreement-study.R` is a **Bayesian agreement probe**: the serialized Julia target lies inside 95% HPD intervals for all 8 covariance elements, all 4 fixed effects, and both per-trait `h²` values, with posterior-mean EBV correlations > 0.9997. The `MCMCglmm` leg is not same-estimand REML parity, and the multivariate capability **stays `partial`** — this is evidence toward, not promotion of, the twin-gated covered gate (#10).

# hsquared 0.1.0

## New features

* **The default `hsquared()` call now fits the v0.1 univariate Gaussian animal model** `y ~ fixed + animal(1 | id, pedigree = ped)` by REML (average-information) through the `HSquared.jl` engine, returning heritability, variance components, breeding values, fixed effects, fitted values, residuals, and diagnostics. Fitting requires a local Julia, the `JuliaCall` package, and an `HSquared.jl` checkout; without them the default call errors with install guidance, and `hs_control(engine = "validate")` validates the contract without fitting. ML is not implemented. `REML = FALSE` is rejected on the default fit and on `engine = "validate"`; the covered claim is this REML estimator, not ML. The fit is validated by known-truth recovery, the published gryphon REML anchor (within the signed-off comparator band), and `sommer` agreement (#6, #7). These engine-recovery checks run locally through the R-to-Julia bridge (a local Julia and `HSquared.jl` checkout are required); public CI exercises the equivalent pure-R REML reference and skip-guards the live-engine tests.
* `animal()` is now exported as an inert formula marker, and `hsquared()` parses the narrow v0.1 formula contract `animal(1 | id, pedigree = ped)` (#4, #6).
* `EBV()` and `BLUP()` now alias `breeding_values()` for `hsquared_fit` objects, and `accuracy()` derives square-root reliability when reliability estimates are present (#5).
* `fitted()` and `residuals()` now work for `hsquared_fit` objects that contain fitted-value predictions and response values (#5).
* Added `fit_diagnostics()` for `hsquared_fit` objects so users and developers can inspect engine, method, target, convergence, optimizer status, iterations, log-likelihood, and dense-validation-path metadata without refitting or implying production support (#5, #6).
* `coef()` now aliases fixed-effect extraction for `hsquared_fit` objects, and `nobs()` reports the result observation count or response-payload fallback when available (#5).
* `hs_control()` has an `engine = "julia"` option that exposes advanced engine targets (supplied-variance Henderson MME, the opt-in sparse REML optimizer, explicit `ai_reml` control, repeatability/two-effect/genomic/SNP-BLUP targets, and the opt-in multivariate target) through a sibling `HSquared.jl` checkout via JuliaCall. The default `engine = "fit"` already fits the v0.1 model via `ai_reml` (non-Gaussian fitting was planned at the time of this entry; see the development-version entry above for the opt-in `target = "nongaussian"` path) (#6).
* `hs_control()` now recognizes an experimental, opt-in `engine_control = list(target = "sparse_reml", initial = ..., iterations = ...)` path that surfaces the Julia-owned `HSquared.fit_sparse_reml()` REML-only sparse optimizer through the local bridge. It is not the default `hsquared()` behaviour, not variance-component estimation in the public R interface, and not a production sparse fitting, AI-REML, or ASReml-parity claim (#6, #7).
* **Experimental, opt-in repeatability (permanent-environment) model.** `hsquared()` now parses `animal(1 | id, pedigree = ped) + permanent(1 | id)` and fits it through `engine_control = list(target = "repeatability")`, surfacing the Julia-owned `HSquared.fit_repeatability_reml()` REML-only optimizer. It returns three variance components (animal, permanent, residual), `repeatability()`, `heritability()`, breeding values, and `permanent_effects()`. The default `engine = "fit"` path stays single-effect and rejects `permanent()` with a pointer to the opt-in target. This is experimental and REML-only; the additive (σ²a) and permanent-environment (σ²pe) variances are identifiable only with repeated records per individual. It is not the default, not ML, and not yet a production or comparator-validated claim.
* **Experimental, opt-in common-environment (two-effect) model.** `hsquared()` now parses `animal(1 | id, pedigree = ped) + common_env(1 | group)` and fits it through `engine_control = list(target = "two_effect")`, surfacing the Julia-owned `HSquared.fit_two_effect_reml()` REML-only optimizer (additive genetic effect + an IID common-environment effect, e.g. litter or cage). It returns three variance components (animal, common_env, residual), `heritability()`, breeding values, and `common_env_effects()`. The default `engine = "fit"` path stays single-effect. This is experimental and REML-only; not the default, not production or comparator-validated. The correlated direct–maternal (2×2 G) model is a separate opt-in target (`target = "direct_maternal"`; now covered at validation scale — see the development-version entry).
* **Experimental, opt-in maternal-genetic two-effect model.** `hsquared()` now parses `animal(1 | id, pedigree = ped) + maternal_genetic(1 | dam)` and fits it through `engine_control = list(target = "two_effect")`: a direct additive genetic effect plus a maternal genetic effect expressed through the dam, both carrying the pedigree relationship (A₂ = pedigree A). It returns three variance components (animal, maternal_genetic, residual), `heritability()`, breeding values, and `maternal_effects()`; the dams must be animals in the `animal()` pedigree. Experimental and REML-only; this two-effect maternal leg treats the direct and maternal effects as INDEPENDENT. The correlated direct–maternal (2×2 G) model is a separate opt-in target (`target = "direct_maternal"`; now covered at validation scale — see the development-version entry).
* **Experimental, opt-in genomic GREML model.** `hsquared()` now parses `genomic(1 | id, Ginv = Ginv)` — a primary genomic effect with a user-supplied genomic relationship inverse — and fits it through `engine_control = list(target = "genomic")`, surfacing `HSquared.fit_ai_reml()` on a Ginv-based spec (genomic REML). It returns the genomic and residual variance components, genomic `heritability()`, and genomic breeding values (GEBVs). You can either supply a precomputed `Ginv`, or pass a raw marker matrix (`genomic(1 | id, markers = M)`) and let the engine build the genomic relationship and its inverse. The default `engine = "fit"` path is the pedigree animal model; `genomic()` is opt-in. Experimental and REML-only; comparator validation and weighted/standardized-marker variants remain planned.
* **Experimental, opt-in single-step model.** `hsquared()` now parses `single_step(1 | id, Hinv = Hinv)` — a primary effect with a user-supplied single-step relationship inverse — and fits it through `engine_control = list(target = "single_step")`, reusing the same `fit_ai_reml`-on-a-supplied-inverse path as `genomic()`. Experimental and REML-only; building `Hinv` from a pedigree and genomic relationship (single-step HBLUP construction) and comparator validation remain planned.
* **Experimental, opt-in SNP-BLUP / RR-BLUP marker-effect model.** `hsquared()` now fits `genomic(1 | id, markers = M)` through `engine_control = list(target = "snp_blup", variance_components = c(sigma_g2 = ..., sigma_e2 = ...))`, surfacing the Julia-owned `HSquared.fit_snp_blup()`. At the supplied genomic/residual variances it estimates per-marker effects — extractable with `marker_effects()` — together with per-individual genomic breeding values and fixed effects (the engine centres the markers, VanRaden method 1). `marker_variance_explained()` reports a descriptive fitted-marker contribution table for this path, computed from effect squared times centered marker variance and normalized across markers; it is not a marker-scan, p-value, or QTL claim. When `variance_components` are supplied this is a supplied-variance solve; **when they are omitted, `hsquared()` now estimates σ²g/σ²e by REML** from the markers (surfacing the Julia-owned `HSquared.fit_snp_blup_reml()`) — so `genomic(1 | id, markers = M)` with `target = "snp_blup"` and no variances fits end-to-end. The REML path reports `variance_components_source = "estimated_snp_blup_reml"` and is live-verified to match a direct engine `fit_snp_blup_reml()` element-wise. Opt-in, not the default; mirrors the twin `V2-SNPBLUP` gate (the GBLUP↔SNP-BLUP genomic-breeding-value equivalence). Weighted/Bayesian marker priors and comparator parity remain planned.
* **Experimental, opt-in multivariate Gaussian animal model.** `hsquared()` now parses `cbind(trait1, trait2, ...) ~ fixed + animal(1 | id, pedigree = ped)` and fits it through `engine_control = list(target = "multivariate")`, surfacing the Julia-owned `HSquared.fit_multivariate_reml()` REML-only dense estimator. It returns G/R covariance matrices, genetic and residual correlations, per-trait heritability, and cross-trait EBVs through `genetic_covariance()`, `residual_covariance()`, `genetic_correlation()`, `residual_correlation()`, `heritability()`, and `breeding_values()`. Missing trait cells are accepted as `NA` and marshal to Julia `NaN`; rank-deficient fixed-effect designs are rejected up front; non-converged fits do not expose `logLik()`/`AIC()`. This is opt-in, dense validation-scale, and partial; t>=2 known-truth recovery, external comparator parity, and long-format/structured covariance grammar remain planned.
* `hsquared()` now fences the reserved multivariate `engine_control$genetic_structure` and `engine_control$rank` fields: `"unstructured"` is accepted for the current opt-in multivariate bridge, while `"diagonal"` was initially guarded and later promoted to the experimental diagonal-G control described above; `"lowrank"`, `"factor_analytic"`, and `rank` still error as planned rather than being silently ignored.
* `formula_status()` and structured-covariance error messages now list the full planned covariance vocabulary (`cov = us()`, `cov = diag()`, `cov = lowrank(K)`, and `cov = fa(K)`) while keeping it planned, not fitted.
* Multivariate `cbind()` responses now require unique, non-empty trait names before fitting, so G/R matrices, per-trait h², EBVs, comparator files, and future wide/long response paths share an unambiguous trait-order contract (#10).
* `G_matrix()` and `R_matrix()` now alias `genetic_covariance()` and `residual_covariance()` for `hsquared_fit` objects, giving multivariate users the familiar G/R matrix names without changing the underlying extractor contract.
* Added a "Reading G matrices" pkgdown article that explains the current G/R matrix extractors, genetic and residual correlations, per-trait h2, cross-trait EBVs, and the boundaries around `P_matrix()`, factor-analytic loadings, and selection-response claims.
* Added a "Genomic prediction" pkgdown article that separates the current opt-in supplied-`Ginv`, marker-built GREML, SNP-BLUP, supplied-`Hinv` single-step, constructed-`Hinv` single-step, animal-only supplied-`Gamma` metafounder, and supplied-`Gamma` `H^Gamma` single-step paths from APY, external metafounder validation, GWAS/QTL/eQTL, and production-comparator work.
* Added a "QTL, GWAS, and eQTL status" pkgdown article that explains the current reserved scan vocabulary, live SNP-BLUP marker effects / descriptive marker variance shares, scale caveats, and validation gates before marker-scan, QTL, GWAS, or eQTL output can be claimed.
* Added an "Inheritance systems roadmap" pkgdown article that gives selfing, clonal, haplodiploid, polyploid, cytoplasmic, imprinting, dominance, epistasis, and custom-kernel examples as planned relationship/precision-kernel work, while keeping current support limited to the v0.1 additive animal model and opt-in standard two-effect slices.
* `prediction_error_variance()` and `reliability()` are now part of the R fitted-object extractor contract. They work for `hsquared_fit` objects containing those result fields (#5, #6).
* Added marker/QTL/eQTL extractor names: `marker_effects()`, `marker_variance_explained()`, `qtl_table()`, `gwas_table()`, `eqtl_table()`, and `lod_scores()`. `marker_effects()` now returns the per-marker effects of an opt-in SNP-BLUP fit, and `marker_variance_explained()` reports descriptive fitted-marker shares for the same path; scan tables and LOD outputs remain output-vocabulary placeholders, with marker-scan, QTL, GWAS, and eQTL fitting still planned (#5, #9).
* Added an explicit supplied-variance Julia bridge target for Henderson MME validation. `hs_control(engine = "julia", engine_control = list(target = "henderson_mme", variance_components = c(sigma_a2 = ..., sigma_e2 = ...)))` calls Julia `henderson_mme()` for tiny validation examples, returning fixed effects, EBVs, fitted values, variance components, and h² without claiming variance-component estimation or production fitting (#6, #7).
* The supplied-variance Henderson MME bridge target now attaches dense validation-path PEV and reliability fields when the sibling `HSquared.jl` checkout exposes applicable `prediction_error_variance()` and `reliability()` methods for `HendersonMMEResult` (#6, #7).
* The experimental local Julia bridge now enriches tiny `hsquared_fit` results with dense validation-path PEV and reliability fields when the sibling `HSquared.jl` checkout exposes `prediction_error_variance()` and `reliability()` (#5, #6).
* The experimental Julia bridge now sends sparse `Matrix::dgCMatrix` random-effect designs through Julia CSC slots instead of densifying `Z` (#6).
* Added an internal tiny animal-model validation fixture that pins R payload ordering, sparse `Z` construction, and live Julia `pedigree_inverse()` agreement for a three-animal Henderson-style pedigree when a sibling `HSquared.jl` checkout is available (#7).
* Added an optional Mrode9/nadiv pedigree-Ainv comparator fixture. When `nadiv` and a sibling `HSquared.jl` checkout are available, local tests compare Julia `pedigree_inverse()` with `nadiv::makeAinv()` for the Mrode9 pedigree (#7).
* Added a Mrode-style supplied-variance validation fixture that pins Ainv, fixed effects, EBVs, fitted values, PEV, reliability, h², ML log-likelihood, and dense/sparse REML log-likelihood against R reference calculations and the sibling `HSquared.jl` checkout when available, without claiming variance-component estimation or full fitted Mrode validation (#7).
* Added an internal supplied-variance Henderson mixed-model-equation validation fixture that compares R reference fixed effects, EBVs, fitted values, and h2 with Julia `henderson_mme()` when a sibling `HSquared.jl` checkout is available (#7).
* Added a tiny supplied-variance REML likelihood validation fixture. When a sibling `HSquared.jl` checkout is available, optional local tests compare Julia dense REML, sparse REML, and ML hand-check targets without claiming sparse optimization or fitted Mrode output validation (#7).
* Expanded `hs_control()` to preserve planned backend and accelerator vocabulary for CPU threads, CUDA, AMDGPU, Metal, and oneAPI. These are control-surface placeholders only; GPU execution remains planned (#3).
* Added `backend_info()` so users and developers can inspect planned backend names while seeing that backend execution is not available yet.
* Added `validation_status()` so users and developers can inspect current validation atoms, planned comparator lanes, and claim boundaries from R.
* `data_status()` now gives users a direct diagnostic view of `hs_data()` component presence, ID overlap, pedigree coverage, and marker-map/genotype-marker alignment status. It is a status helper only and does not fit models (#8).
* Added `formula_status()` so users and developers can inspect parsed, reserved, and planned formula grammar without reading the full roadmap.
* Added `model_spec()` so users and developers can preview the parsed v0.1 animal-model contract, fixed-effect design columns, sparse animal-effect design dimensions, normalized pedigree ordering, and Julia targets without fitting a model.
* Added inert planned formula markers for `genomic()`, `single_step()`, `markers()`, `marker_scan()`, and `qtl_scan()`. The parser now rejects these terms with explicit planned-not-implemented errors instead of treating them as fixed effects.
* Added inert planned formula markers for `permanent()`, `common_env()`, `maternal_genetic()`, `maternal_env()`, `paternal_genetic()`, `paternal_env()`, `cytoplasmic()`, `imprinting()`, `dominance()`, `epistasis()`, `relmat()`, and `precision()`. They reserve Phase 2+ vocabulary only and currently abort as planned, not implemented.
* Added a pkgdown formula grammar roadmap article that separates parsed v0.1 syntax from planned quantitative-genetic, genomic, multivariate, and inheritance syntax.
* Expanded the genomics/QTL/GLLVM/accelerator design plan and pkgdown roadmap with a source-backed CPU/GPU strategy, QTL/eQTL path, sibling-package lessons, and explicit evidence gates.
* Added an internal R-to-Julia bridge payload builder for the v0.1 animal-model contract. It creates `y`, `X`, sparse `Z`, normalized pedigree metadata, and the validated Julia `animal_model_spec()` target used by the experimental Julia engine (#6).
* Added a local-only experimental JuliaCall smoke path for the tiny v0.1 payload when a sibling `HSquared.jl` checkout is available. This validates bridge shape against Julia `pedigree_inverse()` and `fit_animal_model()` but is not yet the public `hsquared()` fitting path (#6).
* Added the first `hsquared_fit` object and extractor contract, including `variance_components()`, `heritability()`, `breeding_values()`, `fixef()`, `ranef()`, `logLik()`, `AIC()`, `predict()`, and `summary()` methods over internal fit objects. These are contract plumbing only until the Julia engine returns real fits (#5).
* Added `hs_data()` as a lightweight R data container for phenotype, pedigree, genotype, marker, expression, annotation, and environment inputs. It records ID maps for future integrated genomic/QTL/eQTL workflows, but does not fit models (#8).
* `hs_data()` now accepts `annotation_id` to check expression feature columns against annotation rows, and `summary()`/`data_status()` report annotation-feature diagnostics without fitting eQTL or omics models (#8).
* `hs_data()` now accepts `environment_id` to check environment/covariate metadata coverage against phenotype records, and `summary()`/`data_status()` report environment-key diagnostics without constructing environmental model terms (#8).
* `summary(hs_data(...))` and `data_status()` now report expression row counts, expression ID counts, feature counts, unnamed feature columns, duplicate feature IDs, and expression component type without fitting eQTL or omics models (#8).
* `summary(hs_data(...))` and `data_status()` now report genotype row counts, genotype ID counts, marker-column counts, unnamed marker columns, duplicate marker columns, missing genotype value counts, and genotype component type without fitting genomic, marker-scan, or QTL models (#8).
* `summary(hs_data(...))` now includes pedigree coverage and parent-link diagnostics when a pedigree component is supplied (#8).
* `summary(hs_data(...))` now includes an ID overlap table with phenotype, pedigree, genotype, expression, and mismatch counts (#8).
* `summary(hs_data(...))` now includes marker-map and genotype-marker alignment diagnostics when marker or genotype components are supplied (#8).
* `hs_data()` now validates supplied marker maps for marker ID, chromosome, and non-negative numeric position columns. This is metadata validation only; genomic and QTL/eQTL fitting remain planned (#8).
* `hs_data()` now checks that genotype marker column names match marker-map IDs exactly when both `genotypes` and `markers` are supplied (#8).
* `model_spec()` and `hsquared()` can now use an `hs_data()` object directly for the v0.1 parser, reading model variables from `phenotypes` and resolving formula components such as `pedigree = pedigree` from the bundle (#8).
* `model_spec()` and `hsquared()` now allow `animal(1 | id)` to use the pedigree stored in `data = hs_data(..., pedigree = ped)`, while ordinary data frames still require explicit `pedigree = ped` (#4, #8).
* Added Phase 0 project operating documentation, an honest placeholder `hsquared()` entry point, and `hs_control()` for planned engine controls.
