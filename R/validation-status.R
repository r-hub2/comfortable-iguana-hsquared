#' Inspect validation evidence status
#'
#' `r lifecycle::badge("experimental")`
#'
#' `validation_status()` reports the current validation atoms and planned
#' comparator lanes for `hsquared`. It is a **developer evidence table**, not
#' the user-facing list of what you may report. For that list, read
#' [Can I fit and report this?](
#' https://itchyshin.github.io/hsquared/articles/current-limits.html).
#' This function does not run validation checks, fit models, or promote any
#' capability to working status.
#'
#' # This table is not the full list of covered routes
#'
#' The rows here are validation *atoms*, grouped by evidence chain rather than
#' by public modelling route, so the `covered` row count is not the number of
#' covered routes. Two capabilities that count as publicly covered have no row
#' of their own yet:
#'
#' - the opt-in reaction-norm model, k = 2:
#'   `animal(rr(x, order = 2) | id, pedigree = ped)` with
#'   `engine_control = list(target = "random_regression")`;
#' - the opt-in correlated direct-maternal 2x2 G model:
#'   `animal(1 | id, pedigree = ped) + maternal_genetic(1 | dam)` with
#'   `engine_control = list(target = "direct_maternal")`.
#'
#' Both are covered at validation scale, opt-in, REML-only, and dense
#' (n <= 1000 or so). Their scope and reporting limits are recorded in the
#' repository capability ledger (`docs/design/capability-status.md`) and, for
#' readers, in [Can I fit and report this?](
#' https://itchyshin.github.io/hsquared/articles/current-limits.html). Read
#' this table together with that enumeration: an absent row is not an absent
#' capability, and this table alone does not tell you everything you may
#' report.
#'
#' # Bridge production fences (0.9)
#'
#' The R-to-Julia bridge is **live at validation scale** for covered default and
#' opt-in targets, but **not production-complete**. Infrastructure capabilities
#' (`R-to-Julia bridge payload`, `opt-in experimental Julia engine`) stay
#' **partial** even when individual targets are covered. Targets split into:
#' default live (`engine = "fit"`), opt-in live (`engine_control$target`),
#' **PATH_ONLY** smoke (C1-ext intervals — not a fit bridge), and **blocked**
#' planned grammar (for example `genetic_structure = "factor_analytic"`).
#' `payload_v2` routing is limited to `direct_maternal` and `multi_effect`.
#' See `docs/design/45-bridge-production-fences-DRAFT.md`.
#'
#' # Capability identifiers versus labels
#'
#' `capability` is a **stable identifier**. Dated evidence records -- comparator
#' runs under `docs/dev-log/comparator-runs/` and check-log entries -- cite these
#' strings verbatim to name the row they report against, so an identifier is not
#' rewritten once evidence points at it. Look rows up by `capability`.
#'
#' `capability_label` is the **current wording**, and is what reader-facing
#' surfaces print. The two differ only where an identifier's wording has been
#' overtaken by a later change; today that is the multivariate row, whose
#' `cbind()` route became default rather than opt-in.
#'
#' @return A data frame of validation status records with class
#'   `"hs_validation_status"`.
#' @examples
#' validation_status()
#' @export
validation_status <- function() {
  out <- data.frame(
    capability = hs_validation_status_capabilities(),
    capability_label = hs_validation_status_labels(),
    phase = hs_validation_status_phases(),
    status = hs_validation_status_status(),
    evidence = hs_validation_status_evidence(),
    claim_boundary = hs_validation_status_boundaries(),
    stringsAsFactors = FALSE
  )
  class(out) <- c("hs_validation_status", class(out))
  out
}

#' @export
print.hs_validation_status <- function(x, ...) {
  cat("<hs_validation_status>\n")
  cat("  validation: status table only; checks are run by tests and CI\n")
  cat("  public claims: developer evidence table, not the user-facing list\n")
  cat("  of what you may report; see Can I fit and report this?\n")
  cat("  this table is not the full list of covered routes -- opt-in\n")
  cat("  random_regression (k = 2) and direct_maternal are covered with no\n")
  cat("  row yet; see the current-limits article\n")
  cat("  shown below: `capability_label` (current wording); look rows up by\n")
  cat("  the stable `capability` id that dated evidence records cite\n")
  out <- x
  class(out) <- setdiff(class(out), "hs_validation_status")
  print.data.frame(
    out[c("capability_label", "phase", "status")],
    row.names = FALSE
  )
  invisible(x)
}

# Capability identifiers are STABLE: dated comparator-run reports and check-log
# entries cite them verbatim to name the row they report against, so renaming an
# id would orphan its evidence. When an id's wording is overtaken by a later
# change, the correction is recorded here as a display alias instead of a
# rename. See docs/dev-log/decisions.md, "2026-09-02: Capability Ids Are
# Historical; Labels Carry Current Wording".
hs_validation_status_label_overrides <- function() {
  c(
    # MV-4 made the `cbind()` route default rather than opt-in, but three dated
    # records (two comparator-run reports, one check-log entry) already cite the
    # "(opt-in)" id. G10 (2026-09-02) covers the R-public t=2 unstructured
    # surface; the historical id is unchanged.
    "experimental multivariate REML estimator (opt-in)" = "experimental multivariate REML estimator (default route)"
  )
}

hs_validation_status_labels <- function() {
  ids <- hs_validation_status_capabilities()
  overrides <- hs_validation_status_label_overrides()
  hit <- match(ids, names(overrides))
  unname(ifelse(is.na(hit), ids, overrides[hit]))
}

hs_validation_status_capabilities <- function() {
  c(
    "tiny deterministic Ainv fixture",
    "Mrode9 pedigree Ainv comparator",
    "supplied-variance Henderson MME fixture",
    "sparse REML likelihood identity",
    "Mrode-style supplied-variance outputs",
    "experimental sparse REML estimator (opt-in)",
    "experimental repeatability estimator (opt-in)",
    "two-effect / arbitrary-N independent-effect estimator (opt-in; covered: common-env + (1|g) iid / A2=I; experimental: maternal / A2=pedigree)",
    "experimental supplied-relationship estimator (opt-in: genomic, single-step)",
    "experimental SNP-BLUP marker-effect model (opt-in; supplied-variance or REML-estimated)",
    "experimental multivariate REML estimator (opt-in)",
    "univariate Gaussian animal-model fit (default path, AI-REML)",
    "external published-REML recovery (gryphon, R reference)",
    "known-truth DGP variance-component recovery (R reference)",
    "Mrode fitted animal-model outputs",
    "ASReml comparison policy",
    "BLUPF90/DMU/WOMBAT comparison policy",
    "XSim simulation truth",
    "genomic and QTL/eQTL validation",
    "GLLVM-style multivariate validation",
    "CPU/GPU backend comparison"
  )
}

hs_validation_status_phases <- function() {
  c(
    rep("Phase 1", 6L),
    rep("Phase 2", 2L),
    rep("Phase 5", 2L),
    "Phase 3",
    rep("Phase 1", 6L),
    rep("Phase 5+", 2L),
    "Phase 6",
    "Phase 7+"
  )
}

hs_validation_status_status <- function() {
  c(
    rep("partial", 7L), # positions 1-7
    "covered", # position 8 = two-effect / arbitrary-N independent-effect estimator (COMMON-ENV + (1|g) iid / A2=I covered; maternal experimental)
    "covered", # position 9 = genomic GREML (0.7 covered validation-scale; opt-in)
    "partial", # position 10 = SNP-BLUP
    "covered", # position 11 = multivariate REML (G10 t=2 unstructured R-public)
    rep("covered", 3L),
    rep("planned", 7L)
  )
}

hs_validation_status_evidence <- function() {
  c(
    paste(
      "Local tests pin R payload ordering, sparse Z, and live Julia",
      "pedigree_inverse() agreement when HSquared.jl is available."
    ),
    paste(
      "Optional local tests compare Julia pedigree_inverse() with",
      "nadiv::makeAinv() for nadiv::Mrode9 when optional dependencies are",
      "available."
    ),
    paste(
      "Local R reference solve and optional live Julia henderson_mme()",
      "comparison for fixed effects, EBVs, fitted values, h2, and optional",
      "dense validation-path PEV/reliability."
    ),
    paste(
      "Optional local tests compare Julia gaussian_loglik() dense REML and",
      "sparse_reml_loglik() on a tiny three-founder fixture at supplied",
      "variance components, with ML and REML hand-check targets."
    ),
    paste(
      "Local R reference tests and optional live Julia tests pin a",
      "Mrode-style supplied-variance fixture for Ainv, fixed effects, EBVs,",
      "fitted values, PEV, reliability, h2, ML logLik, and REML logLik; the",
      "reference EBVs are additionally pinned to the published Mrode Example",
      "3.1 textbook digits (Mrode 2014, p.39) in a CI-runnable anchor."
    ),
    paste(
      "Pure-R control/validator tests plus a skip-guarded live test running",
      "Julia fit_sparse_reml() through the opt-in target = \"sparse_reml\"",
      "bridge; checks positive estimated variances, finite REML logLik, and",
      "h2 in (0, 1) on the Mrode fixture; fit provenance tagged",
      "variance_components_source = \"estimated_sparse_reml\". Cross-checked",
      "against the dense REML optimizer, an independent pure-R REML optimizer,",
      "and the external pedigreemm package (at least as good by REML logLik)."
    ),
    paste(
      "Pure-R control/validator tests plus a skip-guarded live test running",
      "Julia fit_repeatability_reml() through the opt-in target =",
      "\"repeatability\" bridge on a repeated-records fixture; checks three",
      "positive estimated variance components (animal, permanent, residual),",
      "repeatability in (0, 1), heritability in [0, 1), and permanent-",
      "environment effects; fit provenance tagged",
      "variance_components_source = \"estimated_repeatability_reml\". The",
      "permanent-environment effect shares the animal incidence (A2 = I)."
    ),
    paste(
      "Pure-R control/validator tests plus skip-guarded live tests running Julia",
      "fit_two_effect_reml() through the opt-in target = \"two_effect\" bridge on",
      "`animal(1 | id) + common_env(1 | group)` (IID environment, A2 = I) and on",
      "`animal(1 | id) + maternal_genetic(1 | dam)` (maternal genetic effect, A2 =",
      "pedigree A); each checks three positive estimated variance components",
      "(animal, second effect, residual) and the estimated_two_effect_reml",
      "provenance; the common-environment leg additionally checks heritability in",
      "[0, 1) and the common-environment effect predictions. Fit provenance tagged",
      "variance_components_source = \"estimated_two_effect_reml\"."
    ),
    paste(
      "Opt-in experimental only; Julia V2-SSHINV engine-covered does not promote",
      "R-public covered. AGHmatrix/BLUPF90 external comparators owed. Pure-R",
      "control/validator tests plus skip-guarded live tests running Julia",
      "fit_ai_reml() on a user-supplied relationship inverse through the opt-in",
      "bridge: target = \"genomic\" on `genomic(1 | id, Ginv = Ginv)` (a genomic",
      "relationship inverse) and target = \"single_step\" on",
      "`single_step(1 | id, Hinv = Hinv)` (a single-step relationship inverse);",
      "each checks two positive estimated variance components (effect, residual)",
      "and breeding values (the supplied-Ginv genomic leg also checks heritability",
      "in (0, 1)); fit",
      "provenance tagged variance_components_source =",
      "\"estimated_<genomic|single_step>_ai_reml\". For genomic, a skip-guarded",
      "live test also fits from a raw marker matrix (`genomic(1 | id, markers =",
      "M)`): the engine builds G (genomic_relationship_matrix) and its regularized",
      "inverse, then fits. The R single_step bridge also has a construction",
      "target (`single_step(1 | id, pedigree = ped, markers = M)`,",
      "target = \"single_step_construct\"): the engine builds Ainv + dense A from",
      "the pedigree, builds G from genotyped-subset markers, assembles H^-1,",
      "and fits by REML. Live tests cover marker-row reorder invariance,",
      "all-pedigree GEBV labels including ungenotyped animals, differs-from-",
      "pedigree-model behavior, ridge handling, and hs_data() shorthand. The",
      "animal-only supplied-Gamma metafounder path (`metafounder(..., group =,",
      "Gamma =)`, target = \"metafounder\") calls metafounder_animal_model with",
      "supplied sigma_a2/sigma_e2; live tests pin Gamma = 0 reduction to the",
      "ordinary Henderson MME supplied-variance path and nonzero-Gamma prediction",
      "sensitivity with stable labels/dimensions. The",
      "supplied-Gamma H^Gamma path (`single_step(..., group =, Gamma =)`,",
      "target = \"metafounder_single_step\") calls fit_metafounder_single_step_reml;",
      "live tests pin Gamma = 0 reduction to ordinary single-step construction",
      "and nonzero-Gamma prediction sensitivity with stable labels/dimensions."
    ),
    paste(
      "Pure-R control/validator tests plus skip-guarded live tests running Julia",
      "through the opt-in target = \"snp_blup\" bridge on",
      "`genomic(1 | id, markers = M)`: with supplied variances (sigma_g2,",
      "sigma_e2) it runs fit_snp_blup() (provenance variance_components_source =",
      "\"supplied\"); with variances omitted it runs fit_snp_blup_reml(), which",
      "ESTIMATES sigma_g2/sigma_e2 by REML (provenance \"estimated_snp_blup_reml\",",
      "live-verified to match a direct fit_snp_blup_reml call). The engine centers",
      "the markers, solves the RR-BLUP/SNP-BLUP marker model, and returns",
      "per-marker effects (marker_effects()), per-individual genomic breeding",
      "values, and fixed effects. Mirrors the twin V2-SNPBLUP gate (partial),",
      "whose pinned property is the GBLUP<->SNP-BLUP GEBV equivalence."
    ),
    paste(
      "Pure-R parser/validator tests plus skip-guarded live bridge tests for the",
      "Julia fit_multivariate_reml() target: `cbind(y1, y2) ~ animal(1 | id,",
      "pedigree = ped)` builds an NA-preserving Y matrix, sends missing trait",
      "cells as Julia NaN, and returns G0/R0 covariance matrices, genetic and",
      "residual correlations, per-trait h2, cross-trait EBVs, and convergence",
      "diagnostics; fit provenance tagged variance_components_source =",
      "\"estimated_multivariate_reml\". The shared Phase 4 fixture pins R",
      "payload/extractor parity against serialized Julia targets, and an",
      "optional in-suite sommer comparator checks the same fixture's G0,",
      "diag(R0), and diagonal-target h2. Two reproducible R-lane studies now",
      "extend that evidence without promoting the row: a 100-replicate",
      "cold-start t=2 known-truth recovery study has all six G0/R0 elements,",
      "genetic correlation, and both per-trait h2 within bias +/- 2*MCSE",
      "(100/100 converged; EBV accuracy 0.79/0.74), and a full-unstructured",
      "residual sommer comparator reproduces the serialized phase4_multitrait",
      "target's G0/R0/beta/h2/EBV to <= 8e-5 while recovering the off-diagonal",
      "R0 residual covariance. A Bayesian MCMCglmm agreement probe puts the",
      "serialized Julia target inside 95% HPD intervals for all 8 covariance",
      "elements, all 4 fixed effects, and both per-trait h2 values, with",
      "posterior-mean EBV correlations > 0.9997; this is explicitly not a",
      "same-estimand REML comparator. A pure-R CI anchor also reproduces the",
      "published Mrode Example 5.1 multiple-trait supplied-G0/R0 BLUP/MME",
      "fixed effects and animal BLUPs from the LUKE/Masuda reproductions.",
      "NO-ANCHOR DISCLOSURE (Standard-Tier covered-flip gate item 2): Mrode",
      "Example 5.1 is a SUPPLIED-G0/R0 BLUP/MME anchor - it pins fixed effects",
      "and animal BLUPs GIVEN known covariances - so it does NOT anchor the",
      "ESTIMATED G0/R0 that this row's REML claim is about. No published",
      "textbook anchor for estimated multivariate G0/R0 exists; the",
      "estimated-covariance evidence is comparator- and recovery-based",
      "(sommer 4.4.5, blupf90+ 2.60, cold-start recovery), not",
      "textbook-anchored. Any covered flip carries this disclosure, never an",
      "implied Mrode anchor (precedent: docs/design/43-genomic-greml-g0.md).",
      "Mirrors the twin engine multivariate REML validation row (partial)."
    ),
    paste(
      "The default `hsquared()` control fits the v0.1 Gaussian animal model by",
      "average-information REML (`fit_ai_reml`) on the sparse MME through the",
      "bridge. A skip-guarded live test fits the contract through the default",
      "control and checks positive estimated variances, finite REML logLik, and",
      "h2 in (0, 1) on the Mrode fixture; fit provenance tagged",
      "variance_components_source = \"estimated_ai_reml\". The estimator reaches",
      "the same REML optimum as the sparse REML optimizer, and its known-truth",
      "and external-anchor recovery are evidenced in the two rows below. Mirrors",
      "the twin V1-AI-REML gate (covered)."
    ),
    paste(
      "A skip-guarded test recovers the published gryphon birth-weight REML",
      "variance components and h2 (Wilson et al. 2010) with hsquared's",
      "independent pure-R REML reference, and agrees with the external sommer",
      "package within the signed-off band (VC ~1-2%, h2 ~0.01-0.02), on the",
      "gryphon dataset (CRAN package enhancer). The Julia engine (fit_sparse_reml",
      "and fit_ai_reml) also recovers the published estimates within the",
      "signed-off band via supplied A_gryphon (the engine correctly rejects the",
      "pathological raw pedigree); engine-vs-pure-R agreement is to machine",
      "precision. Gryphon is the maintainer (2026-06-13) signed-off V1-MRODE-FIT",
      "anchor and sommer the V1-COMPARATORS comparator."
    ),
    paste(
      "An ADEMP recovery study (data-raw/dgp-recovery-study.R) simulates from",
      "known variance components over a clean pedigree; the engine (ai_reml) is",
      "near-unbiased (0 within bias +/- 2*MCSE for sigma_a2, sigma_e2, h2 over",
      "120 reps, 100% converged), EBVs track true breeding values (acc ~0.74),",
      "and the engine matches the independent pure-R reference to machine",
      "precision. Recovery is near-unbiased across an h2 grid (0.2/0.4/0.6); the",
      "near-boundary cell (h2 = 0.1) shows mild upward bias, 94% convergence, and",
      "5% boundary pinning. Recovery also holds for the contract model with a",
      "fixed effect (y ~ x + animal): h2 near-unbiased and the fixed effect",
      "recovered (b_x 0.99 vs 1.0). Meets the maintainer-signed-off item-3",
      "thresholds (>= 100 reps, 0 within bias +/- 2*MCSE, mean cor(EBV, true) >=",
      "0.5). A skip-guarded pure-R regression test guards a small-N case."
    ),
    "None yet.",
    "None yet.",
    "None yet.",
    "None yet.",
    "None yet.",
    "None yet.",
    "None yet."
  )
}

hs_validation_status_boundaries <- function() {
  c(
    "Internal tiny Ainv atom only; not production fitting or Mrode coverage.",
    "Pedigree inverse comparator only; not fitted Mrode outputs.",
    paste(
      "Supplied-variance BLUP/MME only; not variance-component estimation,",
      "AI-REML, production sparse reliability, or production sparse fitting."
    ),
    paste(
      "Supplied-variance likelihood identity only; not a sparse optimizer,",
      "AI-REML, fitted Mrode output validation, or production sparse fitting."
    ),
    paste(
      "Supplied-variance validation fixture only; not variance-component",
      "estimation, ASReml parity, or general fitted Mrode output coverage."
    ),
    paste(
      "Experimental opt-in path only; Julia-owned estimator that R surfaces;",
      "gated on twin validation_status; not the default, not production sparse",
      "fitting, AI-REML, or ASReml parity."
    ),
    paste(
      "Experimental opt-in path only; Julia-owned REML-only repeatability",
      "optimizer that R surfaces; mirrors the twin V3-REPEAT-REML gate",
      "(partial). Not the default, not ML, not production fitting, and not a",
      "comparator-validated or known-truth-recovery claim. The additive (Va) and",
      "permanent-environment (Vpe) variances are only identifiable with repeated",
      "records per individual."
    ),
    paste(
      "Julia-owned REML-only two-effect optimizer that R surfaces (opt-in, NOT the",
      "default engine='fit' path). The COMMON-ENVIRONMENT leg (additive animal-A +",
      "IID common-environment via common_env(), A2 = I) is COVERED (experimental,",
      "validation-scale, opt-in) - mirrors the twin V3-TWOEFFECT-REML covered gate:",
      "a pre-declared 48-seed bias/MCSE recovery gate PASSED + a blupf90+ same-estimand",
      "REML comparator agrees ~1e-5 (sommer cross-check ~2e-5). The MATERNAL genetic",
      "leg (maternal_genetic(), A2 = pedigree A) uses the SAME estimator with exact",
      "live parity but STAYS EXPERIMENTAL - its own recovery gate + comparator on the",
      "maternal-A2 design are owed. Not ML, not production sparse fitting; the h2/c2",
      "interval is asymptotic/delta-method and NOT coverage-calibrated; the correlated",
      "direct-maternal (2x2 G) model is a SEPARATE opt-in target (target='direct_maternal'),",
      "now COVERED at validation scale (mirrors the twin V4-DIRECT-MATERNAL: pre-declared",
      "48-seed gate PASSED + sommer covm() same-estimand comparator AGREE; Willham",
      "labelled-triple h2, dense n<=1000). The arbitrary-N generalization to",
      "independent (1|g) i.i.d. effects (target='multi_effect') is COVERED on the same",
      "V3-NEFFECT-REML engine gate with exact live R-Julia parity; the animal-block ratio",
      "is narrow-sense h2, other blocks are variance-explained proportions (not",
      "heritabilities), intervals asymptotic/uncalibrated; INDEPENDENT effects only (NOT",
      "correlated / random-regression / non-Gaussian). V8.6 SCALE opt-in:",
      "engine_control scale_method='auto' (default 'dense') routes large problems through",
      "the engine's EXPERIMENTAL matrix-free Monte-Carlo fit; validation scale uses the",
      "sparse-exact estimator (reduces to the covered dense result, live-bridge verified).",
      "The COVERED claim is unchanged (validation-scale, exact); the large-scale matrix-free",
      "path is opt-in experimental."
    ),
    paste(
      "Covered at VALIDATION scale for opt-in genomic GREML (0.7 / owner #7);",
      "Julia-owned REML estimator (fit_ai_reml / fit_gblup_reml on a",
      "relationship-inverse or marker-built Q_lambda) that R surfaces; mirrors",
      "the twin V2-GREML covered gate for the supplied-Ginv estimator, with",
      "marker==Q identity + Totoro exact-G comparator + design-53 SUPERSEDE.",
      "Single-step / metafounder construction paths in this combined row remain",
      "PARTIAL (V2-SSHINV). Genomic GREML accepts a supplied Ginv or a marker",
      "matrix (engine-built VanRaden1 G + ridge 0.01); public estimand is",
      "genomic_variance_ratio on K_lambda (design-51). Default activation,",
      "APY, SNP-BLUP, intervals, and field-panel robustness are out of the",
      "covered claim. Experimental label retained; not the default route, not",
      "ML, not production. NO-ANCHOR: no clean Mrode genomic-h2 pin.",
      "public_covered_count is 7."
    ),
    paste(
      "Experimental opt-in path only; Julia-owned VanRaden method-1 marker model",
      "that R surfaces; mirrors the twin V2-SNPBLUP gate (partial). The user may",
      "supply sigma_g2 and sigma_e2 (fit_snp_blup), or omit them to ESTIMATE them",
      "by REML from the markers (fit_snp_blup_reml). Weighted/Bayesian marker",
      "priors, low-rank m>>n Woodbury solves, and JWAS/sommer/BLUPF90 comparator",
      "parity are planned. Not the default, not comparator-validated."
    ),
    paste(
      "Covered at VALIDATION scale (2026-09-02 maintainer sign-off); Julia-owned dense",
      "multivariate REML estimator that R surfaces; mirrors the twin engine",
      "multivariate REML validation row. ROUTING: a `cbind()` Gaussian response with an `animal()`",
      "term auto-routes to this fitter on the DEFAULT path; no",
      "engine/target argument is required, and the explicit `engine =",
      "\"julia\", target = \"multivariate\"` spelling still works. Covered numeric",
      "claim is scoped to k = 2 unstructured; k >= 3 stays",
      "parseable-and-fittable-but-experimental; genetic_structure =",
      "\"diagonal\" stays experimental at 0.6; \"factor_analytic\" and \"lowrank\"",
      "error as planned on the R bridge (rotation fence). `cbind()` responses with missing",
      "trait cells are supported, but this is REML-only and animal-model-only.",
      "The R lane has cold-start recovery and one reproduced full-unstructured sommer",
      "comparator leg plus a published Mrode-style supplied-variance BLUP/MME",
      "anchor plus a Bayesian MCMCglmm agreement probe. The MCMCglmm leg is",
      "not same-estimand REML parity. The engine multivariate REML validation",
      "row is covered at validation scale (HSquared.jl#161). This R public surface is COVERED at",
      "validation scale after maintainer sign-off (Shinichi, 2026-09-02);",
      "public_covered_count is 7. Twin C8 broader-DGP",
      "confirm is banked (DRAC job 47925486; 16x500; 14/16 pass; fails only",
      "rg_090_rec1/rg_095_rec1; base_inside clean) as characterization; R-lane",
      "MV-5 is SUPERSEDED (A25). The engine blupf90+ same-estimand REML",
      "comparator is discharged. A26 R<->engine element-wise parity is",
      "DISCHARGED LOCALLY (run at a recorded commit, within predeclared",
      "tolerances) but is NOT CI-backed: Tier-1 parity CI is not yet live (no",
      "workflow provisions Julia and the parity legs skip on a Julia-free",
      "runner), so read it as local evidence, not a CI-verified claim.",
      "Darwin A27 SIGNED; Rose A29 CLEAN. Intervals remain uncalibrated",
      "(design-41 section 3 #7 deferred). NO-ANCHOR DISCLOSURE (gate item 2): the",
      "published Mrode Example 5.1 anchor is SUPPLIED-covariance BLUP/MME and",
      "does NOT anchor an ESTIMATED G0/R0; no published textbook anchor for",
      "estimated multivariate G0/R0 exists. MV-1 sommer requires Suggests on",
      "NOT_CRAN (hs_require_suggests; CRAN still skips). The Julia engine",
      "currently inverts Ainv internally, so deep-inbreeding or",
      "high-condition-number pedigrees remain a twin-side hardening item.",
      "Not production, not ASReml-style multivariate fitting."
    ),
    paste(
      "Univariate Gaussian animal model only (single additive genetic effect);",
      "REML only (`REML = FALSE` is rejected on the default fit and on",
      "`engine = \"validate\"`; use `REML = TRUE`). The covered claim is",
      "this REML estimator, not ML. Genomic, repeatability,",
      "two-effect, marker-effect, and non-Gaussian",
      "(poisson(log)/binomial(logit), Laplace marginal likelihood or variational ELBO)",
      "fitting are separate",
      "opt-in experimental targets, not the default. A `cbind()` multivariate",
      "response also routes on the default path, but it is a separate,",
      "covered (validation-scale) sibling row and is NOT part of this univariate covered claim.",
      "Mirrors the twin-owned",
      "V1-AI-REML gate",
      "(covered); not ASReml multi-trait parity."
    ),
    paste(
      "Published-anchor recovery within the maintainer-signed-off band (not",
      "bit-exact); backs the twin-owned V1-MRODE-FIT gate (covered_external).",
      "Gryphon is teaching data; confirm headline numbers before any new",
      "promotion."
    ),
    paste(
      "Known-truth recovery evidence for the default AI-REML estimator; R-lane",
      "ADEMP study via the read-only bridge across an h2 grid (0.2/0.4/0.6) plus",
      "a near-boundary cell (h2 = 0.1) and a fixed-effect model; backs the",
      "twin-owned V1-AI-REML gate (covered). No interval-coverage or",
      "production-robustness claim beyond the studied grid."
    ),
    "Planned; no Mrode fitted-output claim.",
    "Planned; no ASReml parity claim.",
    "Planned; no external production-software parity claim.",
    "Planned; no simulation-recovery claim.",
    "Planned; no genomic, QTL, eQTL, or marker-scan claim.",
    "Planned; no GLLVM-style animal-model validation claim.",
    "Planned; no backend execution, benchmark, or speedup claim."
  )
}
