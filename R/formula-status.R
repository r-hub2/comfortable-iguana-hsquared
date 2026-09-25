#' Inspect formula grammar status
#'
#' `r lifecycle::badge("experimental")`
#'
#' `formula_status()` reports which pieces of the planned `hsquared()` formula
#' language are parsed today, reserved as syntax markers, or still roadmap-only.
#' It is a status table, not a model-fitting helper.
#'
#' The printed header is derived from the rows being printed, so it cannot
#' lag the table. Covered-versus-experimental fences live in
#' `$current_behavior` and are not dumped into the default print.
#'
#' @return A data frame of formula grammar records with class
#'   `"hs_formula_status"`.
#' @examples
#' formula_status()
#' @export
formula_status <- function() {
  cols <- list(
    term = hs_formula_status_terms(),
    category = hs_formula_status_categories(),
    phase = hs_formula_status_phases(),
    syntax_status = hs_formula_status_syntax(),
    fitting_status = hs_formula_status_fitting(),
    current_behavior = hs_formula_status_behavior()
  )
  n <- lengths(cols)
  if (length(unique(n)) != 1L) {
    stop(
      "formula_status() helper vectors are different lengths: ",
      paste(sprintf("%s=%d", names(n), n), collapse = ", "),
      call. = FALSE
    )
  }
  out <- as.data.frame(cols, stringsAsFactors = FALSE)
  class(out) <- c("hs_formula_status", class(out))
  out
}

#' @export
print.hs_formula_status <- function(x, ...) {
  hs_formula_status_print_header(x)
  out <- x
  class(out) <- setdiff(class(out), "hs_formula_status")
  display_cols <- intersect(
    c("term", "syntax_status", "fitting_status"),
    names(out)
  )
  if (length(display_cols) == 0L) {
    display_cols <- names(out)
  }
  print.data.frame(
    out[display_cols],
    row.names = FALSE
  )
  invisible(x)
}

hs_formula_status_short_term <- function(terms) {
  vapply(
    as.character(terms),
    function(term) {
      if (startsWith(term, "cbind(")) {
        return("cbind()")
      }
      if (startsWith(term, "animal(rr(")) {
        return("animal(rr())")
      }
      if (grepl("^animal\\(1 \\|", term)) {
        return("animal()")
      }
      if (grepl("^\\(1 \\|", term)) {
        return("(1 | group)")
      }
      if (grepl("^missing = miss_control", term)) {
        return("miss_control()")
      }
      if (startsWith(term, "mi(")) {
        return("mi()")
      }
      head <- sub("\\(.*", "", term)
      if (nzchar(head) && !identical(head, term)) {
        return(paste0(head, "()"))
      }
      term
    },
    character(1L),
    USE.NAMES = FALSE
  )
}

hs_formula_status_cat_list <- function(label, items) {
  if (!length(items)) {
    return(invisible())
  }
  body <- paste(unique(items), collapse = "; ")
  prefix <- paste0("  ", label, ": ")
  if ((nchar(prefix) + nchar(body)) <= 78L) {
    cat(prefix, body, "\n", sep = "")
    return(invisible())
  }
  cat(prefix, "\n", sep = "")
  cat(
    paste(strwrap(body, width = 70L, prefix = "    "), collapse = "\n"),
    "\n",
    sep = ""
  )
  invisible()
}

hs_formula_status_print_header <- function(x) {
  cat("<hs_formula_status>\n")
  if (all(c("term", "fitting_status") %in% names(x))) {
    fit <- x$fitting_status
    default_idx <- grepl("default", fit, fixed = TRUE)
    opt_idx <- grepl("opt-in", fit, fixed = TRUE)
    if (any(default_idx)) {
      hs_formula_status_cat_list(
        "default",
        hs_formula_status_short_term(x$term[default_idx])
      )
    }
    if (any(opt_idx)) {
      hs_formula_status_cat_list(
        "opt-in",
        hs_formula_status_short_term(x$term[opt_idx])
      )
    }
  }
  if ("syntax_status" %in% names(x)) {
    n_reserved <- sum(x$syntax_status == "reserved", na.rm = TRUE)
    n_planned <- sum(x$syntax_status == "planned", na.rm = TRUE)
    if (n_reserved + n_planned > 0L) {
      cat(
        "  reserved/planned: ",
        n_reserved,
        " reserved, ",
        n_planned,
        " planned (error before fitting)\n",
        sep = ""
      )
    }
  }
  if ("current_behavior" %in% names(x)) {
    cat("  fences: $current_behavior\n")
  }
  invisible()
}

hs_formula_status_terms <- function() {
  c(
    "animal(1 | id, pedigree = ped)",
    "animal(1 | id) with data = hs_data(..., pedigree = ped)",
    "permanent(1 | id)",
    "common_env(1 | group)",
    "maternal_genetic(1 | dam)",
    "(1 | group)",
    "animal(rr(covariate, order = 2) | id, pedigree = ped)",
    "maternal_env(1 | dam)",
    "paternal_genetic(1 | sire, pedigree = ped)",
    "paternal_env(1 | sire)",
    "group(1 | genetic_group)",
    "unknown_parent_group(1 | upg)",
    "metafounder(1 | id, pedigree = ped, group = group, Gamma = Gamma)",
    "inbreeding(1 | id)",
    "cytoplasmic(1 | maternal_line)",
    "imprinting(1 | id, pedigree = ped, parent = \"maternal\")",
    "dominance(1 | id, pedigree = ped)",
    "epistasis(1 | id, pedigree = ped)",
    "relmat(1 | id, K = K)",
    "precision(1 | id, Q = Q)",
    "genomic(1 | id, Ginv = Ginv)",
    "genomic(1 | id, markers = M)",
    "single_step(1 | id, Hinv = Hinv)",
    "single_step(1 | id, pedigree = ped, markers = M)",
    paste0(
      "single_step(1 | id) with data = hs_data(..., pedigree = ped, ",
      "genotypes = M)"
    ),
    paste0(
      "single_step(1 | id, pedigree = ped, markers = M, group = group, ",
      "Gamma = Gamma)"
    ),
    "markers(M, model = \"random\")",
    "marker_scan(M, map = marker_map)",
    "qtl_scan(position, genotype_probs = probs)",
    "cbind(trait1, trait2) ~ animal(1 | id, pedigree = ped)",
    "animal(trait | id, pedigree = ped, cov = us())",
    "animal(trait | id, pedigree = ped, cov = diag())",
    "animal(trait | id, pedigree = ped, cov = lowrank(K = 2))",
    "animal(trait | id, pedigree = ped, cov = fa(K = 2))",
    "missing = miss_control(response = \"include\")",
    "mi(x) with missing = miss_control(predictor = \"model\")"
  )
}

hs_formula_status_categories <- function() {
  c(
    rep("v0.1 animal model", 2L),
    rep("standard quantitative genetics", 12L),
    rep("inheritance and relationship kernels", 6L),
    rep("genomic and marker models", 9L),
    rep("multivariate and factor analytic", 5L),
    rep("missing-data grammar", 2L)
  )
}

hs_formula_status_phases <- function() {
  c(
    rep("Phase 1", 2L),
    rep("Phase 2", 12L),
    rep("Phase 3+", 6L),
    rep("Phase 5", 9L),
    rep("Phase 3-4", 5L),
    rep("Phase 8", 2L)
  )
}

hs_formula_status_syntax <- function() {
  c(
    rep("parsed", 7L),
    rep("reserved", 5L),
    "parsed",
    # inbreeding, cytoplasmic, imprinting, dominance, epistasis (reserved), then
    # the experimental opt-in relmat()/precision() supplied-relationship terms.
    rep("reserved", 5L),
    rep("parsed", 2L),
    rep("parsed", 6L),
    rep("reserved", 3L),
    "parsed",
    # cbind() multivariate is parsed; long-format cov=us/diag/lowrank stay
    # planned; cov=fa is reserved (Julia V4-FA engine-covered, R not fitted);
    # missing-data rows stay planned.
    "planned",
    "planned",
    "planned",
    "reserved",
    "planned",
    "planned"
  )
}

hs_formula_status_fitting <- function() {
  c(
    rep("fitted (v0.1 default)", 2L),
    "fitted (opt-in repeatability)",
    "fitted (opt-in common-environment)",
    "fitted (opt-in maternal)",
    "fitted (opt-in multi-effect)",
    "fitted (opt-in random-regression)",
    rep("not available", 5L),
    "fitted (opt-in supplied-Gamma metafounder)",
    rep("not available", 5L),
    "fitted (opt-in supplied relationship, experimental)",
    "fitted (opt-in supplied precision, experimental)",
    "fitted (opt-in genomic, covered validation-scale)",
    "fitted (opt-in genomic / SNP-BLUP)",
    "fitted (opt-in single-step)",
    "fitted (opt-in single-step construction)",
    "fitted (opt-in single-step bundle construction)",
    "fitted (opt-in supplied-Gamma H^Gamma)",
    rep("not available", 3L),
    "fitted (default route, covered multivariate)",
    rep("not available", 6L)
  )
}

hs_formula_status_behavior <- function() {
  inert_marker_text <- paste(
    "Exported as an inert marker; hsquared() errors as planned, not",
    "implemented."
  )
  c(
    paste(
      "Parsed and fitted by the default v0.1 path (Gaussian animal model,",
      "REML through the HSquared.jl engine)."
    ),
    paste(
      "Fitted by the default v0.1 path when data is an hs_data() bundle with a",
      "pedigree component (Gaussian animal model, REML)."
    ),
    paste(
      "Permanent-environment effect of the opt-in, experimental repeatability",
      "model; requires an animal() term, repeated records, and",
      "engine = \"julia\", target = \"repeatability\". THIS IS THE STANDARD",
      "REPEATED-MEASURES ANIMAL MODEL: with more than one record per individual",
      "and no permanent() term, sigma^2_a ABSORBS the permanent-environment",
      "variance and both it and h2 are inflated (the fit warns). Add",
      "scale_method = \"auto\" to engine_control for any real pedigree -- the",
      "default dense estimator is capped at nobs^2 + nanimals^2 <=",
      "max_dense_cells. Returns variance-component and h2 standard errors on",
      "the sparse route; experimental, REML only, not coverage-calibrated.",
      "cbind(...) + permanent() is a named reject (hsquared#237); see",
      "docs/design/57-mv-pe-cbind-permanent-237.md -- no multi-trait PE fitter yet."
    ),
    paste(
      "Common-environment effect of the opt-in two-effect model (additive",
      "animal genetic + IID common environment, A2 = I); requires an animal()",
      "term and engine = \"julia\", target = \"two_effect\". COVERED at validation",
      "scale (opt-in; NOT the default fit path): mirrors the twin",
      "V3-TWOEFFECT-REML covered gate - a pre-declared 48-seed bias/MCSE recovery",
      "gate PASSED + a blupf90+ same-estimand REML comparator agrees ~1e-5",
      "(sommer cross-check ~2e-5); the h2/c2 interval is asymptotic delta-method",
      "and NOT coverage-calibrated. The MATERNAL genetic leg (maternal_genetic(),",
      "A2 = pedigree, same target) uses the same estimator with exact live parity",
      "but STAYS experimental."
    ),
    paste(
      "Maternal genetic effect. TWO opt-in paths. (1) The EXPERIMENTAL two-effect",
      "model (A2 = pedigree A via the dam, INDEPENDENT of the direct effect);",
      "requires an animal() term and engine = \"julia\", target = \"two_effect\".",
      "(2) The CORRELATED direct-maternal 2x2 G model (target = \"direct_maternal\"),",
      "which estimates the direct-maternal genetic covariance sigma_dm - COVERED at",
      "validation scale (opt-in; NOT the default fit path): engine V4-DIRECT-MATERNAL",
      "covered via a pre-declared 48-seed bias/MCSE gate PASSED (48/48, all four",
      "|bias|<=2*MCSE) + a sommer 4.4.5 covm() same-estimand REML comparator AGREE.",
      "Willham fence: heritability() returns the LABELLED TRIPLE (direct h2_d,",
      "maternal m2, Willham total h2_T, r_am), never a bare scalar; direct h2 is NOT",
      "the total heritability; a negative r_am is real and expected. SCOPE: single",
      "relationship A, dense/validation-scale n<=1000, |r_am|->1 rides on converged.",
      "The two-effect maternal leg (path 1) STAYS experimental."
    ),
    paste(
      "Bare (1 | group) i.i.d. random intercept of the opt-in multi-effect model;",
      "combines with an animal() term (and, for K independent effects, additional",
      "(1 | group) terms) and requires engine = \"julia\", target = \"multi_effect\".",
      "COVERED at validation scale (opt-in; the arbitrary-N INDEPENDENT",
      "generalization of the two-effect model): mirrors the twin V3-NEFFECT-REML",
      "covered gate (48-seed bias/MCSE gate PASSED + a sommer same-estimand REML",
      "comparator) with exact live R-Julia parity. Random slopes (x | group) and",
      "correlated (x || group) terms remain rejected. The animal-block ratio is",
      "narrow-sense h2; other blocks are variance-explained proportions (not",
      "heritabilities). Point estimates plus per-component ratio intervals",
      "(heritability_interval() resolves to the animal ratio; other blocks in",
      "variance_ratio_intervals); intervals are asymptotic delta-method and NOT",
      "coverage-calibrated. INDEPENDENT effects only (NOT correlated /",
      "random-regression / non-Gaussian)."
    ),
    paste(
      "Opt-in random-regression (reaction-norm) model: rr(covariate, order = k) on",
      "the animal() left-hand side fits a k-coefficient normalized-Legendre polynomial",
      "of a within-individual covariate; requires repeated records and engine =",
      "\"julia\", target = \"random_regression\". COVERED at k=2 (linear reaction norm,",
      "intercept + ONE slope, K_g 2x2, Gaussian, homogeneous residual, D=I2",
      "normalized Legendre): random-regression REML validation row covered via pre-declared 48-seed",
      "bias/MCSE gate PASSED + sommer 4.4.5 leg() same-estimand REML comparator",
      "AGREE (<=1.9e-5); live R<->engine parity EXACT (<=1.03e-5); h2(t) <=4.24e-6.",
      "rr_heritability() returns h2(t) as a CURVE, never a scalar; heritability()",
      "is NOT the RR accessor and errors on this result, naming rr_heritability().",
      "PE-overstatement caveat when no permanent-environment term. POINT-ESTIMATE",
      "only (no interval/CI on the K_g curve). k>=3 experimental; (x|g) raw slopes",
      "rejected; PE term and heterogeneous residual deferred."
    ),
    rep(inert_marker_text, 5L),
    paste(
      "Primary metafounder effect of the opt-in, experimental supplied-variance",
      "model; requires `group`, supplied `Gamma`, engine = \"julia\",",
      "target = \"metafounder\", and supplied variance components.",
      "`Gamma` and the variance components are supplied, not estimated."
    ),
    paste(
      "Exported as an inert marker; hsquared() errors as planned, not",
      "implemented. Inbreeding coefficients F are already computed internally",
      "for Ainv construction in the engine; this reserved term is the future",
      "user-facing F-as-effect surface, not yet fittable."
    ),
    rep(inert_marker_text, 4L),
    paste(
      "Primary supplied-relationship effect of the opt-in, EXPERIMENTAL model:",
      "the user supplies a dense symmetric positive-definite relationship",
      "matrix `K` (ID-keyed via row/column names); the parser marshals the",
      "inverse Kinv = solve(K) and fits through the SAME supplied-relationship-",
      "inverse REML path as genomic()'s `Ginv` (engine = \"julia\", target =",
      "\"relmat\"). NOT covered, NOT the default; `K` is supplied provenance,",
      "not estimated."
    ),
    paste(
      "Primary supplied-precision effect of the opt-in, EXPERIMENTAL model: the",
      "user supplies the precision (inverse) `Q` directly (ID-keyed; `Kinv =`",
      "also accepted); fits through the SAME supplied-relationship-inverse REML",
      "path as genomic()'s `Ginv` (engine = \"julia\", target = \"precision\").",
      "NOT covered, NOT the default; `Q` is supplied provenance, not estimated."
    ),
    paste(
      "Primary genomic effect for the narrow opt-in Gaussian REML route; a",
      "supplied `Ginv` is used without alteration and its construction method,",
      "allele frequencies, ridge, and denominator remain unknown. NOT on the",
      "default `engine = \"fit\"` path: ordinary-call / default-route genomic",
      "activation remains held (design-44 G5; BOUNDARY_HOLDOUT_FAIL). The",
      "explicit model requires engine = \"julia\", target = \"genomic\". The",
      "coefficient-scale result is labelled `genomic_variance_ratio`; interval",
      "and SE accessors are unavailable. Covered at validation scale (0.7),",
      "experimental; public_covered_count remains 7."
    ),
    paste(
      "Primary genomic effect for the narrow opt-in Gaussian REML marker route:",
      "sample allele frequencies, unweighted VanRaden method 1, ridge 0.01.",
      "NOT on the default `engine = \"fit\"` path: ordinary-call genomic",
      "activation remains held. Use target = \"genomic\"; target = \"snp_blup\"",
      "remains opt-in. The coefficient-scale result is labelled",
      "`genomic_variance_ratio`; interval and SE accessors are unavailable.",
      "Covered at validation scale (0.7), experimental; public_covered_count",
      "remains 7."
    ),
    paste(
      "Primary single-step effect of the opt-in, experimental supplied-inverse",
      "route; requires a user-supplied ID-keyed `Hinv` and",
      "engine = \"julia\", target = \"single_step\". This route consumes the",
      "relationship precision as supplied: it does not construct H^-1 from",
      "pedigree and markers. To construct H^-1, use the separate",
      "`single_step(1 | id, pedigree = ped, markers = M)` route and",
      "target = \"single_step_construct\"."
    ),
    paste(
      "Primary single-step effect of the opt-in, experimental construction path;",
      "requires `pedigree` + `markers` or an hs_data() bundle carrying both, and",
      "engine = \"julia\", target = \"single_step_construct\". The engine builds",
      "H^-1 from pedigree and genotyped-subset markers at validation scale.",
      "HSquared.jl #295 records AGHmatrix Hmatrix construction AGREE and an",
      "n=240 recovery PASS; that is engine evidence, not an R-public covered",
      "claim. This is separate from the supplied-`Hinv` route, is not a",
      "default-route promotion, and public_covered_count stays 7."
    ),
    paste(
      "Bundle shorthand for the opt-in, experimental single-step construction",
      "path; when data is an hs_data() object with pedigree and genotypes,",
      "`single_step(1 | id)` resolves both from the bundle. Explicit `pedigree`",
      "or `markers` arguments override the bundle. Requires engine = \"julia\",",
      "target = \"single_step_construct\"; a plain data.frame still needs an",
      "explicit pedigree and marker matrix. Engine evidence is HSquared.jl #295;",
      "R is not a default-route / ordinary-path promotion and",
      "public_covered_count stays 7."
    ),
    paste(
      "Primary single-step effect of the opt-in, experimental supplied-Gamma",
      "H^Gamma path; validates ID-keyed metafounder `group`, supplied",
      "symmetric positive-semidefinite `Gamma`, marker ordering, and",
      "`genotyped_rows`, then fits with engine = \"julia\", target =",
      "\"metafounder_single_step\". `Gamma` is supplied, not estimated."
    ),
    rep(inert_marker_text, 3L),
    paste(
      "Covered at validation scale (2026-09-02 maintainer sign-off) multivariate Gaussian animal model;",
      "a `cbind()` response with an `animal()` term routes to the multivariate",
      "fitter on the DEFAULT path (no engine/target argument needed; the explicit",
      "engine = \"julia\", target = \"multivariate\" spelling still works).",
      "Covered numeric claim is scoped to k = 2 unstructured G0/R0;",
      "k >= 3 traits stay parseable-and-fittable-but-experimental;",
      "genetic_structure = \"diagonal\" stays experimental at 0.6.",
      "Missing trait cells are allowed as `NA`. Under",
      "`family = binomial()`, `cbind(successes, failures)` is instead a",
      "binomial-counts GLMM via target = \"nongaussian\" (scalar or varying",
      "positive row totals admitted), not a multivariate Gaussian."
    ),
    paste(
      "Roadmap syntax for long-format unstructured covariance; the current",
      "parser rejects trait and `cov` arguments and points users to the",
      "`cbind()` multivariate path, which fits on the default path."
    ),
    paste(
      "Roadmap syntax for long-format diagonal covariance; the current",
      "parser rejects trait and `cov` arguments. Rotation-free diagonal G0 is",
      "already reachable as engine_control$genetic_structure = \"diagonal\"",
      "on the cbind() multivariate path."
    ),
    paste(
      "Roadmap syntax for long-format low-rank covariance (Lambda Lambda');",
      "the current parser rejects trait and `cov` arguments. Not R-activated."
    ),
    paste(
      "Reserved stub for long-format factor-analytic G (Lambda Lambda' + Psi).",
      "Julia V4-FA is engine-covered (HSquared.jl 60895208 / #300). The R",
      "parser rejects trait and `cov` arguments and does not fit this form.",
      "Not an R-public FA claim. Use cbind() unstructured or",
      "genetic_structure = \"diagonal\". public_covered_count stays 7."
    ),
    paste(
      "Ratified planned missing-response control. Future behavior will keep",
      "rows with missing responses in the model frame and mask their direct",
      "likelihood contribution. Current behavior remains complete-case or",
      "target-specific handling only; no miss_control() function is exported."
    ),
    paste(
      "Ratified planned missing-predictor grammar. Future behavior will allow",
      "one bare missing covariate declared as mi(x), paired with impute =",
      "list(x = x ~ ...), and integrated by a model-based Laplace path. Current",
      "parser behavior remains unsupported; no mi() marker is exported."
    )
  )
}
