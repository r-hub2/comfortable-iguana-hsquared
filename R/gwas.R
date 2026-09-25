#' Post-fit relatedness-corrected marker scan (GWAS)
#'
#' `r lifecycle::badge("experimental")`
#'
#' `gwas()` runs a dense, supplied-variance, relatedness-corrected mixed-model
#' (GLS) Wald marker scan on a fitted Gaussian animal model, reusing the fit's
#' estimated variance components `(sigma^2_a, sigma^2_e)` and pedigree relationship. It is an
#' experimental, validation-scale screen that surfaces the Julia-owned
#' `HSquared.mixed_model_marker_scan()`.
#'
#' **The default `p_value`/`bonferroni_p`/`bh_qvalue` columns are NOT genome-wide
#' calibrated.** They are marker-by-marker Wald (nominal) p-values plus
#' deterministic Bonferroni and Benjamini-Hochberg adjustments over the *supplied*
#' marker set only. Do not report genome-wide significance from those columns.
#'
#' **Genome-wide calibration** is available via `genome_wide = TRUE` (which requires
#' `method = "single"`, the validated fixed-effect scope). It adds a `genome_wide_p`
#' column: the exact per-dataset add-one permutation p-value (the marker is declared
#' genome-wide significant when `genome_wide_p <= alpha`, `alpha = 0.05`). For each
#' analysis the permutation null is rebuilt from the analysed phenotype (`y`
#' permuted conditional on `X`, re-scanned `n_permutations` times), and the
#' genome-wide p is `(1 + #{null max >= observed})/(n_permutations + 1)`. This is the
#' HSquared.jl `genome_wide_marker_scan` engine call, whose family-wise type-I
#' control is validated at validation and production scale (the per-dataset add-one
#' REBUILD gate; the anti-conservative `(1-alpha)` quantile rule is NOT used). The
#' result carries a `calibration` attribute (method `permutation_addone`,
#' `empirical_type1 = NA` because the per-dataset rule's validity is by construction
#' and externally validated, named in `validation_reference`). SCOPE: fixed-effect /
#' intercept-only; the relatedness-corrected mixed-model genome-wide null is NOT yet
#' calibrated, so `genome_wide = TRUE` is rejected for `method = "mixed"`/`"loco"`.
#'
#' A leave-one-group-out (LOCO) scan is available with `method = "loco"` and a
#' `marker_groups` argument: when a marker is tested, the genomic relationship
#' correction is built from the markers **not** in that marker's group (e.g. its
#' chromosome), so the marker's own signal does not leak into the background
#' relationship. The LOCO relationship is **genomic** (VanRaden) while the reused
#' variance components are **pedigree**-estimated -- a scale mismatch that keeps
#' this validation-scale (see `method`).
#'
#' @param object A fitted Gaussian animal model (`hsquared_fit` from the default
#'   pedigree path); its variance components and pedigree relationship are reused
#'   so the scan is conditioned on the same covariance the model was fit under.
#' @param markers A numeric matrix of marker dosages with one row per animal in
#'   the fit's pedigree (in pedigree order) and one column per marker.
#' @param marker_ids Optional marker names; defaults to the `markers` column
#'   names, then to sequential ids.
#' @param method `"mixed"` (default) for the relatedness-corrected mixed-model
#'   (GLS) scan with one whole-pedigree relationship correction across all
#'   markers; `"single"` for the relatedness-**un**corrected single-marker (OLS)
#'   scan (a naive screen useful mainly as a contrast -- it is more inflated by
#'   relatedness than the mixed scan); or `"loco"` for a leave-one-group-out scan
#'   with a per-group **genomic** relationship correction (requires
#'   `marker_groups`). The LOCO scan reuses the pedigree fit's variance components
#'   while correcting with a genomic relationship (a scale mismatch), so it is
#'   validation-scale and uncalibrated like the others.
#' @param marker_groups Required for `method = "loco"` (and only then): a vector
#'   with one group label per marker column (for example a chromosome label).
#'   Markers in a group are tested with a genomic relationship built from all
#'   **other** groups. Needs at least two distinct, non-missing labels.
#' @param genome_wide Logical; if `TRUE` (requires `method = "single"`), add a
#'   genome-wide-calibrated `genome_wide_p` column via the exact per-dataset add-one
#'   permutation rule (see Details). Defaults to `FALSE` (nominal p-values only).
#' @param n_permutations Number of permutations for the genome-wide null when
#'   `genome_wide = TRUE` (default 1000). The add-one floor is `1/(n_permutations+1)`.
#' @param seed Integer RNG seed for the genome-wide permutation null (default 1), so
#'   the `genome_wide_p` column is reproducible for a given Julia version.
#' @param ... Unused.
#'
#' @return An `hs_gwas` data frame with one row per marker: `marker`, `effect`,
#'   `se`, `z`, `chisq`, `p_value`, `bonferroni_p`, `bh_qvalue`, `lod` (and, when
#'   `genome_wide = TRUE`, `genome_wide_p`), carrying a `scan_method` attribute (and,
#'   when `genome_wide = TRUE`, a `calibration` attribute). Its `print()` method
#'   restates the uncalibrated-significance caveat for the nominal columns (and, for
#'   `method = "single"`, the absence of any relatedness correction; for
#'   `method = "loco"`, the genomic-vs-pedigree scale mismatch).
#' @export
gwas <- function(
  object,
  markers,
  marker_ids = NULL,
  method = c("mixed", "single", "loco"),
  marker_groups = NULL,
  genome_wide = FALSE,
  n_permutations = 1000L,
  seed = 1L,
  ...
) {
  UseMethod("gwas")
}

#' @export
gwas.default <- function(
  object,
  markers,
  marker_ids = NULL,
  method = c("mixed", "single", "loco"),
  marker_groups = NULL,
  ...
) {
  stop(
    "`gwas()` requires a fitted Gaussian animal model (`hsquared_fit`).",
    call. = FALSE
  )
}

#' @export
gwas.hsquared_fit <- function(
  object,
  markers,
  marker_ids = NULL,
  method = c("mixed", "single", "loco"),
  marker_groups = NULL,
  genome_wide = FALSE,
  n_permutations = 1000L,
  seed = 1L,
  ...
) {
  method <- match.arg(method)
  genome_wide <- isTRUE(genome_wide)
  if (genome_wide) {
    # The genome-wide permutation calibration is validated for the FIXED-effect
    # scan only (the exact per-dataset add-one rule; HSquared.jl REBUILD gate). The
    # relatedness-corrected mixed-model / LOCO genome-wide null is a different,
    # not-yet-validated calibration, so `genome_wide = TRUE` requires `method = "single"`.
    if (!identical(method, "single")) {
      stop(
        "`genome_wide = TRUE` is supported only for `method = \"single\"` (the ",
        "validated fixed-effect permutation calibration). The relatedness-corrected ",
        "mixed-model / LOCO genome-wide null is not yet calibrated.",
        call. = FALSE
      )
    }
    n_permutations <- as.integer(n_permutations)
    if (length(n_permutations) != 1L || is.na(n_permutations) || n_permutations < 1L) {
      stop("`n_permutations` must be a positive whole number.", call. = FALSE)
    }
    seed <- as.integer(seed)
    if (length(seed) != 1L || is.na(seed)) {
      stop("`seed` must be a single integer.", call. = FALSE)
    }
  }
  hs_validate_gwas_fit(object)
  payload <- object$payload
  vc <- object$result$variance_components
  sigma_a2 <- vc$estimate[vc$component == "animal"][[1L]]
  sigma_e2 <- vc$estimate[vc$component == "residual"][[1L]]

  markers <- hs_validate_gwas_markers(markers, payload)
  marker_ids <- hs_gwas_marker_ids(marker_ids, markers)
  marker_groups <- hs_gwas_marker_groups(marker_groups, markers, method)

  project <- hs_default_julia_project()
  if (!hs_julia_bridge_available(project)) {
    stop(
      "The experimental Julia bridge requires Julia, the `JuliaCall` R ",
      "package, and a local `HSquared.jl` project.",
      call. = FALSE
    )
  }
  markers_rec <- as.matrix(payload$Z %*% markers)

  hs_julia_setup(project)
  JuliaCall::julia_assign("hsq_y", payload$y)
  JuliaCall::julia_assign("hsq_X", payload$X)
  JuliaCall::julia_assign("hsq_markers", markers_rec)
  JuliaCall::julia_assign("hsq_sigma_e2", as.numeric(sigma_e2))
  JuliaCall::julia_assign("hsq_marker_ids", as.character(marker_ids))
  if (identical(method, "mixed")) {
    hs_julia_assign_sparse_csc("hsq_Z", payload$Z)
    JuliaCall::julia_assign("hsq_id", payload$pedigree$id)
    JuliaCall::julia_assign(
      "hsq_sire",
      hs_parent_for_julia(payload$pedigree$sire)
    )
    JuliaCall::julia_assign(
      "hsq_dam",
      hs_parent_for_julia(payload$pedigree$dam)
    )
    JuliaCall::julia_assign("hsq_sigma_a2", as.numeric(sigma_a2))
    scan_cmd <- paste(
      "hsq_ped = HSquared.normalize_pedigree(hsq_id, hsq_sire, hsq_dam);",
      "hsq_Ainv = HSquared.pedigree_inverse(hsq_ped);",
      "hsq_scan = HSquared.mixed_model_marker_scan(",
      "hsq_y, hsq_X, hsq_Z, hsq_Ainv, hsq_markers, hsq_sigma_a2, hsq_sigma_e2;",
      "marker_ids = hsq_marker_ids);"
    )
  } else if (identical(method, "loco")) {
    # leave-one-group-out: precisions from ANIMAL-level markers (n_animals x
    # n_animals, the Ainv slot); the scan tests the RECORD-level markers. Reuses
    # the pedigree fit's variance components (genomic-vs-pedigree scale mismatch
    # is stated in print()/docs).
    hs_julia_assign_sparse_csc("hsq_Z", payload$Z)
    JuliaCall::julia_assign("hsq_markers_animal", markers)
    JuliaCall::julia_assign("hsq_groups", as.character(marker_groups))
    JuliaCall::julia_assign("hsq_sigma_a2", as.numeric(sigma_a2))
    scan_cmd <- paste(
      "hsq_prec = HSquared.loco_relationship_precisions(",
      "hsq_markers_animal, hsq_groups);",
      "hsq_scan = HSquared.loco_mixed_model_marker_scan(",
      "hsq_y, hsq_X, hsq_Z, hsq_prec, hsq_groups, hsq_markers,",
      "hsq_sigma_a2, hsq_sigma_e2; marker_ids = hsq_marker_ids);"
    )
  } else if (genome_wide) {
    # relatedness-UNcorrected fixed-effect scan WITH genome-wide permutation
    # calibration: the exact per-dataset add-one rule (`genome_wide_marker_scan`),
    # the HSquared.jl REBUILD-gate-validated procedure. The RNG is pinned
    # (MersenneTwister(seed)) so the permutation null is reproducible.
    JuliaCall::julia_assign("hsq_nperm", as.integer(n_permutations))
    JuliaCall::julia_assign("hsq_alpha", 0.05)
    JuliaCall::julia_assign("hsq_seed", as.integer(seed))
    scan_cmd <- paste(
      "import Random;",
      "hsq_scan = HSquared.genome_wide_marker_scan(",
      "hsq_y, hsq_X, hsq_markers;",
      "n_permutations = hsq_nperm, alpha = hsq_alpha, sigma_e2 = hsq_sigma_e2,",
      "marker_ids = hsq_marker_ids, rng = Random.MersenneTwister(hsq_seed));"
    )
  } else {
    # relatedness-UNcorrected single-marker (OLS) scan: no Z / Ainv / sigma_a2
    scan_cmd <- paste(
      "hsq_scan = HSquared.single_marker_scan(",
      "hsq_y, hsq_X, hsq_markers;",
      "sigma_e2 = hsq_sigma_e2, marker_ids = hsq_marker_ids);"
    )
  }
  gw_dict <- if (genome_wide) {
    paste(
      ",\"genome_wide_p_values\" => collect(Float64, hsq_scan.genome_wide_p_values),",
      "\"genome_wide_threshold\" => hsq_scan.genome_wide_threshold,",
      "\"genome_wide_p_min\" => hsq_scan.genome_wide_p_min,",
      "\"n_permutations\" => hsq_scan.n_permutations,",
      "\"alpha\" => hsq_scan.alpha,",
      "\"calibration_method\" => string(hsq_scan.calibration.method),",
      "\"marker_panel_mode\" => string(hsq_scan.calibration.marker_panel_mode)"
    )
  } else {
    ""
  }
  JuliaCall::julia_command(paste(
    scan_cmd,
    "hsq_gwas_raw = Dict{String,Any}(",
    "\"marker_ids\" => string.(collect(hsq_scan.marker_ids)),",
    "\"effects\" => collect(Float64, hsq_scan.effects),",
    "\"standard_errors\" => collect(Float64, hsq_scan.standard_errors),",
    "\"z_scores\" => collect(Float64, hsq_scan.z_scores),",
    "\"chisq\" => collect(Float64, hsq_scan.chisq),",
    "\"p_values\" => collect(Float64, hsq_scan.p_values),",
    "\"bonferroni\" => collect(Float64, hsq_scan.bonferroni_p_values),",
    "\"bh\" => collect(Float64, hsq_scan.bh_q_values),",
    "\"lod\" => collect(Float64, hsq_scan.lod_scores)",
    gw_dict,
    ");"
  ))
  raw <- JuliaCall::julia_eval("hsq_gwas_raw")
  if (genome_wide) {
    # Build the calibration metadata for the per-dataset permutation rule. It has
    # no per-call empirical type-I (validity is by construction + externally
    # validated), so `empirical_type1 = NA` and `validation_reference` names the
    # HSquared.jl type-I-control gate. The chi-square genome-wide threshold is
    # carried on the LOD scale.
    raw$calibration <- list(
      calibration_method = raw$calibration_method %||% "permutation_addone",
      threshold_scale = "lod",
      threshold = (raw$genome_wide_threshold %||% NA_real_) / (2 * log(10)),
      alpha = raw$alpha %||% 0.05,
      empirical_type1 = NA_real_,
      marker_panel_mode = "real_panel",
      scan_method = method,
      n_replicates = raw$n_permutations %||% n_permutations,
      seed = seed,
      engine = "HSquared.genome_wide_marker_scan",
      package_version = as.character(utils::packageVersion("hsquared")),
      validation_reference = paste0(
        "HSquared.jl production REBUILD gate ",
        "(sim/phase5_qtl_rebuild_production_gate.jl): per-dataset add-one ",
        "permutation type-I 0.0504/0.0542 at alpha=0.05 (fixed-effect, ",
        "intercept-only); the (1-alpha) quantile rule is anti-conservative (#202)."
      )
    )
  }
  hs_normalize_gwas_result(raw, method = method)
}

# A gwas() fit must be the default univariate Gaussian pedigree animal model:
# it reuses the fit's pedigree relationship and additive/residual variance
# components for a relatedness-corrected scan.
hs_validate_gwas_fit <- function(object) {
  family <- object$spec$family$family %||% "gaussian"
  if (!identical(family, "gaussian")) {
    stop(
      "`gwas()` reuses an additive/residual variance decomposition, so it needs ",
      "a Gaussian animal-model fit; `",
      family,
      "` fits are not supported.",
      call. = FALSE
    )
  }
  payload <- object$payload
  if (
    is.null(payload$pedigree) ||
      is.null(payload$Z) ||
      is.null(payload$y) ||
      is.null(payload$X)
  ) {
    stop(
      "`gwas()` requires a fitted pedigree animal model (the default ",
      "`hsquared(y ~ ... + animal(1 | id, pedigree = ped))` path).",
      call. = FALSE
    )
  }
  vc <- object$result$variance_components
  if (
    is.null(vc) ||
      !is.data.frame(vc) ||
      !all(c("animal", "residual") %in% vc$component)
  ) {
    stop(
      "`gwas()` needs a fit reporting `animal` and `residual` variance ",
      "components; multivariate, genomic, single-step, SNP-BLUP, and ",
      "second-effect fits are not supported.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

hs_validate_gwas_markers <- function(markers, payload) {
  if (is.null(markers)) {
    stop("`markers` is required.", call. = FALSE)
  }
  markers <- as.matrix(markers)
  if (!is.numeric(markers)) {
    stop("`markers` must be a numeric matrix.", call. = FALSE)
  }
  n_animals <- length(payload$pedigree$id)
  if (nrow(markers) != n_animals) {
    stop(
      "`markers` must have one row per animal in the fit's pedigree (",
      n_animals,
      " rows, in pedigree order); got ",
      nrow(markers),
      ".",
      call. = FALSE
    )
  }
  if (ncol(markers) < 1L) {
    stop("`markers` must have at least one marker column.", call. = FALSE)
  }
  if (any(!is.finite(markers))) {
    stop("`markers` must contain only finite dosages.", call. = FALSE)
  }
  storage.mode(markers) <- "double"
  markers
}

hs_gwas_marker_ids <- function(marker_ids, markers) {
  if (!is.null(marker_ids)) {
    if (length(marker_ids) != ncol(markers)) {
      stop(
        "`marker_ids` must have one entry per marker column (",
        ncol(markers),
        ").",
        call. = FALSE
      )
    }
    return(as.character(marker_ids))
  }
  labels <- colnames(markers)
  if (is.null(labels) || length(labels) != ncol(markers)) {
    labels <- as.character(seq_len(ncol(markers)))
  }
  labels
}

# LOCO needs one group label per marker column; the genomic relationship for a
# marker is built from all OTHER groups, so at least two distinct groups are
# required. marker_groups is meaningful only for method = "loco".
hs_gwas_marker_groups <- function(marker_groups, markers, method) {
  if (!identical(method, "loco")) {
    if (!is.null(marker_groups)) {
      stop(
        "`marker_groups` is only used when `method = \"loco\"`.",
        call. = FALSE
      )
    }
    return(NULL)
  }
  if (is.null(marker_groups)) {
    stop(
      "`method = \"loco\"` requires `marker_groups` (one group label per ",
      "marker column, e.g. a chromosome).",
      call. = FALSE
    )
  }
  if (length(marker_groups) != ncol(markers)) {
    stop(
      "`marker_groups` must have one entry per marker column (",
      ncol(markers),
      "); got ",
      length(marker_groups),
      ".",
      call. = FALSE
    )
  }
  if (any(is.na(marker_groups))) {
    stop("`marker_groups` must not contain missing labels.", call. = FALSE)
  }
  groups <- as.character(marker_groups)
  if (any(!nzchar(groups))) {
    stop("`marker_groups` must not contain empty labels.", call. = FALSE)
  }
  if (length(unique(groups)) < 2L) {
    stop(
      "LOCO needs at least two distinct marker groups (a marker's relationship ",
      "correction is built from the other groups).",
      call. = FALSE
    )
  }
  groups
}

hs_normalize_gwas_result <- function(raw, method = "mixed") {
  out <- data.frame(
    marker = as.character(raw$marker_ids),
    effect = as.numeric(raw$effects),
    se = as.numeric(raw$standard_errors),
    z = as.numeric(raw$z_scores),
    chisq = as.numeric(raw$chisq),
    p_value = as.numeric(raw$p_values),
    bonferroni_p = as.numeric(raw$bonferroni),
    bh_qvalue = as.numeric(raw$bh),
    lod = as.numeric(raw$lod),
    stringsAsFactors = FALSE
  )
  if (!is.null(raw$genome_wide_p_values)) {
    # genome-wide-calibrated add-one p (per marker), from the per-dataset
    # permutation null; the marker order matches the scan rows.
    out$genome_wide_p <- as.numeric(raw$genome_wide_p_values)
  }
  class(out) <- c("hs_gwas", "data.frame")
  attr(out, "scan_method") <- method
  calibration <- hs_validate_gwas_calibration_metadata(
    raw$calibration %||% NULL,
    scan_method = method
  )
  if (!is.null(calibration)) {
    attr(out, "calibration") <- calibration
  }
  out
}

hs_gwas_calibration_required_fields <- function() {
  c(
    "calibration_method",
    "threshold_scale",
    "threshold",
    "alpha",
    "empirical_type1",
    "marker_panel_mode",
    "scan_method",
    "n_replicates",
    "seed",
    "engine",
    "package_version"
  )
}

hs_validate_gwas_calibration_metadata <- function(
  calibration,
  scan_method = "mixed"
) {
  if (is.null(calibration)) {
    return(NULL)
  }
  if (!is.list(calibration)) {
    stop("GWAS calibration metadata must be a named list.", call. = FALSE)
  }
  missing <- setdiff(
    hs_gwas_calibration_required_fields(),
    names(calibration)
  )
  if (length(missing) > 0L) {
    stop(
      "GWAS calibration metadata is incomplete; missing field",
      if (length(missing) == 1L) " " else "s ",
      "`",
      paste(missing, collapse = "`, `"),
      "`.",
      call. = FALSE
    )
  }

  method <- hs_scalar_character(calibration$calibration_method)
  if (is.na(method) || !nzchar(method) || identical(method, "none")) {
    stop(
      "GWAS calibration metadata must name an activated calibration method ",
      "other than `none`.",
      call. = FALSE
    )
  }

  scale <- hs_scalar_character(calibration$threshold_scale)
  if (!scale %in% c("p_value", "lod")) {
    stop(
      "GWAS calibration `threshold_scale` must be `p_value` or `lod`.",
      call. = FALSE
    )
  }

  threshold <- hs_scalar_numeric(calibration$threshold)
  if (!is.finite(threshold) || threshold <= 0) {
    stop(
      "GWAS calibration `threshold` must be positive and finite.",
      call. = FALSE
    )
  }
  if (identical(scale, "p_value") && threshold > 1) {
    stop(
      "GWAS calibration p-value thresholds must be between 0 and 1.",
      call. = FALSE
    )
  }

  alpha <- hs_scalar_numeric(calibration$alpha)
  empirical_type1 <- hs_scalar_numeric(calibration$empirical_type1)
  if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("GWAS calibration `alpha` must be between 0 and 1.", call. = FALSE)
  }
  # The per-dataset permutation rule (`permutation_addone`) has no per-call empirical
  # type-I: its validity is by construction (the add-one permutation p is a valid exact
  # test) and is established externally by a validation gate, named in
  # `validation_reference`. For that method `empirical_type1` may be NA; for any other
  # (e.g. fixed-panel-simulation) calibration it must be a measured value in [0, 1].
  if (identical(method, "permutation_addone")) {
    if (!is.na(empirical_type1) &&
      (empirical_type1 < 0 || empirical_type1 > 1)) {
      stop(
        "GWAS calibration `empirical_type1` must be NA or between 0 and 1 for ",
        "the `permutation_addone` method.",
        call. = FALSE
      )
    }
    validation_reference <- hs_scalar_character(calibration$validation_reference)
    if (is.na(validation_reference) || !nzchar(validation_reference)) {
      stop(
        "GWAS calibration method `permutation_addone` requires a non-empty ",
        "`validation_reference` (the external type-I-control evidence), because ",
        "the per-dataset permutation rule has no per-call empirical type-I.",
        call. = FALSE
      )
    }
  } else if (
    !is.finite(empirical_type1) ||
      empirical_type1 < 0 ||
      empirical_type1 > 1
  ) {
    stop(
      "GWAS calibration `empirical_type1` must be between 0 and 1.",
      call. = FALSE
    )
  }

  panel_mode <- hs_scalar_character(calibration$marker_panel_mode)
  if (!panel_mode %in% c("fixed", "fresh", "realistic_ld", "real_panel")) {
    stop(
      "GWAS calibration `marker_panel_mode` must be one of `fixed`, `fresh`, ",
      "`realistic_ld`, or `real_panel`.",
      call. = FALSE
    )
  }

  payload_method <- hs_scalar_character(calibration$scan_method)
  if (!identical(payload_method, scan_method)) {
    stop(
      "GWAS calibration `scan_method` must match the scan result method (`",
      scan_method,
      "`).",
      call. = FALSE
    )
  }

  n_replicates <- hs_scalar_numeric(calibration$n_replicates)
  if (
    !is.finite(n_replicates) ||
      n_replicates < 1 ||
      n_replicates != as.integer(n_replicates)
  ) {
    stop(
      "GWAS calibration `n_replicates` must be a positive whole number.",
      call. = FALSE
    )
  }

  seed <- calibration$seed
  if (length(seed) < 1L || any(is.na(seed))) {
    stop(
      "GWAS calibration `seed` must contain at least one non-missing value.",
      call. = FALSE
    )
  }

  engine <- hs_scalar_character(calibration$engine)
  package_version <- hs_scalar_character(calibration$package_version)
  if (is.na(engine) || !nzchar(engine)) {
    stop("GWAS calibration `engine` must be a non-empty string.", call. = FALSE)
  }
  if (is.na(package_version) || !nzchar(package_version)) {
    stop(
      "GWAS calibration `package_version` must be a non-empty string.",
      call. = FALSE
    )
  }

  list(
    calibration_method = method,
    threshold_scale = scale,
    threshold = threshold,
    alpha = alpha,
    empirical_type1 = empirical_type1,
    marker_panel_mode = panel_mode,
    scan_method = payload_method,
    n_replicates = as.integer(n_replicates),
    seed = seed,
    engine = engine,
    package_version = package_version,
    validation_reference = hs_scalar_character(
      calibration$validation_reference %||% NA_character_
    )
  )
}

hs_scalar_character <- function(x) {
  if (length(x) != 1L) {
    return(NA_character_)
  }
  as.character(x)
}

hs_scalar_numeric <- function(x) {
  if (length(x) != 1L) {
    return(NA_real_)
  }
  as.numeric(x)
}

#' @export
print.hs_gwas <- function(x, ...) {
  method <- attr(x, "scan_method") %||% "mixed"
  if (identical(method, "single")) {
    cat(
      "<hs_gwas> dense single-marker (OLS) scan -- relatedness-UNcorrected\n"
    )
    cat(
      "  EXPERIMENTAL: no pedigree/relatedness correction (a naive contrast to\n"
    )
    cat(
      "  the relatedness-corrected `method = \"mixed\"` scan).\n"
    )
    if ("genome_wide_p" %in% names(x)) {
      cal <- attr(x, "calibration")
      cat(
        "  The `genome_wide_p` column IS genome-wide calibrated: the exact\n"
      )
      cat(
        "  per-dataset add-one permutation p (HSquared.jl `genome_wide_marker_scan`,\n"
      )
      cat(sprintf(
        "  %s, type-I-control validated). Significant when genome_wide_p <= %.2g.\n",
        cal$calibration_method %||% "permutation_addone", cal$alpha %||% 0.05
      ))
      cat(
        "  The nominal `p_value`/`bonferroni_p`/`bh_qvalue` columns remain NOT\n"
      )
      cat("  genome-wide calibrated.\n")
    } else {
      cat(
        "  p-values are NOT genome-wide calibrated (nominal Wald + Bonferroni/BH\n"
      )
      cat(
        "  only). Pass `genome_wide = TRUE` for a calibrated `genome_wide_p`\n"
      )
      cat(
        "  column. Do not report genome-wide significance from the nominal p's.\n"
      )
    }
    out <- x
    class(out) <- "data.frame"
    attr(out, "scan_method") <- NULL
    print(utils::head(out, ...), row.names = FALSE)
    return(invisible(x))
  }
  if (identical(method, "loco")) {
    cat(
      "<hs_gwas> dense leave-one-group-out (LOCO) genomic marker scan\n"
    )
    cat(
      "  EXPERIMENTAL: each marker is corrected by a genomic relationship built\n"
    )
    cat(
      "  from the OTHER marker groups. Variance components are reused from the\n"
    )
    cat(
      "  pedigree fit while the LOCO relationship is genomic (a scale mismatch),\n"
    )
    cat(
      "  so the effect/SE scale is not a calibrated genomic-VC quantity (use the\n"
    )
    cat(
      "  relative ranking, not the magnitudes). p-values are NOT genome-wide\n"
    )
    cat(
      "  calibrated (nominal Wald + Bonferroni/BH; Julia fixed-panel smoke is\n"
    )
    cat(
      "  banked, but no R significance threshold is activated).\n"
    )
    cat("  Do not report genome-wide significance from these.\n")
    out <- x
    class(out) <- "data.frame"
    attr(out, "scan_method") <- NULL
    print(utils::head(out, ...), row.names = FALSE)
    return(invisible(x))
  }
  cat("<hs_gwas> dense supplied-variance relatedness-corrected marker scan\n")
  cat(
    "  EXPERIMENTAL: p-values are NOT genome-wide calibrated. They are nominal\n"
  )
  cat(
    "  Wald p-values + Bonferroni/BH over the supplied markers only (one whole-\n"
  )
  cat(
    "  pedigree correction; no R threshold activation, no permutation, no\n"
  )
  cat(
    "  external comparator; Julia fixed-panel smoke banked in PR #134).\n"
  )
  cat(
    "  `method = \"single\"`/`\"loco\"` give the relatedness-\n"
  )
  cat(
    "  uncorrected and leave-one-group-out variants. Do not report genome-wide\n"
  )
  cat("  significance from these.\n")
  out <- x
  class(out) <- "data.frame"
  print(utils::head(out, ...), row.names = FALSE)
  invisible(x)
}
