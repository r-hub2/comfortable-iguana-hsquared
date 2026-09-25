#' Extract variance components
#'
#' `r lifecycle::badge("experimental")`
#'
#' `variance_components()` is part of the planned v0.1 fitted-object contract.
#' It works for `hsquared_fit` objects that contain a Julia result.
#'
#' @param object A fitted model object.
#' @param ... Reserved for future arguments.
#'
#' @return Variance component results for `hsquared_fit` objects.
#' @export
variance_components <- function(object, ...) {
  UseMethod("variance_components")
}

#' @export
variance_components.default <- function(object, ...) {
  stop(
    "`variance_components()` requires an `hsquared_fit` object. The current ",
    "package only returns these from fitted `hsquared_fit` results.",
    call. = FALSE
  )
}

#' @export
variance_components.hsquared_fit <- function(object, ...) {
  hs_fit_result(object, "variance_components", "variance components")
}

#' Extract heritability estimates
#'
#' `r lifecycle::badge("experimental")`
#'
#' `heritability()` is part of the planned v0.1 fitted-object contract. It
#' works for `hsquared_fit` objects that contain a Julia result.
#'
#' Falconer fence for the opt-in two-effect model (`target = "two_effect"`):
#' the reported number is the **narrow-sense direct heritability**
#' `h2 = sigma_a2 / (sigma_a2 + sigma_2 + sigma_e2)` *within that model* (the
#' additive-genetic variance is divided by the total phenotypic variance,
#' which now includes the second component `sigma_2`). The second component's
#' variance ratio (common-environment `c2` or maternal `m2`) is **not** a
#' heritability and is returned separately by [common_env_proportion()] /
#' [maternal_proportion()].
#'
#' Willham fence for the opt-in direct-maternal correlated model
#' (`target = "direct_maternal"`): `heritability()` returns the **labelled
#' triple** -- direct h2_d, maternal m2, Willham total h2_T, and r_am -- as a
#' data frame (Willham 1963, 1972). `sigma^2_P = sigma_ad + sigma_am + sigma_dm +
#' sigma_e2 = Var(y_i)` (coefficient 1 on sigma_dm). A warning is issued
#' because h2 is denominator-dependent under maternal effects and
#' h2_T < h2_d is expected when r_am < 0. Use [direct_heritability()] or
#' [total_heritability()] for targeted accessors without the warning.
#'
#' Reaction-norm fence for the opt-in random-regression model
#' (`target = "random_regression"`): heritability is a **trajectory** `h2(t)`
#' over the covariate, not a scalar, so `heritability()` **errors** on such a
#' fit and names the implemented accessor. Use [rr_heritability()] for the
#' `h2(t)` curve.
#'
#' A non-converged fit still returns the engine number so you can inspect it,
#' but **warns**: that number is not an estimate. A near-zero value from a
#' failed fit is not evidence that heritability is zero. Use
#' [fit_diagnostics()] before reading any number. `logLik()` already refuses
#' a non-converged fit; `heritability()` keeps the value and shouts instead.
#'
#' @inheritParams variance_components
#'
#' @return Heritability results for `hsquared_fit` objects.
#' @export
heritability <- function(object, ...) {
  UseMethod("heritability")
}

#' @export
heritability.default <- function(object, ...) {
  stop(
    "`heritability()` requires an `hsquared_fit` object. The current package ",
    "only returns these from fitted `hsquared_fit` results.",
    call. = FALSE
  )
}

#' @export
heritability.hsquared_fit <- function(object, ...) {
  hs_warn_if_unusable_fit(object)
  # Willham fence for the direct-maternal correlated model: heritability() on
  # a direct_maternal fit returns the LABELLED TRIPLE -- direct h2_d, maternal
  # m2, and Willham total h2_T -- plus r_am, as a labelled data frame.
  # A warning is issued to signal that h2 is denominator-dependent under
  # maternal effects and that h2_T < h2_d is real and expected when r_am < 0.
  # Use direct_heritability() or total_heritability() for targeted accessors.
  if (hs_fit_is_direct_maternal(object)) {
    vc <- object$result$variance_components
    sigma_ad <- as.numeric(vc$estimate[vc$component == "direct"])
    sigma_am <- as.numeric(vc$estimate[vc$component == "maternal"])
    sigma_dm <- as.numeric(vc$estimate[vc$component == "covariance"])
    sigma_e2 <- as.numeric(vc$estimate[vc$component == "residual"])
    sigma_P <- sigma_ad + sigma_am + sigma_dm + sigma_e2
    h2_d <- if (sigma_P > 0) sigma_ad / sigma_P else NA_real_
    m2 <- if (sigma_P > 0) sigma_am / sigma_P else NA_real_
    h2_T <- if (sigma_P > 0) {
      (sigma_ad + 1.5 * sigma_dm + 0.5 * sigma_am) / sigma_P
    } else {
      NA_real_
    }
    r_am <- as.numeric(object$result$genetic_correlation$estimate[[1L]])
    out <- data.frame(
      component = c("h2_direct", "m2_maternal", "h2_total_willham", "r_am"),
      estimate = c(h2_d, m2, h2_T, r_am),
      stringsAsFactors = FALSE
    )
    attr(out, "interpretation") <- paste(
      "Direct-maternal correlated model (target = \"direct_maternal\"):",
      "heritability() returns the labelled triple (Willham 1963, 1972).",
      "h2_direct = sigma_ad / sigma^2_P (direct narrow-sense heritability);",
      "m2_maternal = sigma_am / sigma^2_P (maternal variance ratio, NOT a",
      "heritability); h2_total_willham = (sigma_ad + 1.5*sigma_dm +",
      "0.5*sigma_am) / sigma^2_P (Willham total heritability - PREDICTS",
      "RESPONSE TO MASS SELECTION; can be LOWER than h2_direct when r_am < 0);",
      "r_am = genetic correlation between direct and maternal effects.",
      "sigma^2_P = sigma_ad + sigma_am + sigma_dm + sigma_e2 (Willham 1972).",
      "h2 is denominator-dependent under maternal effects; compare",
      "(co)variance components rather than h2 across software.",
      "Use direct_heritability() or total_heritability() for targeted accessors."
    )
    warning(
      "`heritability()` on a `direct_maternal` fit returns the labelled ",
      "triple (h2_direct, m2_maternal, h2_total_willham, r_am; Willham 1963, ",
      "1972). h2_total_willham < h2_direct is expected when r_am < 0. ",
      "Use `direct_heritability()` or `total_heritability()` to suppress ",
      "this warning.",
      call. = FALSE
    )
    return(out)
  }
  # Reaction-norm fence: h2 is a trajectory h2(t), not a scalar, on the opt-in
  # random-regression target, so `heritability()` has no scalar answer to give.
  # Name the implemented accessor rather than falling through to the generic
  # missing-field miss, which used to read as "planned v0.1 contract" / "not
  # implemented" on a route that is covered at validation scale.
  if (hs_fit_is_random_regression(object)) {
    stop(
      "`heritability()` is not defined for the opt-in random-regression ",
      "(reaction-norm) model (`target = \"random_regression\"`), where ",
      "heritability is a TRAJECTORY h2(t) over the covariate rather than a ",
      "scalar. Use `rr_heritability()` for the h2(t) curve, and ",
      "`rr_genetic_variance()` / `rr_correlation()` / `rr_eigenfunctions()` ",
      "for the other reaction-norm trajectories.",
      call. = FALSE
    )
  }
  hs_fit_result(object, "heritability", "heritability estimates")
}

hs_fit_is_genomic <- function(object) {
  inherits(object, "hsquared_fit") &&
    identical(object$spec$target, "genomic")
}

#' Extract multivariate covariance and correlation matrices
#'
#' `r lifecycle::badge("experimental")`
#'
#' These extractors return the genetic (`G`) and residual (`R`) covariance or
#' correlation matrices from multivariate `hsquared_fit` objects (a `cbind()`
#' response, which fits on the default path). `G_matrix()` is an applied-workflow alias for
#' `genetic_covariance()`, and `R_matrix()` is an alias for
#' `residual_covariance()`. Use them after checking [fit_diagnostics()] because
#' likelihood-based summaries are intentionally blocked when a multivariate fit
#' has not converged.
#'
#' @inheritParams variance_components
#'
#' @return A numeric matrix for `hsquared_fit` objects that contain the
#'   requested multivariate result field.
#'
#' @examplesIf FALSE
#' fit_mv <- hsquared(
#'   cbind(weight, length) ~ sex + age + animal(1 | id, pedigree = ped),
#'   data = dat,
#'   family = gaussian(),
#'   REML = TRUE,
#'   control = hs_control(
#'     engine = "julia",
#'     engine_control = list(target = "multivariate")
#'   )
#' )
#'
#' fit_diagnostics(fit_mv)
#'
#' genetic_covariance(fit_mv)
#' G_matrix(fit_mv)
#' residual_covariance(fit_mv)
#' R_matrix(fit_mv)
#' genetic_correlation(fit_mv)
#' residual_correlation(fit_mv)
#' heritability(fit_mv)
#' @name multivariate_extractors
NULL

#' @rdname multivariate_extractors
#' @export
genetic_covariance <- function(object, ...) {
  UseMethod("genetic_covariance")
}

#' @export
genetic_covariance.default <- function(object, ...) {
  hs_multivariate_extractor_default("genetic_covariance")
}

#' @export
genetic_covariance.hsquared_fit <- function(object, ...) {
  hs_fit_result(object, "genetic_covariance", "genetic covariance matrix")
}

#' @rdname multivariate_extractors
#' @export
G_matrix <- function(object, ...) {
  UseMethod("G_matrix")
}

#' @export
G_matrix.default <- function(object, ...) {
  hs_multivariate_extractor_default("G_matrix")
}

#' @export
G_matrix.hsquared_fit <- function(object, ...) {
  genetic_covariance(object, ...)
}

#' @rdname multivariate_extractors
#' @export
residual_covariance <- function(object, ...) {
  UseMethod("residual_covariance")
}

#' @export
residual_covariance.default <- function(object, ...) {
  hs_multivariate_extractor_default("residual_covariance")
}

#' @export
residual_covariance.hsquared_fit <- function(object, ...) {
  hs_fit_result(object, "residual_covariance", "residual covariance matrix")
}

#' @rdname multivariate_extractors
#' @export
R_matrix <- function(object, ...) {
  UseMethod("R_matrix")
}

#' @export
R_matrix.default <- function(object, ...) {
  hs_multivariate_extractor_default("R_matrix")
}

#' @export
R_matrix.hsquared_fit <- function(object, ...) {
  residual_covariance(object, ...)
}

#' @rdname multivariate_extractors
#' @export
genetic_correlation <- function(object, ...) {
  UseMethod("genetic_correlation")
}

#' @export
genetic_correlation.default <- function(object, ...) {
  hs_multivariate_extractor_default("genetic_correlation")
}

#' @export
genetic_correlation.hsquared_fit <- function(object, ...) {
  out <- hs_fit_result(
    object,
    "genetic_correlation",
    "genetic correlation matrix"
  )
  # Warn when the direct-maternal genetic correlation is at or near the
  # boundary: |r_am| >= 0.99 signals a poorly identified G_dm (shallow
  # pedigree, too few offspring per dam) or a genuine constraint boundary.
  # A negative r_am is NOT a warning -- it is real and biologically expected.
  if (hs_fit_is_direct_maternal(object) && is.data.frame(out)) {
    r_am <- out$estimate[[1L]]
    if (!is.na(r_am) && abs(r_am) >= 0.99) {
      warning(
        "The direct-maternal genetic correlation (r_am = ",
        round(r_am, 4L),
        ") is at or near the boundary |r_am| = 1. This may indicate a poorly ",
        "identified G_dm (shallow pedigree, few offspring per dam, or collinear ",
        "Zd/Zm columns). Inspect `converged` via `fit_diagnostics()` and check ",
        "the pedigree structure before interpreting this estimate.",
        call. = FALSE
      )
    }
  }
  out
}

#' @rdname multivariate_extractors
#' @export
residual_correlation <- function(object, ...) {
  UseMethod("residual_correlation")
}

#' @export
residual_correlation.default <- function(object, ...) {
  hs_multivariate_extractor_default("residual_correlation")
}

#' @export
residual_correlation.hsquared_fit <- function(object, ...) {
  hs_fit_result(object, "residual_correlation", "residual correlation matrix")
}

hs_multivariate_extractor_default <- function(name) {
  stop(
    "`",
    name,
    "()` requires an `hsquared_fit` object from the multivariate model ",
    "(a `cbind(trait1, trait2, ...)` response with `animal(1 | id, pedigree = ",
    "ped)`, which fits on the default path) or the opt-in direct-maternal ",
    "correlated model (`target = \"direct_maternal\"`).",
    call. = FALSE
  )
}

#' Extract a supplied metafounder Gamma matrix
#'
#' `r lifecycle::badge("experimental")`
#'
#' `gamma_matrix()` returns the supplied metafounder `Gamma` matrix carried by
#' an experimental metafounder or `H^Gamma` single-step `hsquared_fit` object.
#' It is provenance for the fitted relationship, not an estimated parameter.
#'
#' @inheritParams variance_components
#'
#' @return A numeric matrix with metafounder labels when the fitted object
#'   contains a supplied `Gamma` payload.
#' @export
gamma_matrix <- function(object, ...) {
  UseMethod("gamma_matrix")
}

#' @export
gamma_matrix.default <- function(object, ...) {
  stop(
    "`gamma_matrix()` requires an `hsquared_fit` object from an opt-in ",
    "metafounder or `H^Gamma` model.",
    call. = FALSE
  )
}

#' @export
gamma_matrix.hsquared_fit <- function(object, ...) {
  gamma <- object$payload$Gamma
  if (is.null(gamma)) {
    stop(
      "This `hsquared_fit` object does not contain a supplied metafounder ",
      "`Gamma` matrix. Use it only for opt-in metafounder or `H^Gamma` fits.",
      call. = FALSE
    )
  }
  gamma <- as.matrix(gamma)
  storage.mode(gamma) <- "double"
  labels <- object$payload$gamma_labels
  if (
    !is.null(labels) &&
      length(labels) == nrow(gamma) &&
      length(labels) == ncol(gamma)
  ) {
    dimnames(gamma) <- list(as.character(labels), as.character(labels))
  }
  gamma
}

#' Extract supplied metafounder group assignments
#'
#' `r lifecycle::badge("experimental")`
#'
#' `metafounder_groups()` returns the ID-keyed metafounder group assignments
#' carried by an experimental metafounder or `H^Gamma` single-step
#' `hsquared_fit` object. These assignments are provenance for the supplied
#' relationship, not estimated grouping parameters.
#'
#' @inheritParams variance_components
#'
#' @return A data frame with columns `id`, `metafounder_group`, and
#'   `is_metafounder`. Rows are **animals**, keyed by `id`; `is_metafounder`
#'   flags animals whose unknown-parent slot carries a supplied group label,
#'   not metafounders themselves. The metafounders are the pseudo-populations
#'   that are the rows and columns of [gamma_matrix()].
#' @export
metafounder_groups <- function(object, ...) {
  UseMethod("metafounder_groups")
}

#' @export
metafounder_groups.default <- function(object, ...) {
  stop(
    "`metafounder_groups()` requires an `hsquared_fit` object from an opt-in ",
    "metafounder or `H^Gamma` model.",
    call. = FALSE
  )
}

#' @export
metafounder_groups.hsquared_fit <- function(object, ...) {
  group_of <- object$payload$group_of
  if (is.null(group_of)) {
    stop(
      "This `hsquared_fit` object does not contain supplied metafounder ",
      "group assignments. Use it only for opt-in metafounder or `H^Gamma` ",
      "fits.",
      call. = FALSE
    )
  }
  group_of <- as.character(group_of)
  ids <- object$payload$ids %||% names(group_of)
  if (is.null(ids) || length(ids) != length(group_of)) {
    stop(
      "This `hsquared_fit` object has malformed metafounder group metadata; ",
      "the ID and group vectors do not align.",
      call. = FALSE
    )
  }
  data.frame(
    id = as.character(ids),
    metafounder_group = ifelse(nzchar(group_of), group_of, NA_character_),
    is_metafounder = nzchar(group_of),
    stringsAsFactors = FALSE
  )
}

#' Reserved metafounder effect extractor
#'
#' `r lifecycle::badge("experimental")`
#'
#' `metafounder_effects()` reserves the extractor name for future fitted
#' metafounder solutions. Current experimental metafounder and `H^Gamma` fits
#' expose supplied `Gamma` and group-assignment provenance via [gamma_matrix()]
#' and [metafounder_groups()], but the engine does not yet return explicit
#' combined-system metafounder effects for extraction.
#'
#' @inheritParams variance_components
#'
#' @return This extractor currently errors for all objects.
#' @export
metafounder_effects <- function(object, ...) {
  UseMethod("metafounder_effects")
}

#' @export
metafounder_effects.default <- function(object, ...) {
  stop(
    "`metafounder_effects()` requires an `hsquared_fit` object from a future ",
    "metafounder result that returns explicit metafounder solutions.",
    call. = FALSE
  )
}

#' @export
metafounder_effects.hsquared_fit <- function(object, ...) {
  stop(
    "`metafounder_effects()` is reserved but not implemented. Current ",
    "experimental metafounder and `H^Gamma` fits expose supplied provenance ",
    "through `gamma_matrix()` and `metafounder_groups()` only; the engine does ",
    "not yet return explicit metafounder solutions for this extractor.",
    call. = FALSE
  )
}

#' Reserved factor-analytic and G-matrix extractors
#'
#' `r lifecycle::badge("experimental")`
#'
#' These extractor names are reserved for future factor-analytic G-matrix
#' results. The current package can report invariant covariance and correlation
#' matrices from multivariate fits, but it does not yet expose
#' interpreted loadings, uniqueness/specific variance, or latent breeding
#' values. Loading columns (`genetic_loadings()`) and latent breeding values
#' (`latent_breeding_values()`) are rotation-nonunique until a rotation or
#' constraint policy is validated. Future `hsquared_fit` methods reserve
#' `effect` and rotation controls, but these controls currently error rather
#' than implying that loading axes are interpretable. The **rotation-invariant**
#' genetic eigenstructure and evolvability geometry are available now via
#' [eigen_G()] and the [g_matrix_geometry] family.
#'
#' `specific_variance()` (`Psi`, the factor-analytic specific/unique
#' variances) is a **different case**: for `G = Lambda Lambda' + Psi`,
#' rotating `Lambda -> Lambda Q` with `QQ' = I` leaves `Lambda Q Q' Lambda' =
#' Lambda Lambda'` unchanged, so `Psi` is rotation-**invariant** and
#' identified, not rotation-nonunique. The engine payload for structured
#' multivariate fits carries `Psi` as `genetic_uniqueness`, explicitly marked
#' identified, alongside the excluded, rotation-nonidentified loadings. This
#' extractor still errors on the R surface, but for a different reason:
#' `genetic_structure = "factor_analytic"` (and `"lowrank"`) are planned on
#' the R surface and are not yet activated on the R-to-Julia bridge, not
#' because `Psi` is unidentified.
#'
#' @inheritParams variance_components
#'
#' @return These reserved extractors currently error for `hsquared_fit` objects.
#' @name factor_g_extractors
NULL

#' @rdname factor_g_extractors
#' @export
genetic_loadings <- function(object, ...) {
  UseMethod("genetic_loadings")
}

#' @export
genetic_loadings.default <- function(object, ...) {
  hs_factor_g_extractor_default("genetic_loadings")
}

#' @export
genetic_loadings.hsquared_fit <- function(
  object,
  effect = "animal",
  rotate = NULL,
  ...
) {
  hs_factor_g_extractor_planned(
    "genetic_loadings",
    "factor-analytic G-matrix loadings",
    effect = effect,
    rotate = rotate
  )
}

#' @rdname factor_g_extractors
#' @export
specific_variance <- function(object, ...) {
  UseMethod("specific_variance")
}

#' @export
specific_variance.default <- function(object, ...) {
  hs_factor_g_extractor_default("specific_variance")
}

#' @export
specific_variance.hsquared_fit <- function(object, effect = "animal", ...) {
  hs_factor_g_extractor_planned(
    "specific_variance",
    "factor-analytic G-matrix uniqueness / specific variance",
    effect = effect
  )
}

#' @rdname factor_g_extractors
#' @export
latent_breeding_values <- function(object, ...) {
  UseMethod("latent_breeding_values")
}

#' @export
latent_breeding_values.default <- function(object, ...) {
  hs_factor_g_extractor_default("latent_breeding_values")
}

#' @export
latent_breeding_values.hsquared_fit <- function(
  object,
  effect = "animal",
  ...
) {
  hs_factor_g_extractor_planned(
    "latent_breeding_values",
    "latent breeding values from a factor-analytic G matrix",
    effect = effect
  )
}

#' @rdname g_matrix_geometry
#' @export
eigen_G <- function(object, ...) {
  UseMethod("eigen_G")
}

#' @export
eigen_G.default <- function(object, ...) {
  stop("`eigen_G()` requires an `hsquared_fit` object.", call. = FALSE)
}

#' @export
eigen_G.hsquared_fit <- function(object, ...) {
  hs_g_eigen(hs_fit_genetic_G(object))
}

hs_factor_g_extractor_default <- function(name) {
  stop(
    "`",
    name,
    "()` requires an `hsquared_fit` object from a future validated ",
    "factor-analytic or G-matrix result. The current package reserves this ",
    "extractor name but does not expose interpreted factor-analytic G-matrix ",
    "outputs yet.",
    call. = FALSE
  )
}

hs_factor_g_extractor_planned <- function(
  name,
  quantity,
  effect,
  rotate = NULL
) {
  if (!identical(effect, "animal")) {
    stop(
      "`",
      name,
      "()` currently reserves only `effect = \"animal\"`; other effects are ",
      "planned, not implemented.",
      call. = FALSE
    )
  }
  if (!is.null(rotate)) {
    stop(
      "`",
      name,
      "()` rotation controls are planned, not implemented. Loading axes need ",
      "a validated rotation or constraint policy before interpretation.",
      call. = FALSE
    )
  }
  if (identical(name, "specific_variance")) {
    stop(
      "`",
      name,
      "()` for ",
      quantity,
      " is planned, not implemented for `hsquared_fit` objects. Unlike ",
      "loading axes, `Psi` is rotation-INVARIANT and identified (for ",
      "`G = Lambda Lambda' + Psi`, rotating `Lambda -> Lambda Q` leaves ",
      "`Psi` unchanged). It is withheld because `genetic_structure = ",
      "\"factor_analytic\"` (and `\"lowrank\"`) are planned on the R surface ",
      "and not yet activated on the R-to-Julia bridge, not because `Psi` is ",
      "unidentified. Current multivariate fits report invariant ",
      "`genetic_covariance()` and `genetic_correlation()`.",
      call. = FALSE
    )
  }
  stop(
    "`",
    name,
    "()` for ",
    quantity,
    " is planned, not implemented for `hsquared_fit` objects. Current ",
    "multivariate fits report invariant `genetic_covariance()` and ",
    "`genetic_correlation()`; loading axes remain rotation-nonunique until ",
    "validated.",
    call. = FALSE
  )
}

#' Extract repeatability estimates
#'
#' `r lifecycle::badge("experimental")`
#'
#' `repeatability()` reports the repeatability `R = (Va + Vpe) / Vp` of the
#' opt-in, experimental repeatability (permanent-environment) model. It works
#' for `hsquared_fit` objects fitted with
#' `engine_control = list(target = "repeatability")`.
#'
#' Separating `sigma_a2` from `sigma_pe2` needs repeated records per
#' individual; see `?hs_control` (`target = "repeatability"`) for that
#' identifiability requirement.
#'
#' @inheritParams variance_components
#'
#' @return Repeatability results for repeatability `hsquared_fit` objects.
#' @export
repeatability <- function(object, ...) {
  UseMethod("repeatability")
}

#' @export
repeatability.default <- function(object, ...) {
  stop(
    "`repeatability()` requires an `hsquared_fit` object from the opt-in ",
    "repeatability model (`target = \"repeatability\"`).",
    call. = FALSE
  )
}

#' @export
repeatability.hsquared_fit <- function(object, ...) {
  hs_fit_result(object, "repeatability", "repeatability estimates")
}

#' Extract permanent-environment effects
#'
#' `r lifecycle::badge("experimental")`
#'
#' `permanent_effects()` returns the predicted permanent-environment effects of
#' the opt-in, experimental repeatability model.
#'
#' Separating `sigma_a2` from `sigma_pe2` needs repeated records per
#' individual; see `?hs_control` (`target = "repeatability"`) for that
#' identifiability requirement.
#'
#' @inheritParams variance_components
#'
#' @return Permanent-environment effect results for repeatability
#'   `hsquared_fit` objects.
#' @export
permanent_effects <- function(object, ...) {
  UseMethod("permanent_effects")
}

#' @export
permanent_effects.default <- function(object, ...) {
  stop(
    "`permanent_effects()` requires an `hsquared_fit` object from the opt-in ",
    "repeatability model (`target = \"repeatability\"`).",
    call. = FALSE
  )
}

#' @export
permanent_effects.hsquared_fit <- function(object, ...) {
  hs_fit_result(object, "permanent_effects", "permanent-environment effects")
}

#' Extract common-environment effects
#'
#' `r lifecycle::badge("experimental")`
#'
#' `common_env_effects()` returns the predicted common-environment effects of
#' the opt-in, experimental two-effect model (`target = "two_effect"`).
#'
#' @inheritParams variance_components
#'
#' @return Common-environment effect results for two-effect `hsquared_fit`
#'   objects.
#' @export
common_env_effects <- function(object, ...) {
  UseMethod("common_env_effects")
}

#' @export
common_env_effects.default <- function(object, ...) {
  stop(
    "`common_env_effects()` requires an `hsquared_fit` object from the opt-in ",
    "two-effect (common-environment) model (`target = \"two_effect\"`).",
    call. = FALSE
  )
}

#' @export
common_env_effects.hsquared_fit <- function(object, ...) {
  hs_fit_result(object, "common_env_effects", "common-environment effects")
}

#' Extract maternal genetic effects
#'
#' `r lifecycle::badge("experimental")`
#'
#' `maternal_effects()` returns the predicted maternal genetic effects of the
#' opt-in, experimental maternal two-effect model
#' (`target = "two_effect"` with a `maternal_genetic()` term).
#'
#' @inheritParams variance_components
#'
#' @return Maternal genetic effect results for maternal `hsquared_fit` objects.
#' @export
maternal_effects <- function(object, ...) {
  UseMethod("maternal_effects")
}

#' @export
maternal_effects.default <- function(object, ...) {
  stop(
    "`maternal_effects()` requires an `hsquared_fit` object from the opt-in ",
    "maternal two-effect model (`target = \"two_effect\"` with a ",
    "`maternal_genetic()` term).",
    call. = FALSE
  )
}

#' @export
maternal_effects.hsquared_fit <- function(object, ...) {
  hs_fit_result(object, "maternal_effects", "maternal genetic effects")
}

# Falconer interpretation fence attached to the two-effect variance-ratio
# accessors. c2/m2 is a variance ratio (the proportion of phenotypic variance
# from the common-environment / maternal effect), NOT a heritability. In the
# same fit `heritability()` reports narrow-sense h2 = sigma_a2 / (sigma_a2 +
# sigma_2 + sigma_e2) WITHIN this model. Carried as an attribute so the number
# never travels without its interpretation.
hs_two_effect_ratio_fence <- function(kind = c("common_env", "maternal")) {
  kind <- match.arg(kind)
  label <- if (kind == "common_env") {
    "common-environment"
  } else {
    "maternal (dam) genetic"
  }
  paste0(
    "This is a variance ratio: the proportion of phenotypic variance ",
    "attributable to the ",
    label,
    " effect (c2/m2), NOT a heritability. ",
    "In the same fit heritability() reports narrow-sense ",
    "h2 = sigma_a2 / (sigma_a2 + sigma_2 + sigma_e2) WITHIN this two-effect ",
    "model. At a variance-component boundary (sigma -> 0) the ratio is flagged ",
    "and its interval is NA, not a spuriously tight CI."
  )
}

#' Extract the common-environment variance ratio (c2)
#'
#' `r lifecycle::badge("experimental")`
#'
#' `common_env_proportion()` returns the estimated common-environment variance
#' ratio `c2 = sigma_c2 / (sigma_a2 + sigma_c2 + sigma_e2)` from the opt-in,
#' experimental two-effect model
#' (`target = "two_effect"` with a `common_env()` term).
#'
#' Falconer fence: `c2` is a **variance ratio** (the proportion of phenotypic
#' variance from the shared common-environment effect), **not a heritability**.
#' In the same fit, [heritability()] reports narrow-sense
#' `h2 = sigma_a2 / (sigma_a2 + sigma_c2 + sigma_e2)` *within this two-effect
#' model*. The returned data frame carries this note as an `"interpretation"`
#' attribute. For an interval on `c2`, see [common_env_proportion_interval()].
#'
#' @inheritParams variance_components
#'
#' @return A one-row data frame with `term` and `estimate`, plus an
#'   `"interpretation"` attribute, for two-effect `hsquared_fit` objects.
#' @export
common_env_proportion <- function(object, ...) {
  UseMethod("common_env_proportion")
}

#' @export
common_env_proportion.default <- function(object, ...) {
  stop(
    "`common_env_proportion()` requires an `hsquared_fit` object from the ",
    "opt-in two-effect (common-environment) model (`target = \"two_effect\"`).",
    call. = FALSE
  )
}

#' @export
common_env_proportion.hsquared_fit <- function(object, ...) {
  out <- hs_fit_result(
    object,
    "common_env_proportion",
    "a common-environment variance ratio (c2)"
  )
  attr(out, "interpretation") <- hs_two_effect_ratio_fence("common_env")
  out
}

#' Extract the maternal variance ratio (m2)
#'
#' `r lifecycle::badge("experimental")`
#'
#' `maternal_proportion()` returns the estimated maternal variance ratio
#' `m2 = sigma_m2 / (sigma_a2 + sigma_m2 + sigma_e2)` from the opt-in,
#' experimental maternal two-effect model
#' (`target = "two_effect"` with a `maternal_genetic()` term).
#'
#' Falconer fence: `m2` is a **variance ratio** (the proportion of phenotypic
#' variance from the maternal (dam) genetic effect), **not a heritability**. In
#' the same fit, [heritability()] reports the narrow-sense *direct*
#' `h2 = sigma_a2 / (sigma_a2 + sigma_m2 + sigma_e2)` *within this two-effect
#' model*; the direct and maternal genetic effects are modelled as
#' uncorrelated here (no direct-maternal genetic covariance). The returned data
#' frame carries this note as an `"interpretation"` attribute. For an interval
#' on `m2`, see [maternal_proportion_interval()].
#'
#' @inheritParams variance_components
#'
#' @return A one-row data frame with `term` and `estimate`, plus an
#'   `"interpretation"` attribute, for maternal two-effect `hsquared_fit`
#'   objects.
#' @export
maternal_proportion <- function(object, ...) {
  UseMethod("maternal_proportion")
}

#' @export
maternal_proportion.default <- function(object, ...) {
  stop(
    "`maternal_proportion()` requires an `hsquared_fit` object from the opt-in ",
    "maternal two-effect model (`target = \"two_effect\"` with a ",
    "`maternal_genetic()` term).",
    call. = FALSE
  )
}

#' @export
maternal_proportion.hsquared_fit <- function(object, ...) {
  out <- hs_fit_result(
    object,
    "maternal_proportion",
    "a maternal variance ratio (m2)"
  )
  attr(out, "interpretation") <- hs_two_effect_ratio_fence("maternal")
  out
}

#' Extract breeding values
#'
#' `r lifecycle::badge("experimental")`
#'
#' `breeding_values()` is part of the planned v0.1 fitted-object contract. It
#' works for `hsquared_fit` objects that contain a Julia result. `EBV()` and
#' `BLUP()` are aliases for applied quantitative-genetic workflows.
#'
#' @inheritParams variance_components
#'
#' @return Breeding value results for `hsquared_fit` objects.
#' @export
breeding_values <- function(object, ...) {
  UseMethod("breeding_values")
}

#' @export
breeding_values.default <- function(object, ...) {
  hs_breeding_values_default("breeding_values")
}

#' @export
breeding_values.hsquared_fit <- function(object, ...) {
  hs_fit_result(object, "breeding_values", "breeding values")
}

#' @rdname breeding_values
#' @export
EBV <- function(object, ...) {
  UseMethod("EBV")
}

#' @export
EBV.default <- function(object, ...) {
  hs_breeding_values_default("EBV")
}

#' @export
EBV.hsquared_fit <- function(object, ...) {
  breeding_values(object, ...)
}

#' @rdname breeding_values
#' @export
BLUP <- function(object, ...) {
  UseMethod("BLUP")
}

#' @export
BLUP.default <- function(object, ...) {
  hs_breeding_values_default("BLUP")
}

#' @export
BLUP.hsquared_fit <- function(object, ...) {
  breeding_values(object, ...)
}

hs_breeding_values_default <- function(name) {
  stop(
    "`",
    name,
    "()` requires an `hsquared_fit` object. The current package only returns ",
    "these from fitted `hsquared_fit` results.",
    call. = FALSE
  )
}

#' Extract prediction error variances
#'
#' `r lifecycle::badge("experimental")`
#'
#' `prediction_error_variance()` is part of the planned v0.1 fitted-object
#' contract. It returns values only when an `hsquared_fit` object contains a
#' Julia result field for prediction error variances.
#'
#' The engine defines PEV as the animal-block diagonal of the inverse of the
#' precision-scaled Henderson mixed-model-equation coefficient matrix, so
#' `diag(C^-1)` is directly the prediction error variance in variance units,
#' with no separate `sigma_e2` factor to apply.
#'
#' @inheritParams variance_components
#'
#' @return Prediction error variances for `hsquared_fit` objects.
#' @export
prediction_error_variance <- function(object, ...) {
  UseMethod("prediction_error_variance")
}

#' @export
prediction_error_variance.default <- function(object, ...) {
  stop(
    "`prediction_error_variance()` requires an `hsquared_fit` object. The ",
    "current package only returns these from fitted `hsquared_fit` results.",
    call. = FALSE
  )
}

#' @export
prediction_error_variance.hsquared_fit <- function(object, ...) {
  hs_fit_result(
    object,
    "prediction_error_variance",
    "prediction error variances"
  )
}

#' Extract reliability and accuracy estimates
#'
#' `r lifecycle::badge("experimental")`
#'
#' `reliability()` is part of the planned v0.1 fitted-object contract. It
#' returns values only when an `hsquared_fit` object contains a Julia result
#' field for reliability estimates. `accuracy()` returns the square root of
#' reliability for `hsquared_fit` objects.
#'
#' The engine defines reliability as `1 - PEV_i / (sigma_a2 * A_ii)`, using the
#' dense relationship matrix `A = inv(Ainv)` implied by the supplied precision.
#' For a genomic spec (`Ainv = Ginv`) this `A_ii` is `diag(inv(Ginv)) = diag(G)
#' + ridge` (the regularized genomic self-relationship, often not equal to 1),
#' so the genomic ridge perturbs the reported reliability and accuracy.
#' Reliability values are **not clipped**; `accuracy()` **rejects** a
#' reliability value below 0 or above 1 rather than silently clipping it.
#'
#' @inheritParams variance_components
#'
#' @return Reliability estimates for `hsquared_fit` objects.
#' @export
reliability <- function(object, ...) {
  UseMethod("reliability")
}

#' @export
reliability.default <- function(object, ...) {
  stop(
    "`reliability()` requires an `hsquared_fit` object. The current package ",
    "only returns these from fitted `hsquared_fit` results.",
    call. = FALSE
  )
}

#' @export
reliability.hsquared_fit <- function(object, ...) {
  hs_fit_result(object, "reliability", "reliability estimates")
}

#' @rdname reliability
#' @export
accuracy <- function(object, ...) {
  UseMethod("accuracy")
}

#' @export
accuracy.default <- function(object, ...) {
  stop(
    "`accuracy()` requires an `hsquared_fit` object with reliability ",
    "estimates.",
    call. = FALSE
  )
}

#' @export
accuracy.hsquared_fit <- function(object, ...) {
  rel <- reliability(object, ...)
  if (!is.data.frame(rel) || !"value" %in% names(rel)) {
    stop(
      "Reliability results must be a data frame with a `value` column to ",
      "compute accuracy.",
      call. = FALSE
    )
  }
  if (any(is.na(rel$value)) || any(rel$value < 0) || any(rel$value > 1)) {
    stop(
      "Reliability values must be between 0 and 1 to compute accuracy.",
      call. = FALSE
    )
  }
  out <- rel
  out$value <- sqrt(out$value)
  out
}

#' Extract an experimental heritability confidence interval
#'
#' `r lifecycle::badge("experimental")`
#'
#' `heritability_interval()` returns an **experimental** large-sample confidence
#' interval for `h^2`. It is available only when an `hsquared_fit` object
#' contains the interval field, which the default Gaussian animal-model fit
#' (`engine = "fit"`) populates from the engine's
#' `HSquared.heritability_interval()` when a local Julia engine is present and
#' the estimate is interior to `(0, 1)`. On the opt-in two-effect fit it returns
#' the direct-heritability ratio interval (`ratio1`), and on the opt-in
#' multi-effect fit (`target = "multi_effect"`, K >= 3 blocks) it returns the
#' ANIMAL block's ratio interval (the animal additive variance over the total
#' phenotypic variance); the other blocks' variance-ratio intervals are surfaced
#' separately in `fit$result$variance_ratio_intervals`.
#'
#' The interval leg is a REML-only, asymptotic (logit delta-method or profile)
#' approximation returned by the engine. The 2000-rep C1 coverage confirm
#' (HSquared.jl DRAC job **47925485**) places the univariate pedigree h²
#' interval at the **directional-conservative** claim level under doc-34 §4:
#' delta over-covers (worst Ĉ 0.969), profile is in-band (worst Ĉ 0.950), and
#' bootstrap **mildly under-covers versus nominal 0.95** at the governing
#' interior cell (worst Ĉ **0.924** at h²=0.5, still at least 0.90). Those
#' figures were measured on one **interpretable small design** (`q = 120`),
#' at **interior** `h² in {0.3, 0.5, 0.7}`, at the **0.95** level, in DRAC job
#' **47925485**; **unbalanced designs were not in the confirm grid**, so
#' these numbers say nothing about coverage under imbalance. That is **not**
#' coverage-calibrated at nominal and **not** "never under-covers." Profile and
#' bootstrap h² interval legs were measured in the same confirm but are **not**
#' separately surfaced by this extractor; only the engine-returned interval is
#' shown. It is reported as a point estimate plus bounds, not a validated
#' (coverage-calibrated) capability, and remains unreliable near the
#' `h² → 0` boundary. Genomic fits still error (scale-labelled interval not
#' validated). The underlying estimators `V3-TWOEFFECT-REML` /
#' `V3-NEFFECT-REML` are `covered`, but this **interval** is not.
#'
#' @inheritParams variance_components
#'
#' @return A one-row data frame with `estimate`, `lower`, `upper`, `level`, `se`
#'   (`NA` for the profile method), and `method`, for `hsquared_fit` objects that
#'   contain it.
#' @export
heritability_interval <- function(object, ...) {
  UseMethod("heritability_interval")
}

#' @export
heritability_interval.default <- function(object, ...) {
  stop(
    "`heritability_interval()` requires an `hsquared_fit` object. The current ",
    "package only returns these from fitted `hsquared_fit` results.",
    call. = FALSE
  )
}

#' @export
heritability_interval.hsquared_fit <- function(object, ...) {
  if (hs_fit_is_genomic(object)) {
    stop(
      "`heritability_interval()` is not available for genomic fits. The ",
      "genomic variance ratio is defined on the declared relationship scale; ",
      "a scale-labelled, calibrated interval has not been validated.",
      call. = FALSE
    )
  }
  hs_fit_result(
    object,
    "heritability_interval",
    "an experimental heritability confidence interval"
  )
}

#' Extract experimental variance-component and heritability standard errors
#'
#' `r lifecycle::badge("experimental")`
#'
#' `variance_component_standard_errors()` and `heritability_standard_error()`
#' return **experimental** large-sample (delta-method) standard errors derived
#' from the REML average-information matrix. They are available only when an
#' `hsquared_fit` object contains them; the default Gaussian animal-model fit
#' populates them from the engine when a local Julia engine is present and the
#' AI matrix is invertible. When the engine could not compute them, `hsquared()`
#' warns once at fit time and keeps the engine's reason in
#' `attr(fit, "bridge_errors")`, a character vector named by the engine
#' function that failed (HSquared.jl#351); the field is then absent, not `NA`.
#'
#' These mirror the engine row `V1-HERIT-CI` (`partial`): asymptotic,
#' REML-only, and unreliable at small `n` or near a variance-component
#' boundary (where the AI matrix is ill-conditioned and the fields are
#' omitted). The 2000-rep C1 coverage confirm (job **47925485**) measured the
#' **additive-variance delta/Wald** interval implied by these SEs to **under-cover**
#' (0.897 at nominal 0.95, h²=0.5), placing it at **experimental-only**:
#' the SE is a point-estimate reference only, **not a calibrated and not a
#' conservative interval**. (A profile additive-variance interval is
#' `directional-conservative` by the same run — including the h²=0.3 cell
#' 0.963 — but is not surfaced by this delta-method extractor. Additive-variance
#' bootstrap worst-cell Ĉ 0.918 is the same DC mild-under branch as h²
#' bootstrap.) Not coverage-calibrated, not a validated capability.
#'
#' @inheritParams variance_components
#'
#' @return `variance_component_standard_errors()` returns a data frame with
#'   `component` and `se`. `heritability_standard_error()` returns a one-row
#'   data frame with `term` and `se`, matching [heritability()] so the two
#'   merge on `term` (hsquared#236). Only for `hsquared_fit` objects that
#'   contain them. The fit object may still store `heritability_se` as a bare
#'   numeric internally; this extractor is the public shape.
#' @export
variance_component_standard_errors <- function(object, ...) {
  UseMethod("variance_component_standard_errors")
}

#' @export
variance_component_standard_errors.default <- function(object, ...) {
  stop(
    "`variance_component_standard_errors()` requires an `hsquared_fit` object. ",
    "The current package only returns these from fitted `hsquared_fit` results.",
    call. = FALSE
  )
}

#' @export
variance_component_standard_errors.hsquared_fit <- function(object, ...) {
  hs_fit_result(
    object,
    "variance_component_se",
    "experimental variance-component standard errors"
  )
}

#' @rdname variance_component_standard_errors
#' @export
heritability_standard_error <- function(object, ...) {
  UseMethod("heritability_standard_error")
}

#' @export
heritability_standard_error.default <- function(object, ...) {
  stop(
    "`heritability_standard_error()` requires an `hsquared_fit` object. The ",
    "current package only returns these from fitted `hsquared_fit` results.",
    call. = FALSE
  )
}

#' @export
heritability_standard_error.hsquared_fit <- function(object, ...) {
  if (hs_fit_is_genomic(object)) {
    stop(
      "`heritability_standard_error()` is not available for genomic fits. ",
      "A scale-labelled genomic-ratio standard error has not been validated.",
      call. = FALSE
    )
  }
  se <- hs_fit_result(
    object,
    "heritability_se",
    "an experimental heritability standard error"
  )
  hs_as_heritability_se_frame(object, se)
}

# Public SE shape for heritability (hsquared#236): always `data.frame(term, se)`,
# matching `heritability()` so `merge(..., by = "term")` works. Internal storage
# may still be a bare numeric (animal-model and sparse-repeatability bridges).
hs_as_heritability_se_frame <- function(object, se) {
  if (is.data.frame(se)) {
    if (!all(c("term", "se") %in% names(se))) {
      stop(
        "`heritability_se` on the fit must be a numeric or a data frame with ",
        "columns `term` and `se`.",
        call. = FALSE
      )
    }
    out <- se[, c("term", "se"), drop = FALSE]
    out$term <- as.character(out$term)
    out$se <- as.numeric(out$se)
    rownames(out) <- NULL
    return(out)
  }
  term <- "animal"
  h2 <- object$result$heritability
  if (is.data.frame(h2) && "term" %in% names(h2) && nrow(h2) >= 1L) {
    if (nrow(h2) == 1L) {
      term <- as.character(h2$term[[1L]])
    } else if ("animal" %in% h2$term) {
      term <- "animal"
    } else {
      term <- as.character(h2$term[[1L]])
    }
  }
  data.frame(
    term = term,
    se = as.numeric(se)[[1L]],
    stringsAsFactors = FALSE
  )
}

#' Extract an experimental repeatability confidence interval
#'
#' `r lifecycle::badge("experimental")`
#'
#' `repeatability_interval()` returns an **experimental** large-sample (logit
#' delta-method) confidence interval for the repeatability coefficient
#' `t = (Va + Vpe) / Vp` of the opt-in repeatability (permanent-environment)
#' model, available only when the fit contains it.
#'
#' It mirrors the engine row `V3-REPEAT-REML` (`partial`): the engine's
#' repeatability REML estimator and this interval are engine-internal
#' self-consistency tested (recovery of `t` and interval bracketing / range /
#' level-nesting / point-estimate match on seeded fixtures), but there is no
#' external comparator, no `h^2` interval, and no deep-pedigree validation. It is
#' asymptotic, REML-only, unreliable at small `n` or near the (0, 1) boundary
#' (where the engine throws and the field is omitted), and not a validated
#' capability.
#'
#' @inheritParams variance_components
#'
#' @return A one-row data frame with `estimate` (the repeatability `t`), `lower`,
#'   `upper`, `level`, and `se`, for `hsquared_fit` objects that contain it.
#' @export
repeatability_interval <- function(object, ...) {
  UseMethod("repeatability_interval")
}

#' @export
repeatability_interval.default <- function(object, ...) {
  stop(
    "`repeatability_interval()` requires an `hsquared_fit` object from the ",
    "opt-in repeatability model. The current package only returns these from ",
    "fitted `hsquared_fit` results.",
    call. = FALSE
  )
}

#' @export
repeatability_interval.hsquared_fit <- function(object, ...) {
  hs_fit_result(
    object,
    "repeatability_interval",
    "an experimental repeatability confidence interval"
  )
}

#' Extract an experimental common-environment / maternal variance-ratio interval
#'
#' `r lifecycle::badge("experimental")`
#'
#' `common_env_proportion_interval()` and `maternal_proportion_interval()`
#' return an **experimental** large-sample (logit delta-method) confidence
#' interval for the second variance ratio (`c2` / `m2`, `ratio2`) of the opt-in
#' two-effect model, available only when the fit contains it. On the same
#' two-effect fit, [heritability_interval()] returns the matching interval for
#' the direct heritability (`h2`, `ratio1`).
#'
#' This mirrors the engine row `V3-TWOEFFECT-REML`: the interval is the asymptotic
#' delta-method CI built from the two-effect REML observed information (the
#' finite-difference Hessian of the two-effect REML log-likelihood at the
#' optimum). It is **asymptotic, delta-method, REML only, and NOT
#' coverage-calibrated** -- on small samples the REML surface is flat and the
#' interval is unreliable (the parametric bootstrap is the only finite-sample-
#' aware path). No calibrated coverage is claimed.
#'
#' Boundary honesty: when the ratio's variance component sits on the boundary
#' (`sigma -> 0`) it is flagged (`boundary = TRUE`) and `lower`/`upper` are
#' `NA`, never a spuriously tight CI. `lower_clamped`/`upper_clamped` flag when
#' a bound reaches the numerical `(0, 1)` rails.
#'
#' Falconer fence: `c2` / `m2` is a **variance ratio**, not a heritability (see
#' [common_env_proportion()] / [maternal_proportion()]).
#'
#' @inheritParams variance_components
#'
#' @return A one-row data frame with `estimate` (the ratio), `lower`, `upper`,
#'   `level`, `se`, `lower_clamped`, `upper_clamped`, and `boundary`, plus an
#'   `"interpretation"` attribute, for two-effect `hsquared_fit` objects that
#'   contain it.
#' @export
common_env_proportion_interval <- function(object, ...) {
  UseMethod("common_env_proportion_interval")
}

#' @export
common_env_proportion_interval.default <- function(object, ...) {
  stop(
    "`common_env_proportion_interval()` requires an `hsquared_fit` object from ",
    "the opt-in two-effect (common-environment) model.",
    call. = FALSE
  )
}

#' @export
common_env_proportion_interval.hsquared_fit <- function(object, ...) {
  out <- hs_fit_result(
    object,
    "common_env_proportion_interval",
    "an experimental common-environment variance-ratio interval"
  )
  attr(out, "interpretation") <- hs_two_effect_ratio_fence("common_env")
  out
}

#' @rdname common_env_proportion_interval
#' @export
maternal_proportion_interval <- function(object, ...) {
  UseMethod("maternal_proportion_interval")
}

#' @export
maternal_proportion_interval.default <- function(object, ...) {
  stop(
    "`maternal_proportion_interval()` requires an `hsquared_fit` object from ",
    "the opt-in maternal two-effect model.",
    call. = FALSE
  )
}

#' @export
maternal_proportion_interval.hsquared_fit <- function(object, ...) {
  out <- hs_fit_result(
    object,
    "maternal_proportion_interval",
    "an experimental maternal variance-ratio interval"
  )
  attr(out, "interpretation") <- hs_two_effect_ratio_fence("maternal")
  out
}

#' Extract experimental multivariate covariance standard errors
#'
#' `r lifecycle::badge("experimental")`
#'
#' `covariance_standard_errors()` returns **experimental** large-sample
#' (delta-method) standard errors for the multivariate genetic/residual
#' covariance and correlation matrices and per-trait `h^2`, for an opt-in
#' **unstructured** multivariate fit, when the engine returned them.
#'
#' Heavy caveats (engine multivariate REML validation row, `partial`): the strict per-seed
#' recovery gate is still a non-pass (7/12 unstructured seeds in the updated
#' study), but a 12-seed bias/MCSE study (twin `HSquared.jl#78`) shows **no
#' detectable bias** -- all six covariance parameters have `|bias| <= 2*MCSE`
#' (largest 0.84*MCSE), a low-power non-rejection consistent with an unbiased
#' estimator (not a proof), with EBV accuracy ~ 0.90 in both traits; a cold-start
#' replication (#79) reaches the same optimum on all 12 seeds, so it is not a
#' warm-start artifact. The per-seed gate failures reflect sampling variance of
#' the estimated `G` at this design, not a detected bias. These SEs remain
#' asymptotic, REML-only, **unstructured-only** (the engine refuses structured /
#' factor-analytic fits, whose loadings are rotation-nonidentified), omitted at a
#' flat/boundary optimum, not coverage-calibrated, with no external comparator,
#' and not a validated capability.
#'
#' The returned `genetic_correlation` SE is on the correlation (`r`) scale. A
#' naive symmetric interval built as `estimate +/- 1.96 * se` on that scale
#' can leave `(-1, 1)` for strong correlations. The engine instead builds its own
#' genetic-correlation interval on the **Fisher-z scale**
#' (`z = atanh(r)`, `se_z = se_r / (1 - r^2)`, endpoints
#' `tanh(z -+ q * se_z)`), which always lands inside `(-1, 1)`. That
#' range-respecting endpoint is a separate property from calibration: the
#' Fisher-z transform does not make these SEs coverage-calibrated (see above).
#'
#' @inheritParams variance_components
#'
#' @return A named list of standard-error matrices `genetic_covariance`,
#'   `residual_covariance`, `genetic_correlation`, `residual_correlation`, and a
#'   per-trait `heritability` SE vector, for `hsquared_fit` objects that contain
#'   them.
#' @export
covariance_standard_errors <- function(object, ...) {
  UseMethod("covariance_standard_errors")
}

#' @export
covariance_standard_errors.default <- function(object, ...) {
  stop(
    "`covariance_standard_errors()` requires an `hsquared_fit` object from the ",
    "unstructured multivariate model. The current package only returns ",
    "these from fitted `hsquared_fit` results.",
    call. = FALSE
  )
}

#' @export
covariance_standard_errors.hsquared_fit <- function(object, ...) {
  hs_fit_result(
    object,
    "covariance_standard_errors",
    "experimental multivariate covariance standard errors"
  )
}

#' Likelihood-ratio test for genetic covariance structure
#'
#' `r lifecycle::badge("experimental")`
#'
#' `covariance_structure_lrt(constrained, full)` is an **experimental** nested
#' likelihood-ratio test comparing two multivariate fits **on the same
#' data**: a `constrained` genetic structure (currently `genetic_structure =
#' "diagonal"`) against the `full` `"unstructured"` fit. The statistic is
#' `2 * (logLik(full) - logLik(constrained))` on
#' `df = n_genetic_params(full) - n_genetic_params(constrained)`.
#'
#' For diagonal-vs-unstructured the null (off-diagonal genetic covariances = 0)
#' is interior, so the chi^2 reference is exact (`boundary = FALSE`). Structures
#' whose null lies on a rank/PSD boundary (low-rank / factor-analytic) would need
#' a chi^2-mixture correction and are gated out of the R bridge for now.
#'
#' It mirrors the engine multivariate REML validation row (`partial`): asymptotic, REML-only,
#' dense validation-scale, with the multivariate recovery calibration not yet
#' passed -- a reported test, not a validated one. Both fits must be on the same
#' response, fixed effects, and pedigree.
#'
#' @param constrained,full Two `hsquared_fit` objects from the multivariate
#'   model; `full` must nest `constrained` (more genetic
#'   covariance parameters).
#' @param ... Unused.
#'
#' @return A one-row data frame with `statistic`, `df`, `pvalue`, `boundary`, and
#'   the `constrained` / `full` genetic-structure labels.
#' @export
covariance_structure_lrt <- function(constrained, full, ...) {
  if (
    !inherits(constrained, "hsquared_fit") || !inherits(full, "hsquared_fit")
  ) {
    stop(
      "`constrained` and `full` must both be `hsquared_fit` objects from the ",
      "multivariate model.",
      call. = FALSE
    )
  }
  field <- function(fit, nm) {
    v <- fit$result[[nm]]
    if (is.null(v)) {
      stop(
        "This `hsquared_fit` does not carry `",
        nm,
        "`. The covariance-structure LRT needs converged multivariate fits ",
        "that report `loglik` and `n_genetic_params`.",
        call. = FALSE
      )
    }
    v
  }
  ll_c <- as.numeric(field(constrained, "loglik"))
  ll_f <- as.numeric(field(full, "loglik"))
  np_c <- as.integer(field(constrained, "n_genetic_params"))
  np_f <- as.integer(field(full, "n_genetic_params"))
  df <- np_f - np_c
  if (df <= 0L) {
    stop(
      "`full` must have more genetic covariance parameters than `constrained` ",
      "(df = ",
      df,
      "); call as `covariance_structure_lrt(constrained, full)`.",
      call. = FALSE
    )
  }
  sc <- constrained$result$genetic_structure %||% NA_character_
  sf <- full$result$genetic_structure %||% NA_character_
  stat <- 2 * (ll_f - ll_c)
  boundary <- !(identical(sc, "diagonal") && identical(sf, "unstructured"))
  data.frame(
    statistic = stat,
    df = df,
    pvalue = stats::pchisq(max(stat, 0), df = df, lower.tail = FALSE),
    boundary = boundary,
    constrained = sc,
    full = sf,
    stringsAsFactors = FALSE
  )
}

#' Inspect fitted-model diagnostics
#'
#' `r lifecycle::badge("experimental")`
#'
#' `fit_diagnostics()` returns a compact diagnostics table for an
#' `hsquared_fit` object. It is an inspection helper over the current result
#' payload: it does not refit the model, rerun validation checks, or promote an
#' experimental bridge target to production support.
#'
#' Two rows share the word "boundary" but report unrelated things
#' (hsquared#230):
#'
#' * `at_boundary` / `at_boundary_condition` derive from the fitted
#'   **estimates**: a variance component at or near zero, or an inadmissible
#'   negative estimate (`hs_fit_boundary_class()`). Emitted for the Gaussian
#'   animal-model / genomic / single-step routes that report one primary
#'   variance-share layout; the `_condition` row is dropped (not `NA`) when
#'   the fit is interior, since there is nothing extra to say.
#' * `search_boundary` / `search_boundary_condition` mirror the non-Gaussian
#'   bridge's **optimizer** state: `result$boundary` (documented as
#'   `fit$result$boundary` in `?hs_control`) reports whether
#'   `HSquared.jl`'s single-variance Brent search stopped on the rail of its
#'   log-scale bracket (`log(initial$sigma_a2) +/- 6`, HSquared.jl#327) or,
#'   under `restart_check = TRUE`, whether the two-start gap fence fired;
#'   either way the estimate is a function of the search, not the data.
#'   Only fits whose result carries a `boundary` field (the non-Gaussian
#'   route) gain this pair; a Gaussian `engine = "fit"` fit gains neither row
#'   (not even `NA`). Unlike `at_boundary_condition`, `search_boundary_condition`
#'   is always emitted alongside a present `search_boundary` flag, reading
#'   `"interior"` rather than disappearing when the flag is `FALSE` -- the
#'   search-bracket rail is a binary optimizer state with no benign-vs-
#'   inadmissible split to only sometimes report. In practice a
#'   `search_boundary = TRUE` fit is refused before it reaches an
#'   `hsquared_fit` object at all, so a returned fit's `search_boundary` is
#'   always `FALSE`; `TRUE` is reachable only on a hand-built `result` object.
#'
#' @inheritParams variance_components
#'
#' @return A data frame with `metric` and `value` columns and class
#'   `"hs_fit_diagnostics"`.
#' @export
fit_diagnostics <- function(object, ...) {
  UseMethod("fit_diagnostics")
}

#' @export
fit_diagnostics.default <- function(object, ...) {
  stop(
    "`fit_diagnostics()` requires an `hsquared_fit` object.",
    call. = FALSE
  )
}

#' @export
fit_diagnostics.hsquared_fit <- function(object, ...) {
  diagnostics <- object$result$diagnostics %||% list()
  if (!is.list(diagnostics)) {
    diagnostics <- list(diagnostics = diagnostics)
  }
  base <- list(
    engine = object$engine,
    method = object$spec$method %||% diagnostics$method,
    family = object$spec$family$family,
    target = object$spec$target %||%
      diagnostics$target %||%
      "variance_components",
    converged = object$result$converged %||% diagnostics$converged,
    optimizer_status = diagnostics$optimizer_status,
    iterations = diagnostics$iterations,
    loglik = object$result$loglik,
    df = object$result$df,
    nobs = object$result$nobs %||%
      if (!is.null(object$payload$y)) length(object$payload$y) else NULL,
    dense_validation_path = diagnostics$dense_validation_path,
    variance_components_source = diagnostics$variance_components,
    at_boundary = hs_fit_boundary_flag(object),
    at_boundary_condition = hs_fit_boundary_condition_label(
      hs_fit_boundary_class(object)
    ),
    search_boundary = hs_fit_search_boundary_flag(object),
    search_boundary_condition = hs_fit_search_boundary_condition_label(
      hs_fit_search_boundary_flag(object)
    )
  )

  diagnostic_names <- names(diagnostics)
  already_reported <- c(
    "method",
    "target",
    "converged",
    "optimizer_status",
    "iterations",
    "dense_validation_path",
    "variance_components",
    "search_boundary",
    "search_boundary_condition"
  )
  extras <- diagnostics[setdiff(diagnostic_names, already_reported)]
  rows <- c(base, extras)
  rows <- rows[!vapply(rows, is.null, logical(1))]

  out <- data.frame(
    metric = names(rows),
    value = vapply(rows, hs_diagnostic_value, character(1)),
    stringsAsFactors = FALSE
  )
  class(out) <- c("hs_fit_diagnostics", class(out))
  out
}

#' @export
print.hs_fit_diagnostics <- function(x, ...) {
  cat("<hs_fit_diagnostics>\n")
  out <- x
  class(out) <- setdiff(class(out), "hs_fit_diagnostics")
  print.data.frame(out, row.names = FALSE)
  invisible(x)
}

# Classify the variance-component boundary condition, distinguishing a benign
# near-zero boundary from an inadmissible NEGATIVE variance estimate. Returns
# one of "zero" (0 <= min share <= tol), "negative" (min share < 0, an
# inadmissible fit, not a clean boundary), FALSE (interior), or NULL when the
# components are unavailable or the layout is not a single-primary one. A NULL is
# also returned when the variance shape is not a list (e.g. an atomic vector or
# matrix), so a malformed payload drops the row rather than crashing the caller
# with "$ operator is invalid for atomic vectors" (a data.frame IS a list, so
# the normal path proceeds). The negative case is surfaced separately because a
# negative variance is inadmissible: it must not be read as the at/near-zero
# boundary that the v0.1 contract (item 4) treats as a clean h2 -> 0 / 1 edge.
hs_fit_boundary_class <- function(object, tol = 1e-4) {
  vc <- object$result$variance_components
  if (!is.list(vc)) {
    return(NULL)
  }
  if (is.null(vc$estimate) || is.null(vc$component)) {
    return(NULL)
  }
  est <- as.numeric(vc$estimate)
  primary_names <- c("animal", "genomic", "single_step", "metafounder")
  primary <- est[vc$component %in% primary_names]
  if (length(primary) != 1L) {
    return(NULL)
  }
  if (any(!is.finite(est))) {
    return(NULL)
  }
  # A negative component is inadmissible regardless of the total scale, so check
  # the raw sign before normalizing (a negative estimate can drive the total
  # non-positive, which would otherwise drop the row).
  if (min(est) < 0) {
    return("negative")
  }
  total <- sum(est)
  if (total <= 0) {
    return(NULL)
  }
  shares <- est / total
  if (min(shares) <= tol) {
    "zero"
  } else {
    FALSE
  }
}

# Logical boundary flag for the `at_boundary` diagnostics row: TRUE when the fit
# sits at a variance-component boundary. A near-zero ("zero") component and an
# inadmissible negative ("negative") estimate both count as not-interior here;
# the distinct wording for the negative case is carried by
# hs_fit_boundary_class() / the `at_boundary_condition` row. Returns NULL when
# the class is unavailable (row dropped), preserving the unavailable -> NULL
# contract.
#
# A fit is at the boundary when ANY variance component it reports is at/near zero
# relative to the variance total: the primary genetic / effect component
# (h2 -> 0), the residual (sigma_e2 -> 0, i.e. h2 -> 1), or a second effect such
# as permanent, common-environment, or maternal (-> 0). The primary genetic /
# effect component is named differently across targets ("animal" for the
# pedigree animal model, "genomic" for genomic and SNP-BLUP fits, "single_step"
# for single-step fits), so the check first confirms exactly one primary
# component is present. This restricts the flag to the univariate,
# repeatability, two-effect, genomic, and single-step layouts (each has one
# primary component) and leaves multivariate fits unflagged (NULL), since they
# report per-trait "genetic"/"residual" diagonals rather than a single primary
# share. This is the surfacing half of the v0.1 promotion predicate item 4; the
# engine (HSquared.jl) owns boundary-stable optimization.
hs_fit_boundary_flag <- function(object, tol = 1e-4) {
  cls <- hs_fit_boundary_class(object, tol = tol)
  if (is.null(cls)) {
    return(NULL)
  }
  cls %in% c("zero", "negative")
}

# Map a boundary class to the `at_boundary_condition` diagnostics value. Only the
# at/near-zero and inadmissible-negative cases get a row (NULL is dropped by the
# caller), so the condition row appears exactly when it carries information the
# logical `at_boundary` flag cannot: which kind of boundary the fit hit.
hs_fit_boundary_condition_label <- function(cls) {
  if (is.null(cls) || isFALSE(cls)) {
    return(NULL)
  }
  switch(
    cls,
    zero = "at/near zero (clean boundary)",
    negative = "negative (inadmissible variance)",
    NULL
  )
}

# Logical flag for the `search_boundary` diagnostics row (hsquared#230).
# Unrelated to `at_boundary`/`hs_fit_boundary_class()`: `at_boundary` is
# derived from the ESTIMATED variance shares, while `search_boundary` mirrors
# the non-Gaussian bridge's `result$boundary` field (`hs_ng09_boundary()` in
# R/julia-bridge.R), which reports whether HSquared.jl's single-variance Brent
# search stopped on the rail of its log-scale bracket
# (log(initial$sigma_a2) +/- 6, HSquared.jl#327) -- an optimizer-search state,
# not a property of the estimates themselves. Gaussian animal-model fits
# (`engine = "fit"`) never carry a `boundary` field, so this returns NULL for
# them and both `search_boundary` rows are dropped by the same NULL filter
# that drops every other absent diagnostic (no NA row either).
hs_fit_search_boundary_flag <- function(object) {
  boundary <- object$result$boundary
  if (!is.logical(boundary) || length(boundary) != 1L || is.na(boundary)) {
    return(NULL)
  }
  boundary
}

# Map the `search_boundary` flag to the `search_boundary_condition` value.
# Unlike `at_boundary_condition` (only emitted for the at/near-zero and
# negative cases, dropped otherwise), this row is always emitted alongside a
# present `search_boundary` flag: the search-bracket rail is a binary
# optimizer state with no benign-vs-inadmissible split to only-sometimes
# report, so the FALSE case reads "interior" rather than disappearing.
hs_fit_search_boundary_condition_label <- function(flag) {
  if (is.null(flag)) {
    return(NULL)
  }
  if (isTRUE(flag)) {
    return(paste0(
      "single-variance search stopped on the rail of the log-scale bracket ",
      "log(initial$sigma_a2) +/- 6, or (with restart_check = TRUE) the two ",
      "starts disagreed; either way the estimate is a function of the search, ",
      "not the data (HSquared.jl#327)"
    ))
  }
  "interior"
}

hs_diagnostic_value <- function(x) {
  if (length(x) == 0L) {
    return(NA_character_)
  }
  if (is.logical(x)) {
    return(paste(ifelse(x, "TRUE", "FALSE"), collapse = ", "))
  }
  if (is.numeric(x)) {
    return(paste(format(x, digits = 8, trim = TRUE), collapse = ", "))
  }
  if (is.character(x)) {
    return(paste(x, collapse = ", "))
  }
  if (is.factor(x)) {
    return(paste(as.character(x), collapse = ", "))
  }
  if (is.atomic(x)) {
    return(paste(as.character(x), collapse = ", "))
  }
  "<list>"
}

#' Extract planned marker, QTL, GWAS, and eQTL results
#'
#' `r lifecycle::badge("experimental")`
#'
#' These extractor names cover genomic, QTL, GWAS, and eQTL fitted results.
#' They return values only when an `hsquared_fit` object contains the
#' corresponding result field. `marker_effects()` and
#' `marker_variance_explained()` are populated by the opt-in SNP-BLUP path
#' (`target = "snp_blup"`). The variance-explained table is a descriptive
#' fitted-marker share, computed as effect squared times centered marker
#' variance and normalized across markers; it is not a marker-scan p-value,
#' QTL statistic, or causal decomposition under linkage disequilibrium. The
#' `gwas_table()` and `lod_scores()` methods for an already-computed `hs_gwas`
#' object expose the current uncalibrated marker-scan table and marker-level
#' LOD scores. Fit-level QTL/GWAS/eQTL result tables and map-annotated /
#' calibrated scan outputs remain reserved for future results.
#'
#' @inheritParams variance_components
#'
#' @return The requested marker or scan result for `hsquared_fit` objects that
#'   contain the corresponding field.
#' @name marker_extractors
NULL

#' @rdname marker_extractors
#' @export
marker_effects <- function(object, ...) {
  UseMethod("marker_effects")
}

#' @export
marker_effects.default <- function(object, ...) {
  hs_marker_extractor_default("marker_effects")
}

#' @export
marker_effects.hsquared_fit <- function(object, ...) {
  hs_fit_result(object, "marker_effects", "marker effects")
}

#' @rdname marker_extractors
#' @export
marker_variance_explained <- function(object, ...) {
  UseMethod("marker_variance_explained")
}

#' @export
marker_variance_explained.default <- function(object, ...) {
  hs_marker_extractor_default("marker_variance_explained")
}

#' @export
marker_variance_explained.hsquared_fit <- function(object, ...) {
  if (!is.null(object$result$marker_variance_explained)) {
    return(object$result$marker_variance_explained)
  }
  if (
    identical(object$spec$target, "snp_blup") &&
      !is.null(object$result$marker_effects) &&
      !is.null(object$payload$markers)
  ) {
    marker_design <- if (!is.null(object$payload$Z)) {
      as.matrix(object$payload$Z %*% object$payload$markers)
    } else {
      object$payload$markers
    }
    allele_frequencies <- object$result$marker_allele_frequencies
    if (length(allele_frequencies) == 0L) {
      allele_frequencies <- NULL
    }
    return(hs_marker_variance_explained_from_snp_blup(
      effects = object$result$marker_effects$effect,
      markers = marker_design,
      marker_labels = object$result$marker_effects$marker,
      allele_frequencies = allele_frequencies
    ))
  }
  hs_fit_result(
    object,
    "marker_variance_explained",
    "marker variance explained"
  )
}

#' @rdname marker_extractors
#' @export
qtl_table <- function(object, ...) {
  UseMethod("qtl_table")
}

#' @export
qtl_table.default <- function(object, ...) {
  hs_marker_extractor_default("qtl_table")
}

#' @export
qtl_table.hsquared_fit <- function(object, ...) {
  hs_fit_result(object, "qtl_table", "QTL table")
}

#' @rdname marker_extractors
#' @export
gwas_table <- function(object, ...) {
  UseMethod("gwas_table")
}

#' @export
gwas_table.default <- function(object, ...) {
  hs_marker_extractor_default("gwas_table")
}

#' @export
gwas_table.hsquared_fit <- function(object, ...) {
  hs_fit_result(object, "gwas_table", "GWAS table")
}

#' @export
gwas_table.hs_gwas <- function(object, ...) {
  out <- as.data.frame(object)
  attr(out, "scan_method") <- attr(object, "scan_method")
  calibration <- attr(object, "calibration")
  if (!is.null(calibration)) {
    attr(out, "calibration") <- calibration
  }
  out
}

#' @rdname marker_extractors
#' @export
eqtl_table <- function(object, ...) {
  UseMethod("eqtl_table")
}

#' @export
eqtl_table.default <- function(object, ...) {
  hs_marker_extractor_default("eqtl_table")
}

#' @export
eqtl_table.hsquared_fit <- function(object, ...) {
  hs_fit_result(object, "eqtl_table", "eQTL table")
}

#' @rdname marker_extractors
#' @export
lod_scores <- function(object, ...) {
  UseMethod("lod_scores")
}

#' @export
lod_scores.default <- function(object, ...) {
  hs_marker_extractor_default("lod_scores")
}

#' @export
lod_scores.hsquared_fit <- function(object, ...) {
  hs_fit_result(object, "lod_scores", "LOD scores")
}

#' @export
lod_scores.hs_gwas <- function(object, ...) {
  if (!all(c("marker", "lod") %in% names(object))) {
    stop(
      "`lod_scores()` requires an `hs_gwas` object with `marker` and `lod` ",
      "columns.",
      call. = FALSE
    )
  }
  out <- data.frame(
    marker = object$marker,
    lod = object$lod,
    stringsAsFactors = FALSE
  )
  attr(out, "scan_method") <- attr(object, "scan_method")
  calibration <- attr(object, "calibration")
  if (!is.null(calibration)) {
    attr(out, "calibration") <- calibration
  }
  out
}

hs_marker_extractor_default <- function(name) {
  stop(
    "`",
    name,
    "()` requires an `hsquared_fit` object with marker/QTL/eQTL results. ",
    "`marker_effects()` and `marker_variance_explained()` are populated only ",
    "by opt-in SNP-BLUP fits. For a relatedness-corrected post-fit marker scan, ",
    "use `gwas(fit, markers)` (experimental; nominal / Bonferroni / BH ",
    "p-values, NOT genome-wide calibrated), then call `gwas_table(scan)` or ",
    "`lod_scores(scan)` on the returned `hs_gwas` object. Fit-level ",
    "`gwas_table()` plus `qtl_table()` / `eqtl_table()` remain reserved for ",
    "the planned map-annotated scan API and any future calibrated-threshold ",
    "activation.",
    call. = FALSE
  )
}

#' Extract fixed effects
#'
#' `r lifecycle::badge("experimental")`
#'
#' `fixef()` is part of the planned v0.1 fitted-object contract for
#' `hsquared_fit` objects.
#'
#' @inheritParams variance_components
#'
#' @return Fixed-effect estimates for `hsquared_fit` objects.
#' @export
fixef <- function(object, ...) {
  UseMethod("fixef")
}

#' @export
fixef.default <- function(object, ...) {
  stop(
    "`fixef()` requires an `hsquared_fit` object. The current package does ",
    "not provide fixed effects for this object.",
    call. = FALSE
  )
}

#' @export
fixef.hsquared_fit <- function(object, ...) {
  hs_fit_result(object, "fixed_effects", "fixed-effect estimates")
}

#' @export
coef.hsquared_fit <- function(object, ...) {
  fixef(object, ...)
}

#' Extract random effects
#'
#' `r lifecycle::badge("experimental")`
#'
#' `ranef()` is part of the planned v0.1 fitted-object contract for
#' `hsquared_fit` objects.
#'
#' @inheritParams variance_components
#'
#' @return Random-effect estimates for `hsquared_fit` objects.
#' @export
ranef <- function(object, ...) {
  UseMethod("ranef")
}

#' @export
ranef.default <- function(object, ...) {
  stop(
    "`ranef()` requires an `hsquared_fit` object. The current package does ",
    "not provide random effects for this object.",
    call. = FALSE
  )
}

#' @export
ranef.hsquared_fit <- function(object, ...) {
  hs_fit_result(object, "random_effects", "random-effect estimates")
}

#' @export
logLik.hsquared_fit <- function(object, ...) {
  if (identical(object$result$converged, FALSE)) {
    stop(
      "Log-likelihood is unavailable because this `hsquared_fit` object did ",
      "not converge.",
      call. = FALSE
    )
  }
  value <- hs_fit_result(object, "loglik", "log-likelihood")
  out <- value
  class(out) <- "logLik"
  attr(out, "df") <- object$result$df %||% NA_integer_
  attr(out, "nobs") <- object$result$nobs %||% length(object$payload$y)
  out
}

#' @export
AIC.hsquared_fit <- function(object, ..., k = 2) {
  stats::AIC(stats::logLik(object), ..., k = k)
}

#' Block unsupported likelihood-inference helpers
#'
#' These methods intentionally fail with explicit scope messages. Validated
#' `confint()` / `vcov()` / `profile()` / `anova()` surfaces are not
#' implemented. Experimental interval extractors
#' ([heritability_interval()], [variance_component_standard_errors()]) exist
#' when the engine returned those fields; they are not coverage-calibrated.
#' Point estimates remain [variance_components()], [heritability()], and
#' [fit_diagnostics()].
#'
#' @param object An `hsquared_fit` object.
#' @param fitted An `hsquared_fit` object for [stats::profile()].
#' @param ... Reserved for future arguments.
#' @param parm,level Included for compatibility with [stats::confint()].
#'
#' @return These functions always error.
#' @name inference_blocks
NULL

hs_interval_sibling_hint <- function(object) {
  if (hs_fit_is_genomic(object)) {
    return(paste0(
      "`heritability_interval()` is not available for genomic fits ",
      "(no scale-labelled, calibrated interval). ",
      "`variance_component_standard_errors()` is experimental and only ",
      "present when the engine returned it."
    ))
  }
  if (hs_fit_is_random_regression(object)) {
    return(paste0(
      "On random-regression fits use `rr_heritability()` for the h2(t) ",
      "trajectory; there is no scalar `heritability_interval()`."
    ))
  }
  has_hi <- !is.null(object$result$heritability_interval)
  has_se <- !is.null(object$result$variance_component_se)
  if (has_hi && has_se) {
    return(paste0(
      "Use the experimental `heritability_interval()` and ",
      "`variance_component_standard_errors()` extractors ",
      "(not coverage-calibrated)."
    ))
  }
  if (has_hi) {
    return(paste0(
      "Use the experimental `heritability_interval()` extractor ",
      "(not coverage-calibrated)."
    ))
  }
  if (has_se) {
    return(paste0(
      "Use the experimental `variance_component_standard_errors()` extractor ",
      "(not coverage-calibrated)."
    ))
  }
  paste0(
    "Experimental `heritability_interval()` and ",
    "`variance_component_standard_errors()` exist when the engine returned ",
    "those fields; they are not coverage-calibrated."
  )
}

#' @rdname inference_blocks
#' @export
confint.hsquared_fit <- function(object, parm, level = 0.95, ...) {
  stop(
    "Validated confidence intervals for `hsquared_fit` quantities are not ",
    "implemented. ",
    hs_interval_sibling_hint(object),
    " Point estimates remain `variance_components()`, `heritability()`, and ",
    "`fit_diagnostics()`.",
    call. = FALSE
  )
}

#' @rdname inference_blocks
#' @export
vcov.hsquared_fit <- function(object, ...) {
  stop(
    "A validated estimator variance-covariance matrix is not implemented ",
    "for `hsquared_fit` objects. ",
    hs_interval_sibling_hint(object),
    call. = FALSE
  )
}

#' @rdname inference_blocks
#' @export
profile.hsquared_fit <- function(fitted, ...) {
  stop(
    "Profile-likelihood intervals for `hsquared_fit` objects are planned, not ",
    "implemented. v0.1 reports point estimates and convergence diagnostics ",
    "only.",
    call. = FALSE
  )
}

#' @rdname inference_blocks
#' @export
anova.hsquared_fit <- function(object, ...) {
  stop(
    "Likelihood-ratio / ANOVA comparison for `hsquared_fit` objects is planned, ",
    "not implemented. v0.1 exposes `logLik()` and `AIC()` for converged fits, ",
    "but does not yet validate LRT guidance.",
    call. = FALSE
  )
}

#' @importFrom stats nobs
#' @export
nobs.hsquared_fit <- function(object, ...) {
  n <- object$result$nobs %||%
    if (!is.null(object$payload$y)) length(object$payload$y) else NULL
  if (is.null(n)) {
    stop(
      "This `hsquared_fit` object does not contain number-of-observations ",
      "metadata.",
      call. = FALSE
    )
  }
  as.integer(n)
}

# Detect the opt-in multivariate target, whose fitted object stores `Y`
# (a multi-trait matrix) and emits no single-vector predictions. The response-
# scale `predict()`/`fitted()`/`residuals()` contract is univariate-only in
# v0.1, so these block on that target with a target-named scope message rather
# than the generic missing-field miss from `hs_fit_result()`.
hs_fit_is_multivariate <- function(object) {
  identical(object$spec$target, "multivariate") || !is.null(object$payload$Y)
}

hs_block_multivariate_response_scale <- function(name) {
  stop(
    "`",
    name,
    "()` on the response scale is univariate-only in v0.1 and is not defined ",
    "for the multivariate model (a `cbind()` response), which ",
    "fits multiple traits jointly. Use `breeding_values()`, ",
    "`genetic_covariance()`, and `residual_covariance()` for multivariate ",
    "results.",
    call. = FALSE
  )
}

#' Response-scale prediction helpers
#'
#' `predict()`, `fitted()`, and `residuals()` are part of the planned v0.1
#' fitted-object contract for univariate `hsquared_fit` objects. They are
#' univariate-only: the multivariate model (a `cbind()` response)
#' fits multiple traits jointly and is intentionally out of v0.1 response-scale
#' scope, so these methods stop with a scope message pointing to
#' `breeding_values()`, `genetic_covariance()`, and `residual_covariance()`.
#'
#' @inheritParams variance_components
#'
#' @return Response-scale predictions, fitted values, or residuals for
#'   univariate `hsquared_fit` objects.
#' @name response_scale_methods
NULL

#' @rdname response_scale_methods
#' @export
predict.hsquared_fit <- function(object, ...) {
  if (hs_fit_is_multivariate(object)) {
    hs_block_multivariate_response_scale("predict")
  }
  hs_fit_result(object, "predictions", "predictions")
}

#' @rdname response_scale_methods
#' @export
fitted.hsquared_fit <- function(object, ...) {
  if (hs_fit_is_multivariate(object)) {
    hs_block_multivariate_response_scale("fitted")
  }
  predictions <- stats::predict(object, ...)
  if (is.data.frame(predictions) && ".fitted" %in% names(predictions)) {
    return(predictions$.fitted)
  }
  predictions
}

#' @rdname response_scale_methods
#' @export
residuals.hsquared_fit <- function(object, ...) {
  if (hs_fit_is_multivariate(object)) {
    hs_block_multivariate_response_scale("residuals")
  }
  response <- object$payload$y
  if (is.null(response)) {
    stop(
      "This `hsquared_fit` object does not contain response values.",
      call. = FALSE
    )
  }
  fitted_values <- as.numeric(stats::fitted(object, ...))
  response <- as.numeric(response)
  if (length(response) != length(fitted_values)) {
    stop(
      "Response and fitted values must have the same length to compute ",
      "residuals.",
      call. = FALSE
    )
  }
  response - fitted_values
}

# ---------------------------------------------------------------------------
# Random-regression (reaction-norm) extractors
# ---------------------------------------------------------------------------

# Detect the opt-in random-regression target and fetch its result block. The
# trajectory extractors compute v_g(t)/h^2(t)/correlation in R from the supplied
# K_g and the recorded standardization bounds, so they need no live Julia.
hs_fit_is_random_regression <- function(object) {
  identical(object$spec$target, "random_regression") ||
    !is.null(object$result$coefficient_covariance)
}

hs_require_random_regression <- function(object, name) {
  if (!inherits(object, "hsquared_fit")) {
    stop(
      "`",
      name,
      "()` requires an `hsquared_fit` object from the opt-in random-regression ",
      "model (`target = \"random_regression\"`).",
      call. = FALSE
    )
  }
  if (!hs_fit_is_random_regression(object)) {
    stop(
      "`",
      name,
      "()` requires a fit from the opt-in random-regression model ",
      "(`target = \"random_regression\"`), fitted with an ",
      "`animal(rr(covariate, order = k) | id, pedigree = ped)` term.",
      call. = FALSE
    )
  }
  invisible(object)
}

# Resolve the covariate points the trajectory is evaluated at. `at = NULL`
# (the default) builds an evenly spaced grid of `n` points spanning the recorded
# covariate range; otherwise `at` is taken on the ORIGINAL covariate scale and
# validated against that range (#213 -- the Julia engine refuses a
# standardized t outside [-1, 1], so no R extractor may silently clamp one
# either). Returns the original-scale points and the matching standardized
# t in [-1, 1].
hs_rr_eval_points <- function(object, at, n = 25L) {
  rr <- object$result$random_regression
  if (is.null(rr)) {
    stop(
      "This random-regression `hsquared_fit` is missing its covariate ",
      "standardization metadata.",
      call. = FALSE
    )
  }
  if (is.null(at)) {
    at <- seq(rr$lower, rr$upper, length.out = n)
    t_std <- hs_standardize_covariate(at, rr$lower, rr$upper)
    return(list(at = at, t = t_std, covariate = rr$covariate, order = rr$order))
  }
  at <- as.numeric(at)
  if (length(at) == 0L || anyNA(at) || any(!is.finite(at))) {
    stop("`at` must be finite covariate values.", call. = FALSE)
  }
  t_std <- hs_standardize_covariate(at, rr$lower, rr$upper)
  out_of_range <- t_std < -1 - 1e-10 | t_std > 1 + 1e-10
  if (any(out_of_range)) {
    hs_abort_out_of_range(
      "`at` must lie inside the fitted covariate range [",
      format(rr$lower, trim = TRUE),
      ", ",
      format(rr$upper, trim = TRUE),
      "] (covariate `",
      rr$covariate,
      "`); got ",
      paste(format(at[out_of_range], trim = TRUE), collapse = ", "),
      "."
    )
  }
  list(at = at, t = t_std, covariate = rr$covariate, order = rr$order)
}

# Per-point additive genetic variance v_g(t) = phi(t)' K_g phi(t), clamped to >= 0
# (mirrors `HSquared.rr_genetic_variance`).
hs_rr_variance_values <- function(K_g, t_std, order) {
  phi <- hs_legendre_design(t_std, order)
  vapply(
    seq_along(t_std),
    function(i) {
      max(0, drop(phi[i, , drop = FALSE] %*% K_g %*% t(phi[i, , drop = FALSE])))
    },
    numeric(1L)
  )
}

#' Random-regression (reaction-norm) extractors
#'
#' `r lifecycle::badge("experimental")`
#'
#' These extractors summarize an opt-in, **experimental** random-regression
#' (reaction-norm) fit (`target = "random_regression"`), fitted with an
#' `animal(rr(covariate, order = k) | id, pedigree = ped)` term. The model
#' estimates a `k x k` genetic covariance matrix `K_g` among an animal's
#' normalized-Legendre random-regression coefficients plus a single homogeneous
#' residual variance.
#'
#' * `rr_covariance()` returns the estimated `k x k` coefficient genetic
#'   covariance matrix `K_g`.
#' * `random_coefficients()` returns the per-animal predicted Legendre
#'   coefficients (long format: `id`, `coefficient`, `value`).
#' * `rr_genetic_variance()` returns the additive genetic variance trajectory
#'   `v_g(t) = phi(t)' K_g phi(t)` across covariate points.
#' * `rr_heritability()` returns the heritability trajectory
#'   `h^2(t) = v_g(t) / (v_g(t) + sigma_e^2)`. Because the residual is
#'   homogeneous and there is no permanent-environment term yet, this can
#'   OVERSTATE `h^2(t)` for repeated-records designs (test-day, growth curves).
#' * `rr_correlation()` returns the genetic correlation surface among the
#'   covariate points.
#' * `rr_eigenfunctions()` returns the eigen-decomposition of the coefficient
#'   covariance `K_g` evaluated as functions of the covariate: the eigenvalues
#'   (additive genetic variance carried by each axis), the proportion of genetic
#'   variance explained, the sign-canonicalized eigen-coefficients, and the
#'   eigenfunctions `psi_j(t) = phi(t)' v_j` over the covariate grid. These are
#'   **rotation-invariant** functionals of `K_g` (mirroring the engine's
#'   `rr_eigenfunctions`); no raw, rotation-arbitrary loadings are returned.
#'
#' The trajectories are computed in R from the estimated `K_g` and the recorded
#' covariate standardization range; `at` is supplied on the ORIGINAL covariate
#' scale (defaulting to a grid over the fitted range) and re-standardized to
#' `[-1, 1]` internally, matching the Julia engine's basis convention.
#'
#' @section Out-of-range `at`:
#' `at` must lie inside the fitted covariate range (`object$result$random_regression$lower`
#' to `$upper`); a value outside that range errors with a message naming the
#' fitted range, instead of silently returning the nearest endpoint's value.
#' This matches the Julia engine (`HSquared.legendre_basis`), which refuses a
#' standardized covariate outside `[-1, 1]` rather than clamping it (#213).
#'
#' @param object A random-regression `hsquared_fit` object.
#' @param at Covariate values on the original scale at which to evaluate the
#'   trajectory. `NULL` (the default) uses an evenly spaced grid over the fitted
#'   covariate range. A supplied value outside the fitted covariate range
#'   errors (see the "Out-of-range `at`" section below).
#' @param n Number of grid points used when `at = NULL`.
#' @param ... Reserved for future arguments.
#'
#' @return `rr_covariance()` returns a numeric matrix; `random_coefficients()`
#'   returns a data frame; the trajectory extractors return a data frame with a
#'   `covariate` column and the evaluated `value`s.
#'
#' @examplesIf FALSE
#' fit_rr <- hsquared(
#'   weight ~ sex + animal(rr(age, order = 2) | id, pedigree = ped),
#'   data = long_records,
#'   family = gaussian(),
#'   REML = TRUE,
#'   control = hs_control(
#'     engine = "julia",
#'     engine_control = list(target = "random_regression")
#'   )
#' )
#'
#' rr_covariance(fit_rr)
#' random_coefficients(fit_rr)
#' rr_genetic_variance(fit_rr)
#' rr_heritability(fit_rr)
#' rr_correlation(fit_rr, at = c(1, 3, 5))
#' rr_eigenfunctions(fit_rr)
#' @name random_regression_extractors
NULL

#' @rdname random_regression_extractors
#' @export
rr_covariance <- function(object, ...) {
  UseMethod("rr_covariance")
}

#' @export
rr_covariance.default <- function(object, ...) {
  hs_require_random_regression(object, "rr_covariance")
}

#' @export
rr_covariance.hsquared_fit <- function(object, ...) {
  hs_require_random_regression(object, "rr_covariance")
  hs_fit_result(
    object,
    "coefficient_covariance",
    "the random-regression coefficient genetic covariance matrix"
  )
}

#' @rdname random_regression_extractors
#' @export
random_coefficients <- function(object, ...) {
  UseMethod("random_coefficients")
}

#' @export
random_coefficients.default <- function(object, ...) {
  hs_require_random_regression(object, "random_coefficients")
}

#' @export
random_coefficients.hsquared_fit <- function(object, ...) {
  hs_require_random_regression(object, "random_coefficients")
  hs_fit_result(
    object,
    "random_coefficients",
    "the predicted random-regression coefficients"
  )
}

#' @rdname random_regression_extractors
#' @export
rr_genetic_variance <- function(object, at = NULL, n = 25L, ...) {
  UseMethod("rr_genetic_variance")
}

#' @export
rr_genetic_variance.default <- function(object, at = NULL, n = 25L, ...) {
  hs_require_random_regression(object, "rr_genetic_variance")
}

#' @export
rr_genetic_variance.hsquared_fit <- function(object, at = NULL, n = 25L, ...) {
  hs_require_random_regression(object, "rr_genetic_variance")
  K_g <- rr_covariance(object)
  pts <- hs_rr_eval_points(object, at, n)
  data.frame(
    covariate = pts$at,
    value = hs_rr_variance_values(K_g, pts$t, pts$order),
    stringsAsFactors = FALSE
  )
}

#' @rdname random_regression_extractors
#' @export
rr_heritability <- function(object, at = NULL, n = 25L, ...) {
  UseMethod("rr_heritability")
}

#' @export
rr_heritability.default <- function(object, at = NULL, n = 25L, ...) {
  hs_require_random_regression(object, "rr_heritability")
}

#' @export
rr_heritability.hsquared_fit <- function(object, at = NULL, n = 25L, ...) {
  hs_require_random_regression(object, "rr_heritability")
  K_g <- rr_covariance(object)
  sigma_e2 <- hs_fit_result(
    object,
    "residual_variance",
    "the random-regression residual variance"
  )
  pts <- hs_rr_eval_points(object, at, n)
  vg <- hs_rr_variance_values(K_g, pts$t, pts$order)
  data.frame(
    covariate = pts$at,
    value = vg / (vg + as.numeric(sigma_e2)),
    stringsAsFactors = FALSE
  )
}

#' @rdname random_regression_extractors
#' @export
rr_correlation <- function(object, at = NULL, n = 25L, ...) {
  UseMethod("rr_correlation")
}

#' @export
rr_correlation.default <- function(object, at = NULL, n = 25L, ...) {
  hs_require_random_regression(object, "rr_correlation")
}

#' @export
rr_correlation.hsquared_fit <- function(object, at = NULL, n = 25L, ...) {
  hs_require_random_regression(object, "rr_correlation")
  K_g <- rr_covariance(object)
  pts <- hs_rr_eval_points(object, at, n)
  phi <- hs_legendre_design(pts$t, pts$order)
  G <- phi %*% K_g %*% t(phi)
  G <- 0.5 * (G + t(G))
  d <- diag(G)
  if (any(d <= 0)) {
    stop(
      "`rr_correlation()` is undefined: at least one covariate point has ",
      "non-positive genetic variance. Evaluate at points inside the genetic ",
      "variance support, or inspect `rr_genetic_variance()`.",
      call. = FALSE
    )
  }
  corr <- stats::cov2cor(G)
  dimnames(corr) <- list(
    format(pts$at, trim = TRUE),
    format(pts$at, trim = TRUE)
  )
  corr
}

#' @rdname random_regression_extractors
#' @export
rr_eigenfunctions <- function(object, at = NULL, n = 25L, ...) {
  UseMethod("rr_eigenfunctions")
}

#' @export
rr_eigenfunctions.default <- function(object, at = NULL, n = 25L, ...) {
  hs_require_random_regression(object, "rr_eigenfunctions")
}

#' @export
rr_eigenfunctions.hsquared_fit <- function(object, at = NULL, n = 25L, ...) {
  hs_require_random_regression(object, "rr_eigenfunctions")
  K_g <- rr_covariance(object)
  e <- hs_g_eigen(K_g) # eigenvalues (desc) + sign-canonicalized eigenvectors
  pts <- hs_rr_eval_points(object, at, n)
  phi <- hs_legendre_design(pts$t, pts$order)
  psi <- phi %*% e$vectors # eigenfunctions over the covariate grid
  total <- sum(e$values)
  prop <- if (total > 0) e$values / total else rep(0, length(e$values))
  k <- length(e$values)
  list(
    covariate = pts$at,
    eigenvalues = e$values,
    variance_explained = prop,
    eigen_coefficients = e$vectors,
    eigenfunctions = data.frame(
      covariate = rep(pts$at, times = k),
      axis = rep(seq_len(k), each = length(pts$at)),
      value = as.numeric(psi),
      stringsAsFactors = FALSE
    )
  )
}

# ---------------------------------------------------------------------------
# Opt-in direct-maternal correlated model extractors (Phase 4)
# ---------------------------------------------------------------------------

# Detect a direct-maternal fit. Checks both the stored spec target and a
# result field that only this model populates.
hs_fit_is_direct_maternal <- function(object) {
  identical(object$spec$target, "direct_maternal") ||
    (!is.null(object$result$direct_variance) &&
      !is.null(object$result$partner_variance))
}

hs_require_direct_maternal <- function(object, name) {
  if (!inherits(object, "hsquared_fit")) {
    stop(
      "`",
      name,
      "()` requires an `hsquared_fit` object from the opt-in ",
      "direct-maternal correlated model (`target = \"direct_maternal\"`).",
      call. = FALSE
    )
  }
  if (!hs_fit_is_direct_maternal(object)) {
    stop(
      "`",
      name,
      "()` requires a fit from the opt-in direct-maternal ",
      "correlated model (`target = \"direct_maternal\"`), fitted with ",
      "`animal(1 | id, pedigree = ped) + maternal_genetic(1 | dam)`.",
      call. = FALSE
    )
  }
  invisible(object)
}

#' Direct-maternal correlated model extractors
#'
#' `r lifecycle::badge("experimental")`
#'
#' These extractors summarize an opt-in, **covered at validation scale**
#' direct-maternal correlated model (`target = "direct_maternal"`), fitted with
#' `animal(1 | id, pedigree = ped) + maternal_genetic(1 | dam)`.
#' The model estimates a 2x2 genetic covariance matrix G_dm between the direct
#' additive and maternal additive effects, plus a residual variance.
#'
#' Willham fence: the direct-maternal model (Willham 1963, 1972) distinguishes
#' the **direct** narrow-sense heritability `h2_d = sigma_ad / sigma^2_P` from
#' the maternal variance ratio `m2 = sigma_am / sigma^2_P`, the Willham total
#' (selection-response) heritability
#' `h2_T = (sigma_ad + 1.5*sigma_dm + 0.5*sigma_am) / sigma^2_P`, and the
#' genetic covariance `sigma_dm`.
#' `sigma^2_P = sigma_ad + sigma_am + sigma_dm + sigma_e2 = Var(y_i)` for a
#' non-inbred base (coefficient 1 on `sigma_dm` because
#' `2 * A[i,dam] = 2 * (1/2) = 1`). **h2 is denominator-dependent under
#' maternal effects; compare (co)variance components, not h2 values, across
#' software (ASReml/BLUPF90/WOMBAT/sommer/MCMCglmm all leave sigma^2_P to the
#' user).** A **negative** genetic correlation `r_am` is real and biologically
#' expected in many livestock traits; it reflects an antagonistic
#' direct-maternal relationship and does NOT indicate a model failure.
#'
#' Note: the 2x2 G_dm formulation is due to Willham (1963, 1972), not
#' Falconer (1965). Falconer's single-m model fixes `r_am = +/-1` and is
#' a special case (Bijma 2011).
#'
#' `direct_heritability()` returns the direct narrow-sense heritability with
#' an interpretation attribute. `total_heritability()` returns Willham's
#' selection-response heritability `h2_T`. `direct_variance()`,
#' `partner_variance()`, and `direct_maternal_covariance()` return the raw
#' variance-component scalars. `genetic_correlation()` on a direct-maternal
#' fit returns the between-effect correlation r_am.
#'
#' @param object An `hsquared_fit` from `target = "direct_maternal"`.
#' @param ... Not used.
#'
#' @return
#' - `direct_heritability()`: a one-row data frame with `term` and `estimate`,
#'   plus an `"interpretation"` attribute.
#' - `total_heritability()`: a one-row data frame with `term = "total_willham"`
#'   and `estimate`, plus an `"interpretation"` attribute.
#' - `direct_variance()`, `partner_variance()`, `direct_maternal_covariance()`:
#'   single numerics (the raw variance component).
#'
#' @seealso [variance_components()], [genetic_correlation()], [maternal_effects()],
#'   [total_heritability()]
#' @examplesIf FALSE
#' fit_dm <- hsquared(
#'   y ~ 1 + animal(1 | id, pedigree = ped) + maternal_genetic(1 | dam),
#'   data = dat,
#'   family = gaussian(),
#'   REML = TRUE,
#'   control = hs_control(
#'     engine = "julia",
#'     engine_control = list(target = "direct_maternal")
#'   )
#' )
#'
#' direct_heritability(fit_dm)    # h2_d (direct narrow-sense)
#' total_heritability(fit_dm)     # h2_T (Willham selection-response)
#' genetic_correlation(fit_dm)    # r_am (may be negative)
#' direct_variance(fit_dm)        # sigma_ad
#' partner_variance(fit_dm)       # sigma_am
#' direct_maternal_covariance(fit_dm)  # sigma_dm
#' variance_components(fit_dm)    # all four components
#' maternal_effects(fit_dm)       # maternal EBVs
#' @name direct_maternal_extractors
NULL

#' @rdname direct_maternal_extractors
#' @export
direct_heritability <- function(object, ...) {
  UseMethod("direct_heritability")
}

#' @export
direct_heritability.default <- function(object, ...) {
  stop(
    "`direct_heritability()` requires an `hsquared_fit` object from the opt-in ",
    "direct-maternal correlated model (`target = \"direct_maternal\"`).",
    call. = FALSE
  )
}

#' @export
direct_heritability.hsquared_fit <- function(object, ...) {
  hs_require_direct_maternal(object, "direct_heritability")
  out <- hs_fit_result(object, "heritability", "direct heritability estimate")
  attr(out, "interpretation") <- paste(
    "Direct narrow-sense heritability from the direct-maternal correlated",
    "model: h2_d = sigma_ad / sigma^2_P, where",
    "sigma^2_P = sigma_ad + sigma_am + sigma_dm + sigma_e2 = Var(y_i)",
    "(Willham 1963, 1972; coefficient 1 on sigma_dm because",
    "2*A[i,dam] = 2*(1/2) = 1 for a non-inbred base).",
    "h2 is denominator-dependent under maternal effects; hsquared defines",
    "sigma^2_P to include sigma_dm (Willham 1972). Compare (co)variance",
    "components rather than h2 values across software.",
    "A negative genetic correlation r_am is real and expected in many",
    "livestock traits. Use total_heritability() for Willham's selection-",
    "response h2_T, genetic_correlation() for r_am, and",
    "partner_variance() for sigma_am."
  )
  out
}

#' @rdname direct_maternal_extractors
#' @export
direct_variance <- function(object, ...) {
  UseMethod("direct_variance")
}

#' @export
direct_variance.default <- function(object, ...) {
  stop(
    "`direct_variance()` requires an `hsquared_fit` from the opt-in ",
    "direct-maternal correlated model (`target = \"direct_maternal\"`).",
    call. = FALSE
  )
}

#' @export
direct_variance.hsquared_fit <- function(object, ...) {
  hs_require_direct_maternal(object, "direct_variance")
  as.numeric(object$result$direct_variance)
}

#' @rdname direct_maternal_extractors
#' @export
partner_variance <- function(object, ...) {
  UseMethod("partner_variance")
}

#' @export
partner_variance.default <- function(object, ...) {
  stop(
    "`partner_variance()` requires an `hsquared_fit` from the opt-in ",
    "direct-maternal correlated model (`target = \"direct_maternal\"`).",
    call. = FALSE
  )
}

#' @export
partner_variance.hsquared_fit <- function(object, ...) {
  hs_require_direct_maternal(object, "partner_variance")
  as.numeric(object$result$partner_variance)
}

#' @rdname direct_maternal_extractors
#' @export
direct_maternal_covariance <- function(object, ...) {
  UseMethod("direct_maternal_covariance")
}

#' @export
direct_maternal_covariance.default <- function(object, ...) {
  stop(
    "`direct_maternal_covariance()` requires an `hsquared_fit` from the opt-in ",
    "direct-maternal correlated model (`target = \"direct_maternal\"`).",
    call. = FALSE
  )
}

#' @export
direct_maternal_covariance.hsquared_fit <- function(object, ...) {
  hs_require_direct_maternal(object, "direct_maternal_covariance")
  as.numeric(object$result$covariance)
}

#' @rdname direct_maternal_extractors
#' @export
total_heritability <- function(object, ...) {
  UseMethod("total_heritability")
}

#' @export
total_heritability.default <- function(object, ...) {
  stop(
    "`total_heritability()` requires an `hsquared_fit` object from the opt-in ",
    "direct-maternal correlated model (`target = \"direct_maternal\"`).",
    call. = FALSE
  )
}

#' @export
total_heritability.hsquared_fit <- function(object, ...) {
  hs_require_direct_maternal(object, "total_heritability")
  vc <- object$result$variance_components
  sigma_ad <- as.numeric(vc$estimate[vc$component == "direct"])
  sigma_am <- as.numeric(vc$estimate[vc$component == "maternal"])
  sigma_dm <- as.numeric(vc$estimate[vc$component == "covariance"])
  sigma_e2 <- as.numeric(vc$estimate[vc$component == "residual"])
  sigma_P <- sigma_ad + sigma_am + sigma_dm + sigma_e2
  h2_T <- if (sigma_P > 0) {
    (sigma_ad + 1.5 * sigma_dm + 0.5 * sigma_am) / sigma_P
  } else {
    NA_real_
  }
  out <- data.frame(
    term = "total_willham",
    estimate = h2_T,
    stringsAsFactors = FALSE
  )
  attr(out, "interpretation") <- paste(
    "Willham total (selection-response) heritability:",
    "h2_T = (sigma_ad + 1.5*sigma_dm + 0.5*sigma_am) / sigma^2_P",
    "(Willham 1963, 1972; coefficients 1, 1.5, 0.5 for direct,",
    "covariance, maternal terms). PREDICTS RESPONSE TO MASS SELECTION.",
    "h2_T CAN BE LOWER than direct h2_d when r_am < 0 (antagonistic",
    "direct-maternal covariance) - this is real and expected.",
    "sigma^2_P = sigma_ad + sigma_am + sigma_dm + sigma_e2 = Var(y_i).",
    "Do NOT confuse with 'total heritable variance' (coefficients 1, 2, 1)",
    "which is a different quantity. Use direct_heritability() for h2_d",
    "and genetic_correlation() for r_am."
  )
  out
}
