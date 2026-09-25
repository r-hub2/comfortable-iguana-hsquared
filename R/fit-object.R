hs_new_fit <- function(
  call = NULL,
  spec,
  payload,
  result,
  engine = "HSquared.jl",
  version = utils::packageVersion("hsquared")
) {
  if (!is.list(spec)) {
    stop("`spec` must be a list.", call. = FALSE)
  }
  if (!is.list(payload)) {
    stop("`payload` must be a list.", call. = FALSE)
  }
  if (!is.list(result)) {
    stop("`result` must be a list.", call. = FALSE)
  }

  structure(
    list(
      call = call,
      spec = spec,
      payload = payload,
      result = result,
      engine = engine,
      version = as.character(version)
    ),
    class = "hsquared_fit"
  )
}

hs_fit_target_label <- function(object) {
  target <- object$spec$target
  if (is.null(target) || !nzchar(as.character(target)[[1L]])) {
    return("default animal")
  }
  as.character(target)[[1L]]
}

hs_fit_formula_label <- function(object) {
  if (inherits(object$spec$formula, "formula")) {
    return(paste(deparse(object$spec$formula), collapse = " "))
  }
  if (is.language(object$call) && length(object$call) >= 2L) {
    return(paste(deparse(object$call[[2L]]), collapse = " "))
  }
  NULL
}

hs_fit_result_sibling <- function(name) {
  switch(
    name,
    qtl_table = ,
    gwas_table = ,
    eqtl_table = ,
    lod_scores = paste0(
      " Use `gwas(fit, markers)` for a marker scan, then `gwas_table(scan)` ",
      "or `lod_scores(scan)` on the returned `hs_gwas` object."
    ),
    ""
  )
}

hs_fit_result <- function(object, name, label) {
  if (!inherits(object, "hsquared_fit")) {
    stop("`object` must be an `hsquared_fit` object.", call. = FALSE)
  }

  value <- object$result[[name]]
  if (is.null(value)) {
    stop(
      "This `hsquared_fit` object (target = \"",
      hs_fit_target_label(object),
      "\") does not contain ",
      label,
      ". The current result payload did not provide this field.",
      hs_fit_result_sibling(name),
      call. = FALSE
    )
  }

  value
}

hs_fit_not_converged <- function(object) {
  if (identical(object$result$converged, FALSE)) {
    return(TRUE)
  }
  identical(object$result$diagnostics$optimizer_status, "not_converged")
}

# Students copy README, get h2 ~ 0 from a failed n = 4 fit, and believe it.
# logLik() refuses a non-converged fit; heritability() and print() must
# warn at the same bar so a near-zero value is never a silent "result".
hs_warn_if_unusable_fit <- function(object, what = "heritability") {
  if (!hs_fit_not_converged(object)) {
    return(invisible(FALSE))
  }
  boundary <- isTRUE(hs_fit_boundary_flag(object))
  extra <- if (boundary) {
    paste0(
      " A variance component is also at or near a boundary, so a ",
      "near-zero value is a failed-fit artefact, not evidence that ",
      "heritability is zero."
    )
  } else {
    " A near-zero value is not evidence that heritability is zero."
  }
  warning(
    "This `hsquared_fit` object did not converge. The ",
    what,
    " number is not an estimate; do not report it.",
    extra,
    " Inspect `fit_diagnostics(fit)` before reading any number.",
    call. = FALSE
  )
  invisible(TRUE)
}

hs_print_fit_peek <- function(x) {
  h2 <- x$result$heritability
  if (is.null(h2)) {
    return(invisible(x))
  }
  if (is.data.frame(h2) && "estimate" %in% names(h2)) {
    labels <- if ("term" %in% names(h2)) {
      h2$term
    } else if ("component" %in% names(h2)) {
      h2$component
    } else {
      seq_len(nrow(h2))
    }
    peek <- paste0(
      labels,
      "=",
      format(signif(as.numeric(h2$estimate), 4)),
      collapse = ", "
    )
    cat("  heritability: ", peek, "\n", sep = "")
  } else if (is.numeric(h2) && length(h2) == 1L) {
    cat("  heritability: ", format(signif(as.numeric(h2), 4)), "\n", sep = "")
  }
  invisible(x)
}

#' @export
print.hsquared_fit <- function(x, ...) {
  method <- x$spec$method %||% "unknown"
  family <- x$spec$family$family %||% "unknown"
  converged <- x$result$converged

  cat("<hsquared_fit>\n")
  formula_txt <- hs_fit_formula_label(x)
  if (!is.null(formula_txt)) {
    cat("  formula: ", formula_txt, "\n", sep = "")
  }
  cat("  target: ", hs_fit_target_label(x), "\n", sep = "")
  cat("  engine: ", x$engine %||% "unknown", "\n", sep = "")
  cat("  family: ", family, "\n", sep = "")
  cat("  method: ", method, "\n", sep = "")
  if (!is.null(converged)) {
    cat("  converged: ", isTRUE(converged), "\n", sep = "")
  }
  if (hs_fit_not_converged(x)) {
    cat("  heritability: not reportable (fit did not converge)\n")
    hs_warn_if_unusable_fit(x)
  } else {
    hs_print_fit_peek(x)
  }
  boundary <- x$result$genomic_boundary
  if (!is.null(boundary)) {
    cat("  genomic boundary status: ", boundary$status, "\n", sep = "")
    if (boundary$status %in% c("boundary_lower", "boundary_upper")) {
      cat(
        "  genomic ratio (scientific endpoint): ",
        boundary$profile_ratio,
        "\n",
        sep = ""
      )
      cat(
        "  genomic ratio (numerical MME): ",
        boundary$numerical_ratio,
        "\n",
        sep = ""
      )
    }
  }
  invisible(x)
}

#' @export
summary.hsquared_fit <- function(object, ...) {
  structure(
    list(
      call = object$call,
      engine = object$engine,
      method = object$spec$method,
      family = object$spec$family,
      variance_components = object$result$variance_components,
      heritability = object$result$heritability,
      fixed_effects = object$result$fixed_effects,
      diagnostics = object$result$diagnostics,
      genomic_boundary = object$result$genomic_boundary,
      converged = object$result$converged,
      at_boundary = hs_fit_boundary_flag(object),
      at_boundary_class = hs_fit_boundary_class(object),
      # Experimental uncertainty surfaces (engine rows V1-HERIT-CI /
      # V3-REPEAT-REML, partial); present only when the engine returned them.
      heritability_interval = object$result$heritability_interval,
      heritability_se = object$result$heritability_se,
      variance_component_se = object$result$variance_component_se,
      repeatability_interval = object$result$repeatability_interval
    ),
    class = "summary_hsquared_fit"
  )
}

#' @export
print.summary_hsquared_fit <- function(x, ...) {
  cat("<summary_hsquared_fit>\n")
  cat("  engine: ", x$engine %||% "unknown", "\n", sep = "")
  cat("  method: ", x$method %||% "unknown", "\n", sep = "")
  if (!is.null(x$converged)) {
    cat("  converged: ", isTRUE(x$converged), "\n", sep = "")
  }
  if (!is.null(x$genomic_boundary)) {
    cat(
      "  genomic boundary status: ",
      x$genomic_boundary$status,
      "\n",
      sep = ""
    )
    if (x$genomic_boundary$status %in% c("boundary_lower", "boundary_upper")) {
      cat(
        "  genomic ratio (scientific endpoint): ",
        x$genomic_boundary$profile_ratio,
        "\n",
        sep = ""
      )
      cat(
        "  genomic ratio (numerical MME): ",
        x$genomic_boundary$numerical_ratio,
        "\n",
        sep = ""
      )
    }
  }
  if (isTRUE(x$at_boundary)) {
    if (identical(x$at_boundary_class, "negative")) {
      cat(
        "  at boundary: TRUE (a variance component is NEGATIVE; the fit is ",
        "inadmissible, not a clean boundary, and heritability cannot be read ",
        "as an ordinary interior estimate)\n",
        sep = ""
      )
    } else {
      cat(
        "  at boundary: TRUE (a variance component is at/near zero; ",
        "heritability is a boundary estimate, not an ordinary interior one)\n",
        sep = ""
      )
    }
  }
  hs_print_uncertainty(x)
  invisible(x)
}

# Print the experimental uncertainty surfaces (CIs / SEs) when present.
# Governing evidence is HSquared.jl DRAC fir job 47925485 (2000-rep C1
# confirm; doc-34 sec. 4; unratified until Layer B land). Precursor jobs
# 46853279 / 47870067 do not govern. h^2 is directional-conservative,
# leg-specific: delta over-covers; profile in-band; bootstrap worst-cell
# C 0.924 (mild under vs 0.95, >= 0.90). Do not write "never under-covers."
# sigma^2_a profile / bootstrap are directional-conservative (not point).
# sigma^2_a delta/Wald under-covers (0.897) -> experimental-only; the old
# "delta_z ~0.92, approximately nominal" note is superseded and must not return.
# Repeatability is not in C1; its label stays experimental / no direction.
# Labelled experimental so they are never read as validated.
hs_print_uncertainty <- function(x) {
  fmt <- function(v) format(signif(as.numeric(v), 4))
  if (!is.null(x$heritability_se) || !is.null(x$heritability_interval)) {
    cat(
      "  heritability uncertainty (experimental; directional-conservative / not coverage-calibrated; asymptotic REML):\n"
    )
    if (!is.null(x$heritability_se)) {
      cat("    SE: ", fmt(x$heritability_se), "\n", sep = "")
    }
    if (!is.null(x$heritability_interval)) {
      ci <- x$heritability_interval
      cat(
        sprintf(
          "    %g%% CI: [%s, %s] (%s)\n",
          100 * as.numeric(ci$level),
          fmt(ci$lower),
          fmt(ci$upper),
          ci$method
        )
      )
    }
  }
  if (!is.null(x$variance_component_se)) {
    cat(
      "  variance-component SEs (experimental; NOT coverage-calibrated -- the sigma_a2 delta/Wald leg under-covers, ~0.897 at nominal 0.95; asymptotic REML; point estimate +/- SE for reference only, not a calibrated or conservative interval):\n"
    )
    se <- x$variance_component_se
    for (i in seq_len(nrow(se))) {
      cat("    ", se$component[i], ": ", fmt(se$se[i]), "\n", sep = "")
    }
  }
  if (!is.null(x$repeatability_interval)) {
    ci <- x$repeatability_interval
    cat("  repeatability uncertainty (experimental; asymptotic REML):\n")
    cat(
      sprintf(
        "    %g%% CI: [%s, %s]\n",
        100 * as.numeric(ci$level),
        fmt(ci$lower),
        fmt(ci$upper)
      )
    )
  }
  invisible(x)
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
