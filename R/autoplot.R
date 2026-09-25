#' ggplot2 visualizations for hsquared results
#'
#' `r lifecycle::badge("experimental")`
#'
#' `autoplot()` methods render the quantitative-genetic results an
#' `hsquared_fit` (or a `gwas()` scan) carries as `ggplot2` objects, in the
#' style of the `brms`/`bayesplot` ecosystem and consistent with the sister
#' packages `drmTMB`/`gllvmTMB`. They are **uncertainty-first**: where the fit
#' carries the experimental standard errors / reliabilities, the figures show
#' them (clearly labelled experimental and asymptotic).
#'
#' Available `type`s for `autoplot.hsquared_fit()`:
#'
#' * `"variance"` (default) -- a horizontal forest of the variance components
#'   and per-trait `h^2`, each with approximate 95% intervals (`+/- 1.96 * SE`)
#'   when the fit carries the experimental standard errors.
#' * `"breeding_values"` -- a sorted caterpillar of the estimated breeding
#'   values, with `+/- 1.96 * sqrt(PEV)` bands when prediction error variances
#'   are available, faceted by trait for multivariate fits.
#' * `"g_matrix"` -- a **rotation-invariant** genetic-correlation heatmap of the
#'   estimated `G` for multivariate fits (correlations are invariant to the
#'   factor rotation; raw loadings are never plotted -- the ratified cross-lane
#'   convention). Off-diagonal cells involving a low-`h^2` trait are flagged as
#'   imprecise (threshold `low_h2`, default `0.1`).
#' * `"g_geometry"` -- a scree of the **rotation-invariant** genetic
#'   eigenstructure (eigenvalues = variance per genetic axis, with percent
#'   variance explained) for multivariate fits. Axis directions / loadings are
#'   never drawn (rotation-arbitrary; span-ambiguous under repeated eigenvalues).
#' * `"reaction_norm"` -- for random-regression fits, the genetic-variance and
#'   heritability trajectories across the covariate (faceted). The heritability
#'   trajectory carries the same caveat as [rr_heritability()]: with a
#'   homogeneous residual and no permanent-environment term it can overstate
#'   `h^2(t)` for repeated-records designs.
#' * `"rr_eigenfunctions"` -- for random-regression fits, the rotation-invariant
#'   eigenfunctions `psi_j(t)` of `K_g` as covariate functions (faceted by axis,
#'   labelled by percent genetic variance). Signs are arbitrary and the curves are
#'   span-ambiguous under repeated eigenvalues (do not over-read).
#' * `"rr_surface"` -- for random-regression fits, the genetic covariance
#'   surface `S(s, t) = phi(s)' K_g phi(t)` over the covariate grid as a heatmap
#'   (pass `correlation = TRUE` for the genetic-correlation surface). Supplied-`K_g`
#'   descriptive.
#'
#' `autoplot()` on a `gwas()` scan (`hs_gwas`) draws `type = "manhattan"`
#' (default) or `type = "qq"` (observed vs expected `-log10(p)` with a `y = x`
#' null and the genomic-inflation `lambda_GC` as a diagnostic). Both carry the
#' EXPERIMENTAL, NOT-genome-wide-calibrated caveat. The Julia lane has a
#' fixed-panel calibration smoke harness, but R significance thresholds are not
#' activated.
#'
#' The figure helpers are deliberately modular (each takes a tidy data frame and
#' returns a `ggplot`) so they can be factored into a shared visualization
#' package later.
#'
#' @param object An `hsquared_fit` or `hs_gwas` object.
#' @param type Which figure to draw (see Details).
#' @param ... Figure-specific options passed through: `low_h2` (the
#'   genetic-correlation heatmap flags off-diagonal cells involving a trait with
#'   `h^2 < low_h2` as imprecise; default `0.1`), `at`/`n` (the reaction-norm
#'   covariate evaluation points), and `correlation` (`rr_surface`: draw the
#'   genetic-correlation surface instead of the covariance surface).
#'
#' @return A `ggplot` object.
#' @name hsquared-autoplot
NULL

#' @importFrom ggplot2 autoplot
#' @importFrom ggplot2 .data
#' @export
ggplot2::autoplot

#' hsquared ggplot2 theme
#'
#' `r lifecycle::badge("experimental")`
#'
#' A light, publication-oriented `ggplot2` theme shared by the hsquared
#' `autoplot()` figures. Exported so users can restyle or extend the figures.
#'
#' @param base_size Base font size.
#' @return A `ggplot2` theme object.
#' @export
theme_hsquared <- function(base_size = 12) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(colour = "grey35"),
      strip.text = ggplot2::element_text(face = "bold")
    )
}

# --- modular figure builders (tidy data frame -> ggplot) -------------------

# Horizontal forest: data frame with `term`, `estimate`, optional `lo`/`hi`,
# and an optional `panel` grouping for faceting.
hs_gg_forest <- function(
  df,
  xlab = "estimate",
  title = NULL,
  subtitle = NULL,
  zero_line = TRUE
) {
  df$term <- factor(df$term, levels = rev(unique(df$term)))
  has_ci <- all(c("lo", "hi") %in% names(df)) && any(is.finite(df$lo))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$estimate, y = .data$term))
  if (zero_line) {
    p <- p +
      ggplot2::geom_vline(
        xintercept = 0,
        linetype = 3,
        colour = "grey70"
      )
  }
  if (has_ci) {
    p <- p +
      ggplot2::geom_errorbar(
        ggplot2::aes(xmin = .data$lo, xmax = .data$hi),
        orientation = "y",
        width = 0.18,
        na.rm = TRUE,
        colour = "grey45"
      )
  }
  p <- p + ggplot2::geom_point(size = 2.6, colour = "#2c6fbb")
  if ("panel" %in% names(df)) {
    p <- p +
      ggplot2::facet_grid(
        rows = ggplot2::vars(.data$panel),
        scales = "free",
        space = "free_y"
      )
  }
  p +
    ggplot2::labs(x = xlab, y = NULL, title = title, subtitle = subtitle) +
    theme_hsquared()
}

# --- autoplot.hsquared_fit -------------------------------------------------

#' @rdname hsquared-autoplot
#' @exportS3Method ggplot2::autoplot
autoplot.hsquared_fit <- function(
  object,
  type = c(
    "variance",
    "breeding_values",
    "g_matrix",
    "g_geometry",
    "reaction_norm",
    "rr_eigenfunctions",
    "rr_surface"
  ),
  ...
) {
  hs_require_ggplot2()
  type <- match.arg(type)
  switch(
    type,
    variance = hs_autoplot_variance(object, ...),
    breeding_values = hs_autoplot_breeding_values(object, ...),
    g_matrix = hs_autoplot_g_matrix(object, ...),
    g_geometry = hs_autoplot_g_geometry(object, ...),
    reaction_norm = hs_autoplot_reaction_norm(object, ...),
    rr_eigenfunctions = hs_autoplot_rr_eigenfunctions(object, ...),
    rr_surface = hs_autoplot_rr_surface(object, ...)
  )
}

hs_autoplot_variance <- function(object, ...) {
  # Auto-detect the engine Set-B `variance_components_plot_data` payload; else
  # assemble from the stored extractors (recompute fallback).
  pd <- object$result$variance_components_plot_data
  if (!is.null(pd) && all(c("term", "estimate") %in% names(pd))) {
    # The payload is already shaped to the hs_gg_forest contract
    # (term/estimate/lo/hi/panel + interval_status). NaN -> NA so ggplot draws no
    # whisker where the interval is unavailable; lo/hi are RAW (unclamped).
    n <- length(pd$term)
    lo <- if (!is.null(pd$lo)) as.numeric(pd$lo) else rep(NA_real_, n)
    hi <- if (!is.null(pd$hi)) as.numeric(pd$hi) else rep(NA_real_, n)
    lo[is.nan(lo)] <- NA_real_
    hi[is.nan(hi)] <- NA_real_
    df <- data.frame(
      term = as.character(pd$term),
      estimate = as.numeric(pd$estimate),
      panel = if (!is.null(pd$panel)) {
        as.character(pd$panel)
      } else {
        "variance components"
      },
      lo = lo,
      hi = hi,
      stringsAsFactors = FALSE
    )
    # Binary by design: the v1 engine contract emits only "experimental_asymptotic"
    # or "none" (likelihood.jl); a future third status would be a known follow-up
    # to map explicitly, not a silent relabel.
    experimental <- !is.null(pd$interval_status) &&
      !identical(as.character(pd$interval_status), "none")
    # The [0,1] boundary annotation applies to h^2 only (a variance whisker may
    # cross zero, which is expected/honest, not a boundary crossing). NB the
    # engine's h^2 row is a logit-delta interval (always in (0,1)), so on real
    # engine output this never fires -- it is a defensive guard kept symmetric with
    # the recompute path, which uses raw natural-scale bounds that genuinely can.
    is_h2 <- df$panel == "heritability"
    boundary <- any(is_h2 & is.finite(df$lo) & (df$lo < 0 | df$hi > 1))
  } else {
    vc <- object$result$variance_components
    if (is.null(vc) || !all(c("component", "estimate") %in% names(vc))) {
      stop(
        "This `hsquared_fit` has no variance-component estimates to plot.",
        call. = FALSE
      )
    }
    se <- object$result$variance_component_se
    vc_df <- data.frame(
      term = as.character(vc$component),
      estimate = as.numeric(vc$estimate),
      panel = "variance components",
      stringsAsFactors = FALSE
    )
    vc_df$lo <- NA_real_
    vc_df$hi <- NA_real_
    experimental <- FALSE
    if (!is.null(se) && all(c("component", "se") %in% names(se))) {
      idx <- match(vc_df$term, as.character(se$component))
      sev <- as.numeric(se$se)[idx]
      vc_df$lo <- vc_df$estimate - 1.96 * sev
      vc_df$hi <- vc_df$estimate + 1.96 * sev
      experimental <- any(is.finite(sev))
    }

    # heritability panel (with SE if available)
    her <- object$result$heritability
    if (!is.null(her) && all(c("term", "estimate") %in% names(her))) {
      h_df <- data.frame(
        term = paste0("h\u00b2[", as.character(her$term), "]"),
        estimate = as.numeric(her$estimate),
        panel = "heritability",
        lo = NA_real_,
        hi = NA_real_,
        stringsAsFactors = FALSE
      )
      # `heritability_se` is a scalar SE for univariate fits and may be a
      # per-trait vector for multivariate; recycle a scalar across the h^2 rows.
      hse <- object$result$heritability_se
      if (!is.null(hse)) {
        hsev <- as.numeric(if (is.data.frame(hse)) hse$se else hse)
        if (length(hsev) == 1L) {
          hsev <- rep(hsev, nrow(h_df))
        }
        if (length(hsev) == nrow(h_df)) {
          # Surface the RAW asymptotic bounds and annotate when they cross [0, 1];
          # do not silently clamp (plotting standard 24 section 2; mirrors the engine's
          # boundary-throw discipline). A whisker crossing 0/1 is the honest signal
          # that h^2 is imprecise.
          h_df$lo <- h_df$estimate - 1.96 * hsev
          h_df$hi <- h_df$estimate + 1.96 * hsev
          experimental <- experimental || any(is.finite(hsev))
        }
      }
      boundary <- any(is.finite(h_df$lo) & (h_df$lo < 0 | h_df$hi > 1))
      df <- rbind(vc_df, h_df)
    } else {
      boundary <- FALSE
      df <- vc_df
    }
  }

  sub <- if (experimental) {
    paste0(
      "experimental +/- 1.96 SE (asymptotic, REML, not coverage-calibrated)",
      if (isTRUE(boundary)) "; h^2 CI crosses the [0,1] boundary" else ""
    )
  } else {
    NULL
  }
  hs_attach_meta(
    hs_gg_forest(
      df,
      xlab = "estimate",
      title = "Variance components and heritability",
      subtitle = sub,
      zero_line = TRUE
    ),
    type = "variance",
    interval_status = if (experimental) "experimental_asymptotic" else "none",
    notes = "variance components + per-trait h^2; SEs asymptotic/REML, not coverage-calibrated"
  )
}

# Extract a tidy breeding-values frame (id, value, pev[, trait]) from an attached
# engine `breeding_values_plot_data` payload. Rename-robust (value/breeding_value;
# pev/prediction_error_variance). Returns NULL when no usable payload is present,
# so the caller falls back to the `breeding_values()` recompute path.
hs_breeding_values_from_payload <- function(payload) {
  if (is.null(payload)) {
    return(NULL)
  }
  payload <- as.list(payload)
  pick <- function(...) {
    for (nm in c(...)) {
      if (!is.null(payload[[nm]])) {
        return(payload[[nm]])
      }
    }
    NULL
  }
  id <- pick("id", "ids")
  value <- pick("value", "values", "breeding_value", "breeding_values")
  if (is.null(id) || is.null(value) || length(id) != length(value)) {
    return(NULL)
  }
  pev <- pick("pev", "prediction_error_variance")
  trait <- pick("trait")
  out <- data.frame(
    id = as.character(id),
    value = as.numeric(value),
    stringsAsFactors = FALSE
  )
  if (!is.null(pev) && length(pev) == nrow(out)) {
    out$pev <- as.numeric(pev)
  }
  if (!is.null(trait) && length(trait) == nrow(out)) {
    out$trait <- trait
  }
  out
}

hs_autoplot_breeding_values <- function(object, ...) {
  # Auto-detect the engine `breeding_values_plot_data` payload when the bridge
  # attaches it to the fit (id/value/pev, rename-robust); otherwise assemble from
  # the `breeding_values()` + `prediction_error_variance()` extractors (recompute
  # fallback). The live R<->engine parity test (`test-plot-data-parity.R`) pins
  # that the recompute matches the engine.
  bv <- hs_breeding_values_from_payload(
    object$result$breeding_values_plot_data
  )
  from_payload <- !is.null(bv)
  if (!from_payload) {
    bv <- tryCatch(breeding_values(object), error = function(e) NULL)
    if (is.null(bv) || !all(c("id", "value") %in% names(bv))) {
      stop(
        "This `hsquared_fit` has no breeding values to plot.",
        call. = FALSE
      )
    }
  }
  bv <- bv[is.finite(bv$value), , drop = FALSE]
  has_trait <- "trait" %in% names(bv)
  # PEV band: carried in the payload (already on `bv$pev`) or looked up from the
  # `prediction_error_variance()` extractor on the recompute path.
  if (from_payload) {
    if (is.null(bv$pev)) {
      bv$pev <- NA_real_
    }
  } else {
    pev <- object$result$prediction_error_variance
    if (!is.null(pev) && all(c("id", "value") %in% names(pev))) {
      key <- if (has_trait && "trait" %in% names(pev)) {
        paste(bv$id, bv$trait)
      } else {
        bv$id
      }
      pkey <- if (has_trait && "trait" %in% names(pev)) {
        paste(pev$id, pev$trait)
      } else {
        pev$id
      }
      bv$pev <- as.numeric(pev$value)[match(key, pkey)]
    } else {
      bv$pev <- NA_real_
    }
  }
  # rank within trait (or overall)
  split_key <- if (has_trait) bv$trait else rep("1", nrow(bv))
  bv$rank <- stats::ave(
    bv$value,
    split_key,
    FUN = function(v) rank(v, ties.method = "first")
  )
  bv$lo <- bv$value - 1.96 * sqrt(bv$pev)
  bv$hi <- bv$value + 1.96 * sqrt(bv$pev)
  has_band <- any(is.finite(bv$lo))

  p <- ggplot2::ggplot(bv, ggplot2::aes(x = .data$rank, y = .data$value)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 3, colour = "grey70")
  if (has_band) {
    p <- p +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$lo, ymax = .data$hi),
        fill = "#2c6fbb",
        alpha = 0.18,
        na.rm = TRUE
      )
  }
  p <- p + ggplot2::geom_point(size = 1.2, colour = "#2c6fbb", na.rm = TRUE)
  if (has_trait) {
    p <- p + ggplot2::facet_wrap(ggplot2::vars(.data$trait), scales = "free")
  }
  sub <- if (has_band) {
    "sorted EBVs with +/- 1.96 sqrt(PEV) bands"
  } else {
    "sorted EBVs (no prediction error variances available)"
  }
  hs_attach_meta(
    p +
      ggplot2::labs(
        x = "rank",
        y = "breeding value",
        title = "Estimated breeding values",
        subtitle = sub
      ) +
      theme_hsquared(),
    type = "breeding_values",
    interval_status = if (has_band) "pev_band" else "none",
    notes = "sorted EBVs; bands are +/- 1.96 sqrt(PEV) when available"
  )
}

hs_autoplot_g_matrix <- function(object, low_h2 = 0.1, ...) {
  # Auto-detect the engine's `genetic_correlation_plot_data` payload when the
  # bridge attaches it to the fit, otherwise assemble from the stored extractors
  # (recompute fallback). The live R<->engine parity test
  # (test-plot-data-parity.R) is the drift guard: the engine preparer's
  # `genetic_correlations` is `D^-1 G D^-1`, the same correlation
  # `genetic_correlation()` reads from the fit. Both paths feed one tidy shape.
  # Consume the engine payload only when it is present AND honestly claims
  # rotation invariance (this figure stamps rotation_status="rotation_invariant",
  # so a payload that says otherwise must not pass silently); else recompute.
  pd <- object$result$genetic_correlation_plot_data
  h2 <- NULL
  rg <- NULL
  if (
    !is.null(pd) &&
      !is.null(pd$genetic_correlations) &&
      !isFALSE(pd$rotation_invariant)
  ) {
    rg <- hs_as_square_matrix(pd$genetic_correlations)
    if (!is.null(rg)) {
      if (!is.null(pd$traits)) {
        dimnames(rg) <- list(as.character(pd$traits), as.character(pd$traits))
      }
      if (!is.null(pd$heritabilities)) {
        h2 <- as.numeric(pd$heritabilities)
      }
    }
  }
  if (is.null(rg)) {
    rg <- tryCatch(genetic_correlation(object), error = function(e) NULL)
    h2 <- NULL
  }
  if (is.null(rg) || !is.matrix(rg) || nrow(rg) < 2L) {
    stop(
      "`type = \"g_matrix\"` needs a multivariate fit with a genetic ",
      "correlation matrix (rotation-invariant). Fit with ",
      "`engine_control = list(target = \"multivariate\")`.",
      call. = FALSE
    )
  }
  traits <- rownames(rg)
  if (is.null(traits)) {
    # match the engine preparer's default labels (evolvability.jl: `trait_$(i)`)
    traits <- paste0("trait_", seq_len(nrow(rg)))
  }

  # Per-trait h^2 aligned to `traits`: payload first, else the fit's heritability
  # extractor; NA where unavailable so the low-h^2 flag degrades gracefully.
  if (is.null(h2)) {
    her <- object$result$heritability
    if (!is.null(her) && all(c("term", "estimate") %in% names(her))) {
      h2 <- as.numeric(her$estimate)[match(traits, as.character(her$term))]
    }
  }
  if (!is.null(h2) && length(h2) != length(traits)) {
    h2 <- NULL # mismatched length: do not flag rather than mis-flag
  }

  df <- expand.grid(
    row = factor(traits, levels = traits),
    col = factor(traits, levels = rev(traits)),
    stringsAsFactors = FALSE
  )
  df$value <- as.numeric(rg[cbind(
    match(df$row, traits),
    match(df$col, traits)
  )])
  df$label <- formatC(df$value, format = "f", digits = 2)

  # Low-h^2 flag (plotting standard 24 section 2): an off-diagonal correlation that
  # involves a low-h^2 trait is imprecise -> mark it. The diagonal (==1) is
  # definitional, so it is never flagged. `intToUtf8` builds the glyphs from code
  # points (ASCII source, no backslash-escape ambiguity).
  dagger <- intToUtf8(0x2020)
  flagged_any <- FALSE
  if (!is.null(h2) && any(is.finite(h2))) {
    h2_row <- h2[match(as.character(df$row), traits)]
    h2_col <- h2[match(as.character(df$col), traits)]
    df$low <- as.character(df$row) != as.character(df$col) &
      ((is.finite(h2_row) & h2_row < low_h2) |
        (is.finite(h2_col) & h2_col < low_h2))
    df$label <- ifelse(df$low, paste0(df$label, dagger), df$label)
    flagged_any <- any(df$low)
  } else {
    df$low <- FALSE
  }

  base_sub <- "rotation-invariant; raw loadings are never plotted"
  sub <- if (flagged_any) {
    paste0(
      base_sub,
      "; ",
      dagger,
      " involves a low-h",
      intToUtf8(0x00b2),
      " (< ",
      low_h2,
      ") trait ",
      intToUtf8(0x2014),
      " correlation imprecise"
    )
  } else {
    base_sub
  }

  hs_attach_meta(
    ggplot2::ggplot(
      df,
      ggplot2::aes(x = .data$row, y = .data$col, fill = .data$value)
    ) +
      ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
      ggplot2::geom_text(ggplot2::aes(label = .data$label), size = 3.4) +
      ggplot2::scale_fill_gradient2(
        low = "#b2182b",
        mid = "white",
        high = "#2166ac",
        midpoint = 0,
        limits = c(-1, 1),
        name = "genetic\ncorrelation"
      ) +
      ggplot2::coord_equal() +
      ggplot2::labs(
        x = NULL,
        y = NULL,
        title = "Genetic correlation (G)",
        subtitle = sub
      ) +
      theme_hsquared(),
    type = "g_matrix",
    rotation_status = "rotation_invariant",
    notes = paste0(
      "genetic correlations only; raw factor loadings are never plotted",
      if (flagged_any) "; low-h2 cells flagged (imprecise)" else ""
    )
  )
}

# Scree of the genetic eigenstructure: rotation-invariant eigenvalues (variance
# per genetic axis) with variance-explained labels. Axis DIRECTIONS / loadings are
# never drawn -- they are rotation-arbitrary and span-ambiguous under repeated
# eigenvalues (plotting standard 24 section 2 g_geometry; the cataloged figure is the
# scree, not a loadings biplot).
hs_autoplot_g_geometry <- function(object, ...) {
  # Auto-detect the engine `genetic_pca_plot_data` payload; else recompute from
  # the fit's G via eigen_G().
  pd <- object$result$genetic_pca_plot_data
  ev <- NULL
  ve <- NULL
  axis <- NULL
  if (
    !is.null(pd) &&
      !is.null(pd$eigenvalues) &&
      !isFALSE(pd$rotation_invariant) &&
      # section 3-enforced for g_geometry: a payload that disclaims eigenstructure status
      # (i.e. carries loadings/directions) must NOT be drawn as a scree.
      !isFALSE(pd$is_eigenstructure_not_loadings)
  ) {
    ev <- as.numeric(pd$eigenvalues)
    if (!is.null(pd$variance_explained)) {
      ve <- as.numeric(pd$variance_explained)
    }
    if (!is.null(pd$axis_labels)) {
      axis <- as.character(pd$axis_labels)
    }
  }
  if (is.null(ev)) {
    eg <- tryCatch(eigen_G(object), error = function(e) NULL)
    if (is.null(eg) || is.null(eg$values) || length(eg$values) < 2L) {
      stop(
        "`type = \"g_geometry\"` needs a multivariate fit with a genetic ",
        "covariance matrix (`engine_control = list(target = \"multivariate\")`).",
        call. = FALSE
      )
    }
    ev <- as.numeric(eg$values)
  }
  if (is.null(ve) || length(ve) != length(ev)) {
    total <- sum(ev)
    ve <- if (is.finite(total) && total > 0) {
      ev / total
    } else {
      rep(NA_real_, length(ev))
    }
  }
  if (is.null(axis) || length(axis) != length(ev)) {
    axis <- paste0("PC", seq_along(ev))
  }
  df <- data.frame(
    axis = factor(axis, levels = axis),
    eigenvalue = ev,
    variance_explained = ve,
    stringsAsFactors = FALSE
  )
  # The recompute path is PSD-gated (eigen_G via hs_fit_genetic_G throws on a
  # non-PSD G), but an engine payload can carry a (small) negative eigenvalue.
  # With a negative eigenvalue the percent-variance shares are not a clean
  # partition, so omit the % labels and flag the non-PSD G honestly rather than
  # print a meaningless "116%".
  non_psd <- any(is.finite(df$eigenvalue) & df$eigenvalue < 0)
  if (non_psd) {
    df$label <- ""
  } else {
    df$label <- ifelse(
      is.finite(df$variance_explained),
      paste0(
        formatC(100 * df$variance_explained, format = "f", digits = 0),
        "%"
      ),
      ""
    )
  }
  sub <- if (non_psd) {
    paste(
      "rotation-invariant eigenvalues; G is non-positive-definite",
      "(variance shares omitted); axis directions / loadings are never shown"
    )
  } else {
    paste(
      "rotation-invariant eigenvalues (% variance explained);",
      "axis directions / loadings are never shown"
    )
  }
  hs_attach_meta(
    ggplot2::ggplot(
      df,
      ggplot2::aes(x = .data$axis, y = .data$eigenvalue)
    ) +
      ggplot2::geom_col(fill = "#2c6fbb", width = 0.7) +
      ggplot2::geom_text(
        ggplot2::aes(label = .data$label),
        vjust = -0.4,
        size = 3.2,
        na.rm = TRUE
      ) +
      ggplot2::labs(
        x = "genetic axis",
        y = "eigenvalue (variance)",
        title = "Genetic eigenstructure (G)",
        subtitle = sub
      ) +
      theme_hsquared(),
    type = "g_geometry",
    rotation_status = "rotation_invariant",
    notes = paste0(
      "rotation-invariant eigenstructure (variance per genetic axis); ",
      "axis directions/loadings not shown",
      if (non_psd) "; non-positive-definite G (variance shares omitted)" else ""
    )
  )
}

# Reaction-norm eigenfunctions psi_j(t): the rotation-invariant eigen-decomposition
# of K_g drawn as covariate functions, faceted by axis with each curve's genetic
# variance share. Eigenfunction SIGNS are arbitrary and the curves are span-ambiguous
# under repeated eigenvalues (plotting standard 24 section 2) -- caveated, never over-read.
hs_rr_use_default_payload_grid <- function(at, n) {
  nn <- suppressWarnings(as.integer(n))
  is.null(at) && length(nn) == 1L && !is.na(nn) && identical(nn, 25L)
}

hs_autoplot_rr_eigenfunctions <- function(object, at = NULL, n = 25L, ...) {
  pd <- object$result$rr_eigenfunctions_plot_data
  ve <- NULL
  if (
    hs_rr_use_default_payload_grid(at, n) &&
      !is.null(pd) &&
      !is.null(pd$eigenfunctions) &&
      !is.null(pd$covariate) &&
      !isFALSE(pd$rotation_invariant)
  ) {
    ef_mat <- as.matrix(pd$eigenfunctions) # m x k
    cov <- as.numeric(pd$covariate)
    if (!is.null(pd$variance_explained)) {
      ve <- as.numeric(pd$variance_explained)
    }
    k <- ncol(ef_mat)
    df <- data.frame(
      covariate = rep(cov, times = k),
      axis = rep(seq_len(k), each = length(cov)),
      value = as.numeric(ef_mat),
      stringsAsFactors = FALSE
    )
  } else {
    ef <- rr_eigenfunctions(object, at = at, n = n)
    df <- ef$eigenfunctions
    ve <- ef$variance_explained
  }
  k <- max(df$axis)
  if (is.null(ve) || length(ve) != k) {
    ve <- rep(NA_real_, k)
  }
  lab <- ifelse(
    is.finite(ve),
    sprintf(
      "axis %d (%s%%)",
      seq_len(k),
      formatC(100 * ve, format = "f", digits = 0)
    ),
    sprintf("axis %d", seq_len(k))
  )
  df$axis_label <- factor(lab[df$axis], levels = lab)
  hs_attach_meta(
    ggplot2::ggplot(
      df,
      ggplot2::aes(x = .data$covariate, y = .data$value)
    ) +
      ggplot2::geom_hline(yintercept = 0, linetype = 3, colour = "grey70") +
      ggplot2::geom_line(colour = "#2c6fbb", linewidth = 0.8) +
      ggplot2::facet_wrap(
        ggplot2::vars(.data$axis_label),
        ncol = 1L,
        scales = "free_y"
      ) +
      ggplot2::labs(
        x = "covariate",
        y = expression(psi[j](t)),
        title = "Reaction-norm eigenfunctions",
        subtitle = paste(
          "rotation-invariant (% genetic variance per axis);",
          "eigenfunction signs are arbitrary, span-ambiguous under",
          "repeated eigenvalues"
        )
      ) +
      theme_hsquared(),
    type = "rr_eigenfunctions",
    rotation_status = "rotation_invariant",
    notes = paste(
      "rotation-invariant eigenfunctions of supplied K_g;",
      "signs arbitrary; span-ambiguous under repeated eigenvalues"
    )
  )
}

# Reaction-norm genetic covariance (or correlation) surface over the covariate
# grid: S(s, t) = phi(s)' K_g phi(t). Supplied-K_g descriptive and
# rotation-invariant (basis-rotation-invariant genetic covariance function).
hs_autoplot_rr_surface <- function(
  object,
  at = NULL,
  n = 25L,
  correlation = FALSE,
  ...
) {
  pd <- object$result$rr_covariance_surface_plot_data
  if (
    hs_rr_use_default_payload_grid(at, n) &&
      !is.null(pd) &&
      !is.null(pd$surface) &&
      !is.null(pd$covariate)
  ) {
    cov <- as.numeric(pd$covariate)
    surf <- as.matrix(pd$surface)
    is_corr <- isTRUE(pd$is_correlation)
  } else {
    k_g <- rr_covariance(object) # rejects non-RR fits with a clear message
    pts <- hs_rr_eval_points(object, at, n)
    phi <- hs_legendre_design(pts$t, pts$order)
    surf <- phi %*% k_g %*% t(phi)
    surf <- 0.5 * (surf + t(surf))
    cov <- pts$at
    is_corr <- isTRUE(correlation)
    if (is_corr) {
      d <- diag(surf)
      if (any(d <= 0)) {
        stop(
          "`type = \"rr_surface\"` correlation is undefined: a covariate point ",
          "has non-positive genetic variance. Use `correlation = FALSE` or ",
          "evaluate inside the genetic-variance support.",
          call. = FALSE
        )
      }
      surf <- stats::cov2cor(surf)
    }
  }
  m <- length(cov)
  df <- expand.grid(i = seq_len(m), j = seq_len(m))
  df$covariate_i <- cov[df$i]
  df$covariate_j <- cov[df$j]
  df$value <- as.numeric(surf[cbind(df$i, df$j)])
  fill_scale <- if (is_corr) {
    ggplot2::scale_fill_gradient2(
      low = "#b2182b",
      mid = "white",
      high = "#2166ac",
      midpoint = 0,
      limits = c(-1, 1),
      name = "genetic\ncorrelation"
    )
  } else {
    ggplot2::scale_fill_gradient2(
      low = "#b2182b",
      mid = "white",
      high = "#2166ac",
      midpoint = 0,
      name = "genetic\ncovariance"
    )
  }
  hs_attach_meta(
    ggplot2::ggplot(
      df,
      ggplot2::aes(
        x = .data$covariate_i,
        y = .data$covariate_j,
        fill = .data$value
      )
    ) +
      ggplot2::geom_tile() +
      fill_scale +
      ggplot2::coord_equal() +
      ggplot2::labs(
        x = "covariate",
        y = "covariate",
        title = paste0(
          "Reaction-norm genetic ",
          if (is_corr) "correlation" else "covariance",
          " surface"
        ),
        subtitle = "supplied-K_g descriptive; rotation-invariant"
      ) +
      theme_hsquared(),
    type = "rr_surface",
    rotation_status = "rotation_invariant",
    notes = paste0(
      "supplied-K_g descriptive genetic ",
      if (is_corr) "correlation" else "covariance",
      " surface; rotation-invariant"
    )
  )
}

hs_autoplot_reaction_norm <- function(object, at = NULL, n = 25L, ...) {
  # Auto-detect the engine `rr_genetic_variance_plot_data` payload (covariate +
  # genetic-variance + heritability trajectories) when the user has not asked for a
  # custom grid; else recompute from K_g via rr_genetic_variance()/rr_heritability()
  # (which also reject non-RR fits with a clear message). Rename-robust: accept
  # either `value` (the #93-agreed field) or the engine's legacy
  # `genetic_variance`.
  pd <- object$result$rr_genetic_variance_plot_data
  gv_payload <- NULL
  if (!is.null(pd)) {
    gv_payload <- if (!is.null(pd$value)) pd$value else pd$genetic_variance
  }
  if (
    hs_rr_use_default_payload_grid(at, n) &&
      !is.null(pd) &&
      !is.null(pd$covariate) &&
      !is.null(gv_payload) &&
      !is.null(pd$heritability)
  ) {
    df <- rbind(
      data.frame(
        covariate = as.numeric(pd$covariate),
        value = as.numeric(gv_payload),
        panel = "genetic variance",
        stringsAsFactors = FALSE
      ),
      data.frame(
        covariate = as.numeric(pd$covariate),
        value = as.numeric(pd$heritability),
        panel = "heritability",
        stringsAsFactors = FALSE
      )
    )
  } else {
    vg <- rr_genetic_variance(object, at = at, n = n)
    h2 <- rr_heritability(object, at = at, n = n)
    df <- rbind(
      data.frame(
        covariate = vg$covariate,
        value = vg$value,
        panel = "genetic variance",
        stringsAsFactors = FALSE
      ),
      data.frame(
        covariate = h2$covariate,
        value = h2$value,
        panel = "heritability",
        stringsAsFactors = FALSE
      )
    )
  }
  df$panel <- factor(df$panel, levels = c("genetic variance", "heritability"))
  hs_attach_meta(
    ggplot2::ggplot(df, ggplot2::aes(x = .data$covariate, y = .data$value)) +
      ggplot2::geom_line(colour = "#2c6fbb", linewidth = 0.8) +
      ggplot2::geom_point(size = 1.4, colour = "#2c6fbb") +
      ggplot2::facet_wrap(
        ggplot2::vars(.data$panel),
        scales = "free_y",
        ncol = 1L
      ) +
      ggplot2::labs(
        x = "covariate",
        y = NULL,
        title = "Reaction norm (random regression)",
        subtitle = paste(
          "genetic-variance and heritability trajectories;",
          "h2(t) experimental (homogeneous residual, no PE term)"
        )
      ) +
      theme_hsquared(),
    type = "reaction_norm",
    interval_status = "descriptive",
    rotation_status = "rotation_invariant",
    notes = "supplied-K_g descriptive trajectories; h2(t) can overstate without a PE term"
  )
}

# --- autoplot.hs_gwas (Manhattan) ------------------------------------------

#' @rdname hsquared-autoplot
#' @exportS3Method ggplot2::autoplot
autoplot.hs_gwas <- function(object, type = c("manhattan", "qq"), ...) {
  hs_require_ggplot2()
  type <- match.arg(type)
  method <- attr(object, "scan_method") %||% "mixed"
  df <- as.data.frame(object)
  if (!all(c("p_value") %in% names(df))) {
    stop(
      "This `hs_gwas` object has no `p_value` column to plot.",
      call. = FALSE
    )
  }
  switch(
    type,
    manhattan = hs_autoplot_manhattan(df, method),
    qq = hs_autoplot_qq(df, method)
  )
}

# Honest-status suffix when the scan ignored the pedigree relationship.
hs_gwas_method_note <- function(method) {
  if (identical(method, "single")) {
    "; relatedness-UNcorrected (single-marker, OLS)"
  } else if (identical(method, "loco")) {
    "; LOCO genomic correction (pedigree-estimated VCs)"
  } else {
    ""
  }
}

hs_autoplot_manhattan <- function(df, method = "mixed") {
  df$index <- seq_len(nrow(df))
  df$neglog10p <- -log10(pmax(df$p_value, .Machine$double.xmin))
  m <- nrow(df)
  bonf <- -log10(0.05 / m)
  hs_attach_meta(
    ggplot2::ggplot(df, ggplot2::aes(x = .data$index, y = .data$neglog10p)) +
      ggplot2::geom_hline(
        yintercept = bonf,
        linetype = 2,
        colour = "#b2182b"
      ) +
      ggplot2::geom_point(size = 1.3, colour = "#2c6fbb") +
      ggplot2::annotate(
        "text",
        x = 1,
        y = bonf,
        vjust = -0.5,
        hjust = 0,
        size = 3,
        colour = "#b2182b",
        label = "Bonferroni 0.05"
      ) +
      ggplot2::labs(
        x = "marker",
        y = expression(-log[10](p)),
        title = "Marker scan (Manhattan)",
        subtitle = paste0(
          "EXPERIMENTAL: nominal Wald p-values, NOT genome-wide ",
          "calibrated (no R threshold activation)",
          hs_gwas_method_note(method)
        )
      ) +
      theme_hsquared(),
    type = "manhattan",
    source = "gwas",
    interval_status = "uncalibrated",
    notes = paste0(
      "nominal Wald p-values; Bonferroni line is visual only; not genome-wide ",
      "calibrated; Julia fixed-panel smoke is not an R threshold activation",
      hs_gwas_method_note(method)
    )
  )
}

# QQ of the scan p-values against the uniform null, with the genomic-inflation
# lambda_GC as a diagnostic annotation. Pure-R from the p-values (no engine).
hs_autoplot_qq <- function(df, method = "mixed") {
  p <- pmax(as.numeric(df$p_value), .Machine$double.xmin)
  m <- length(p)
  qq <- data.frame(
    expected = -log10(stats::ppoints(m)),
    observed = -log10(sort(p))
  )
  lambda_gc <- stats::median(stats::qchisq(1 - p, df = 1)) /
    stats::qchisq(0.5, df = 1)
  lab <- formatC(lambda_gc, format = "f", digits = 2)
  hs_attach_meta(
    ggplot2::ggplot(
      qq,
      ggplot2::aes(x = .data$expected, y = .data$observed)
    ) +
      ggplot2::geom_abline(
        slope = 1,
        intercept = 0,
        linetype = 2,
        colour = "#b2182b"
      ) +
      ggplot2::geom_point(size = 1.3, colour = "#2c6fbb") +
      ggplot2::labs(
        x = expression(expected ~ -log[10](p)),
        y = expression(observed ~ -log[10](p)),
        title = "Marker scan (QQ)",
        subtitle = paste0(
          "EXPERIMENTAL: nominal Wald p-values, NOT genome-wide calibrated ",
          "(no R threshold activation); genomic inflation lambda_GC = ",
          lab,
          " (diagnostic only; >1 may reflect structure/polygenicity, ",
          "not corrected)",
          hs_gwas_method_note(method)
        )
      ) +
      theme_hsquared(),
    type = "qq",
    source = "gwas",
    interval_status = "uncalibrated",
    notes = paste0(
      "QQ of nominal Wald p-values; y=x is the null; lambda_GC = ",
      lab,
      " is diagnostic only; not genome-wide calibrated"
    )
  )
}

# --- recovery forest (validation-study data frame) -------------------------

#' Forest plot of a known-truth recovery study
#'
#' `r lifecycle::badge("experimental")`
#'
#' Visualizes a bias +/- 2*MCSE recovery study (e.g.
#' `data-raw/multivariate-recovery-study.R`): each target's bias with its
#' +/- 2*MCSE interval and a zero-bias reference line. Targets whose interval
#' covers zero show "no detectable bias".
#'
#' @param data A data frame with `target`, `bias`, and `mcse` columns.
#' @return A `ggplot` object.
#' @export
hs_recovery_forest <- function(data) {
  hs_require_ggplot2()
  if (!all(c("target", "bias", "mcse") %in% names(data))) {
    stop(
      "`data` must have `target`, `bias`, and `mcse` columns.",
      call. = FALSE
    )
  }
  df <- data.frame(
    term = as.character(data$target),
    estimate = as.numeric(data$bias),
    lo = as.numeric(data$bias) - 2 * as.numeric(data$mcse),
    hi = as.numeric(data$bias) + 2 * as.numeric(data$mcse),
    stringsAsFactors = FALSE
  )
  df$covers0 <- df$lo <= 0 & df$hi >= 0
  hs_attach_meta(
    hs_gg_forest(
      df,
      xlab = "bias (mean(hat) - truth)",
      title = "Known-truth recovery (bias +/- 2 MCSE)",
      subtitle = "intervals covering 0 indicate no detectable bias",
      zero_line = TRUE
    ),
    type = "recovery_forest",
    source = "study",
    interval_status = "mcse_band",
    notes = "bias +/- 2 MCSE; interval covering 0 = no detectable bias"
  )
}

# --- helper ----------------------------------------------------------------

# Coerce a bridge-marshalled correlation field to a square matrix. JuliaCall may
# hand a NamedTuple matrix field back as a column-major numeric vector rather than
# a 2D matrix; restore p x p from a perfect-square length. Returns NULL when the
# input cannot be a >=2x2 square matrix, so the caller falls back to recompute.
hs_as_square_matrix <- function(x) {
  if (is.matrix(x)) {
    return(x)
  }
  v <- suppressWarnings(as.numeric(x))
  n <- sqrt(length(v))
  if (length(v) >= 4L && is.finite(n) && n == floor(n)) {
    return(matrix(v, nrow = n, ncol = n))
  }
  NULL
}

hs_require_ggplot2 <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "The hsquared plotting layer requires the `ggplot2` package.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# Attach the cross-lane honest-status meta (the `13-plotting-layer.md` section 3
# contract, mirroring gllvmTMB/drmTMB): machine-readable type / source /
# interval / rotation status + caveat notes, alongside the human-readable
# subtitle, so the honest-status guardrails are inspectable, not just printed.
hs_attach_meta <- function(
  p,
  type,
  source = "fit",
  interval_status = "none",
  rotation_status = "not_applicable",
  notes = NULL
) {
  attr(p, "hsquared_meta") <- list(
    type = type,
    source = source,
    interval_status = interval_status,
    rotation_status = rotation_status,
    notes = notes
  )
  p
}
