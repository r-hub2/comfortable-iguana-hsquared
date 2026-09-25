#' Diagnostic plots for a fitted animal model
#'
#' A base-graphics diagnostic plot for `hsquared_fit` objects. Two panels are
#' available:
#'
#' * `type = "variance"` (default) plots the estimated variance components as
#'   points. When the fit carries the **experimental** variance-component
#'   standard errors (see [variance_component_standard_errors()]), it adds
#'   approximate `+/- 1.96 * SE` whiskers and labels the panel experimental;
#'   those intervals are asymptotic, REML-only, and not coverage-calibrated. The
#'   governing coverage study found under-coverage for the sigma^2_a delta/Wald
#'   interval (worst interior coverage 0.897 at h^2 = 0.5), so no conservative
#'   direction is claimed for the plotted variance-component SE.
#' * `type = "residuals"` plots residuals against fitted values (with a zero
#'   reference line), when the fit carries fitted values and a response.
#'
#' No plotting dependency is added; this uses base graphics only.
#'
#' @param x An `hsquared_fit` object.
#' @param type Which panel to draw: `"variance"` or `"residuals"`.
#' @param ... Passed to the underlying base-graphics call.
#'
#' @return `x`, invisibly.
#' @export
plot.hsquared_fit <- function(x, type = c("variance", "residuals"), ...) {
  type <- match.arg(type)
  if (identical(type, "variance")) {
    hs_plot_variance(x, ...)
  } else {
    hs_plot_residuals(x, ...)
  }
  invisible(x)
}

hs_plot_variance <- function(x, ...) {
  vc <- x$result$variance_components
  if (is.null(vc) || !all(c("component", "estimate") %in% names(vc))) {
    stop(
      "This `hsquared_fit` has no variance-component estimates to plot.",
      call. = FALSE
    )
  }
  est <- as.numeric(vc$estimate)
  comp <- as.character(vc$component)
  n <- length(est)

  se <- x$result$variance_component_se
  has_se <- !is.null(se) && all(c("component", "se") %in% names(se))
  lo <- hi <- NULL
  if (has_se) {
    idx <- match(comp, as.character(se$component))
    se_v <- as.numeric(se$se)[idx]
    lo <- est - 1.96 * se_v
    hi <- est + 1.96 * se_v
  }

  ylim <- range(c(est, lo, hi, 0), na.rm = TRUE)
  main <- if (has_se) {
    "Variance components (experimental +/- 1.96 SE; not coverage-calibrated; asymptotic REML)"
  } else {
    "Variance components"
  }
  graphics::plot(
    seq_len(n),
    est,
    xlim = c(0.5, n + 0.5),
    ylim = ylim,
    xaxt = "n",
    xlab = "",
    ylab = "variance",
    main = main,
    pch = 19,
    ...
  )
  graphics::axis(1, at = seq_len(n), labels = comp)
  if (has_se) {
    ok <- !is.na(lo) & !is.na(hi)
    graphics::segments(seq_len(n)[ok], lo[ok], seq_len(n)[ok], hi[ok])
  }
  graphics::abline(h = 0, col = "grey70", lty = 3)
  invisible(x)
}

hs_plot_residuals <- function(x, ...) {
  fitted <- x$result$predictions$.fitted
  y <- x$payload$y
  if (is.null(fitted) || is.null(y) || length(fitted) != length(y)) {
    stop(
      "This `hsquared_fit` has no fitted values and response to plot ",
      "residuals from.",
      call. = FALSE
    )
  }
  resid <- as.numeric(y) - as.numeric(fitted)
  graphics::plot(
    as.numeric(fitted),
    resid,
    xlab = "fitted",
    ylab = "residual",
    main = "Residuals vs fitted",
    pch = 19,
    ...
  )
  graphics::abline(h = 0, col = "grey70", lty = 3)
  invisible(x)
}
