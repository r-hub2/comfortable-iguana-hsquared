#' G-matrix geometry and evolvability
#'
#' `r lifecycle::badge("experimental")`
#'
#' These extractors summarise the genetic variance-covariance matrix `G` of an
#' multivariate `hsquared_fit` (a `cbind()` response, which fits on the default
#' path) through its
#' **rotation-invariant** geometry, following Hansen & Houle (2008). They are
#' defined on `G` itself (not on factor loadings), so they are well defined for
#' any multivariate fit -- unstructured, diagonal, or (when bridged) low-rank /
#' factor-analytic -- and do not depend on a loading rotation convention.
#'
#' - `eigen_G()` returns the genetic eigenstructure: `values` (the variance along
#'   each genetic principal axis, descending) and `vectors` (the genetic
#'   principal components, sign-canonicalised so the largest-magnitude element of
#'   each is positive).
#' - `g_max()` returns the leading genetic axis: its `eigenvalue` and
#'   `eigenvector` (the direction of maximum evolvability).
#' - `mean_evolvability()` is the average evolvability over random selection
#'   directions, `tr(G) / t`.
#' - `evolvability(fit, direction)` is `e(beta) = beta'G beta` (unit `beta`): the additive
#'   genetic variance available to directional selection along `direction`.
#' - `variance_along_gradient(fit, direction, normalize)` is `beta'G beta` for a unit
#'   `direction` (`normalize = TRUE`, equal to `evolvability()`) or for the raw
#'   `direction` (`normalize = FALSE`).
#' - `respondability(fit, direction)` is `||G beta||`: the length of the response.
#' - `conditional_evolvability(fit, direction)` is `1 / (beta' G^-1 beta)`: evolvability
#'   when all other directions are held under stabilising selection (requires a
#'   positive-definite `G`).
#' - `autonomy(fit, direction)` is `conditional_evolvability / evolvability`
#'   in `[0, 1]` (requires a positive-definite `G`).
#'
#' These are computed in R from the fitted `genetic_covariance(fit)` and match
#' the engine's `evolvability.jl` definitions (verified by a skip-guarded live
#' parity test). They are descriptive geometry of the estimated `G`; they carry
#' the same experimental, REML-only, not-coverage-calibrated status as the
#' multivariate fit itself and report no standard errors.
#'
#' @param object An `hsquared_fit` from the multivariate model.
#' @param direction A numeric vector of selection gradients, one per trait. It
#'   is normalised to unit length internally; only its direction matters.
#' @param ... Unused.
#'
#' @return A scalar for the directional metrics; a list for `eigen_G()`
#'   (`values`, `vectors`) and `g_max()` (`eigenvalue`, `eigenvector`).
#'
#' @references Hansen, T. F., & Houle, D. (2008). Measuring and comparing
#'   evolvability and constraint in multivariate characters. *Journal of
#'   Evolutionary Biology*, 21(5), 1201-1219.
#' @name g_matrix_geometry
NULL

# Extract and validate the genetic covariance G of a multivariate fit. Mirrors
# the engine's symmetric-PSD (or PD) check before computing geometry.
hs_fit_genetic_G <- function(object, require_pd = FALSE) {
  if (!inherits(object, "hsquared_fit")) {
    stop("`object` must be an `hsquared_fit` object.", call. = FALSE)
  }
  G <- object$result$genetic_covariance
  if (is.null(G)) {
    stop(
      "G-matrix geometry requires a multivariate fit that reports a ",
      "`genetic_covariance()` matrix -- a `cbind(trait1, trait2, ...)` response ",
      "with `animal(1 | id, pedigree = ped)`.",
      call. = FALSE
    )
  }
  G <- as.matrix(G)
  storage.mode(G) <- "double"
  if (nrow(G) != ncol(G)) {
    stop("The genetic covariance matrix must be square.", call. = FALSE)
  }
  if (any(!is.finite(G))) {
    stop("The genetic covariance matrix must be finite.", call. = FALSE)
  }
  gscale <- max(1, max(abs(G)))
  if (!isTRUE(all.equal(G, t(G), tolerance = 1e-8 * gscale))) {
    stop("The genetic covariance matrix must be symmetric.", call. = FALSE)
  }
  G <- (G + t(G)) / 2
  ev <- eigen(G, symmetric = TRUE, only.values = TRUE)$values
  escale <- max(1, max(abs(ev)))
  if (min(ev) < -1e-8 * escale) {
    stop(
      "The genetic covariance matrix must be positive semidefinite.",
      call. = FALSE
    )
  }
  if (require_pd && min(ev) <= 1e-10 * escale) {
    stop(
      "This metric requires a positive-definite genetic covariance matrix (it ",
      "inverts G); a singular / reduced-rank G has no conditional evolvability ",
      "or autonomy.",
      call. = FALSE
    )
  }
  G
}

hs_normalize_direction <- function(direction, t) {
  if (!is.numeric(direction) || length(direction) != t) {
    stop(
      "`direction` must be a numeric vector with one entry per trait (length ",
      t,
      ").",
      call. = FALSE
    )
  }
  if (any(!is.finite(direction))) {
    stop("`direction` must contain only finite values.", call. = FALSE)
  }
  nrm <- sqrt(sum(direction^2))
  if (nrm <= 0) {
    stop("`direction` must be a nonzero direction.", call. = FALSE)
  }
  direction / nrm
}

hs_g_eigen <- function(G) {
  e <- eigen(G, symmetric = TRUE)
  ord <- order(e$values, decreasing = TRUE)
  values <- e$values[ord]
  vectors <- e$vectors[, ord, drop = FALSE]
  for (j in seq_len(ncol(vectors))) {
    k <- which.max(abs(vectors[, j]))
    if (vectors[k, j] < 0) {
      vectors[, j] <- -vectors[, j]
    }
  }
  list(values = values, vectors = vectors)
}

#' @rdname g_matrix_geometry
#' @export
evolvability <- function(object, direction, ...) {
  UseMethod("evolvability")
}

#' @export
evolvability.default <- function(object, direction, ...) {
  stop("`evolvability()` requires an `hsquared_fit` object.", call. = FALSE)
}

#' @export
evolvability.hsquared_fit <- function(object, direction, ...) {
  G <- hs_fit_genetic_G(object)
  b <- hs_normalize_direction(direction, nrow(G))
  max(0, as.numeric(crossprod(b, G %*% b)))
}

#' @rdname g_matrix_geometry
#' @export
respondability <- function(object, direction, ...) {
  UseMethod("respondability")
}

#' @export
respondability.default <- function(object, direction, ...) {
  stop("`respondability()` requires an `hsquared_fit` object.", call. = FALSE)
}

#' @export
respondability.hsquared_fit <- function(object, direction, ...) {
  G <- hs_fit_genetic_G(object)
  b <- hs_normalize_direction(direction, nrow(G))
  sqrt(sum((G %*% b)^2))
}

#' @rdname g_matrix_geometry
#' @export
conditional_evolvability <- function(object, direction, ...) {
  UseMethod("conditional_evolvability")
}

#' @export
conditional_evolvability.default <- function(object, direction, ...) {
  stop(
    "`conditional_evolvability()` requires an `hsquared_fit` object.",
    call. = FALSE
  )
}

#' @export
conditional_evolvability.hsquared_fit <- function(object, direction, ...) {
  G <- hs_fit_genetic_G(object, require_pd = TRUE)
  b <- hs_normalize_direction(direction, nrow(G))
  1 / as.numeric(crossprod(b, solve(G, b)))
}

#' @rdname g_matrix_geometry
#' @export
autonomy <- function(object, direction, ...) {
  UseMethod("autonomy")
}

#' @export
autonomy.default <- function(object, direction, ...) {
  stop("`autonomy()` requires an `hsquared_fit` object.", call. = FALSE)
}

#' @export
autonomy.hsquared_fit <- function(object, direction, ...) {
  G <- hs_fit_genetic_G(object, require_pd = TRUE)
  b <- hs_normalize_direction(direction, nrow(G))
  e <- as.numeric(crossprod(b, G %*% b))
  cond <- 1 / as.numeric(crossprod(b, solve(G, b)))
  cond / e
}

#' @rdname g_matrix_geometry
#' @export
mean_evolvability <- function(object, ...) {
  UseMethod("mean_evolvability")
}

#' @export
mean_evolvability.default <- function(object, ...) {
  stop(
    "`mean_evolvability()` requires an `hsquared_fit` object.",
    call. = FALSE
  )
}

#' @export
mean_evolvability.hsquared_fit <- function(object, ...) {
  G <- hs_fit_genetic_G(object)
  sum(diag(G)) / nrow(G)
}

#' @rdname g_matrix_geometry
#' @export
g_max <- function(object, ...) {
  UseMethod("g_max")
}

#' @export
g_max.default <- function(object, ...) {
  stop("`g_max()` requires an `hsquared_fit` object.", call. = FALSE)
}

#' @export
g_max.hsquared_fit <- function(object, ...) {
  e <- hs_g_eigen(hs_fit_genetic_G(object))
  list(eigenvalue = e$values[[1L]], eigenvector = e$vectors[, 1L])
}

#' @rdname g_matrix_geometry
#' @param normalize For `variance_along_gradient()`, whether to scale
#'   `direction` to unit length first (`TRUE`, the default, which then equals
#'   `evolvability()`); `FALSE` uses `direction` as given, returning the genetic
#'   variance in the raw (un-normalised) selection gradient `beta'G beta`.
#' @export
variance_along_gradient <- function(object, direction, normalize = TRUE, ...) {
  UseMethod("variance_along_gradient")
}

#' @export
variance_along_gradient.default <- function(
  object,
  direction,
  normalize = TRUE,
  ...
) {
  stop(
    "`variance_along_gradient()` requires an `hsquared_fit` object.",
    call. = FALSE
  )
}

#' @export
variance_along_gradient.hsquared_fit <- function(
  object,
  direction,
  normalize = TRUE,
  ...
) {
  G <- hs_fit_genetic_G(object)
  t <- nrow(G)
  if (isTRUE(normalize)) {
    b <- hs_normalize_direction(direction, t)
  } else {
    if (
      !is.numeric(direction) ||
        length(direction) != t ||
        any(!is.finite(direction))
    ) {
      stop(
        "`direction` must be a finite numeric vector with one entry per trait ",
        "(length ",
        t,
        ").",
        call. = FALSE
      )
    }
    b <- direction
  }
  max(0, as.numeric(crossprod(b, G %*% b)))
}
