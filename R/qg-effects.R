#' Quantitative-genetic formula markers
#'
#' `r lifecycle::badge("experimental")`
#'
#' These functions provide readable formula vocabulary for standard
#' quantitative-genetic, parental, inheritance, and custom-kernel models. Called
#' directly they are inert (they return `NULL`); they take meaning only inside an
#' [hsquared()] formula. `permanent()`, `common_env()`, and `maternal_genetic()`
#' now fit through an opt-in, experimental engine path (`engine = "julia"`,
#' REML-only, not the default, mirroring a `partial` validation gate) as the
#' second random effect alongside `animal()`. `metafounder()` now fits through an
#' opt-in, supplied-variance, experimental engine path (`engine = "julia"`,
#' `target = "metafounder"`); `Gamma` and the variance components are supplied,
#' not estimated. `relmat()` and `precision()` now fit through an opt-in,
#' EXPERIMENTAL supplied-relationship primary engine path (`engine = "julia"`,
#' REML-only, not the default, mirroring a `partial` validation gate):
#' `relmat(1 | id, K = K)` supplies a dense symmetric positive-definite
#' relationship/covariance matrix K (the parser marshals the inverse
#' `Kinv = solve(K)`; `relmat()` also accepts a directly supplied inverse via
#' `Kinv =`/`Q =`, no solve), and `precision(1 | id, Q = Q)` supplies the
#' precision (inverse) directly (`Kinv =` also accepted); both fit through the
#' SAME supplied-relationship-inverse path as `genomic()`'s `Ginv`, and the
#' supplied matrix is provenance, not an estimate. The remaining markers
#' (paternal/maternal-environment, dominance, epistasis, cytoplasmic, imprinting,
#' genetic groups / unknown-parent-groups, and inbreeding) are still inert syntax
#' reservations that the parser rejects with a planned-not-implemented message.
#'
#' Some markers use generic names (e.g. `group()`, `inbreeding()`). They are
#' formula-only tokens detected by call head and are not meant to be called
#' directly, so attaching `hsquared` alongside a package that exports a
#' same-named function (e.g. `pedigreemm::inbreeding()`) may print a masking
#' message; this is expected and harmless.
#'
#' @param formula A random-effect expression such as `1 | id`.
#' @param pedigree A pedigree data frame for future parental and relationship
#'   effects.
#' @param D,Dinv Dominance relationship or precision matrices.
#' @param E,Einv Epistatic relationship or precision matrices.
#' @param K,Kinv,Q User-supplied relationship or precision matrices.
#' @param group Animal-to-metafounder group assignment, aligned to the normalized
#'   pedigree IDs for the opt-in supplied-variance `metafounder()` bridge.
#' @param Gamma A supplied metafounder relationship matrix (an `m`-by-`m`
#'   covariance over the `m` metafounder pseudo-populations; Legarra et al.
#'   2015). Supplied, not estimated; used by the opt-in supplied-variance
#'   `metafounder()` bridge.
#' @param parent Planned parent-of-origin side for imprinting effects.
#' @param ... Reserved for future syntax.
#'
#' @return `NULL`, invisibly. Calls are interpreted by [hsquared()] when they
#'   appear inside model formulas.
#' @name qg_effect_markers
NULL

#' @rdname qg_effect_markers
#' @export
permanent <- function(formula, ...) {
  invisible(NULL)
}

#' @rdname qg_effect_markers
#' @export
common_env <- function(formula, ...) {
  invisible(NULL)
}

#' @rdname qg_effect_markers
#' @export
maternal_genetic <- function(formula, pedigree = NULL, ...) {
  invisible(NULL)
}

#' @rdname qg_effect_markers
#' @export
maternal_env <- function(formula, ...) {
  invisible(NULL)
}

#' @rdname qg_effect_markers
#' @export
paternal_genetic <- function(formula, pedigree = NULL, ...) {
  invisible(NULL)
}

#' @rdname qg_effect_markers
#' @export
paternal_env <- function(formula, ...) {
  invisible(NULL)
}

#' @rdname qg_effect_markers
#' @export
cytoplasmic <- function(formula, ...) {
  invisible(NULL)
}

#' @rdname qg_effect_markers
#' @export
imprinting <- function(
  formula,
  pedigree = NULL,
  parent = c("maternal", "paternal"),
  ...
) {
  invisible(NULL)
}

#' @rdname qg_effect_markers
#' @export
dominance <- function(formula, pedigree = NULL, D = NULL, Dinv = NULL, ...) {
  invisible(NULL)
}

#' @rdname qg_effect_markers
#' @export
epistasis <- function(formula, pedigree = NULL, E = NULL, Einv = NULL, ...) {
  invisible(NULL)
}

#' @rdname qg_effect_markers
#' @export
relmat <- function(formula, K = NULL, Kinv = NULL, Q = NULL, ...) {
  invisible(NULL)
}

#' @rdname qg_effect_markers
#' @export
precision <- function(formula, Q = NULL, ...) {
  invisible(NULL)
}

#' @rdname qg_effect_markers
#' @export
group <- function(formula, ...) {
  invisible(NULL)
}

#' @rdname qg_effect_markers
#' @export
unknown_parent_group <- function(formula, ...) {
  invisible(NULL)
}

#' @rdname qg_effect_markers
#' @export
metafounder <- function(
  formula,
  pedigree = NULL,
  group = NULL,
  Gamma = NULL,
  ...
) {
  invisible(NULL)
}

#' @rdname qg_effect_markers
#' @export
inbreeding <- function(formula, ...) {
  invisible(NULL)
}
