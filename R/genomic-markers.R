#' Genomic and QTL formula markers
#'
#' `r lifecycle::badge("experimental")`
#'
#' These functions provide readable formula vocabulary for genomic, single-step,
#' marker-effect, GWAS, and QTL/eQTL models. Called directly they are inert (they
#' return `NULL`); they take meaning only inside an [hsquared()] formula.
#' `genomic()` fits a narrow Gaussian REML random-intercept model only through
#' the explicit `engine = "julia", target = "genomic"` path. `single_step()`
#' also remains opt-in and experimental.
#' The reported genomic coefficient is `sigma_g2 / (sigma_g2 + sigma_e2)` on
#' the declared relationship scale, not generally an average marginal
#' phenotypic-variance fraction or pedigree/population heritability.
#' Accepted genomic forms are `genomic(1 | id, Ginv = Ginv)` or
#' `genomic(1 | id, markers = M)` (GREML, or SNP-BLUP via
#' `target = "snp_blup"`), and `single_step(1 | id, Hinv = Hinv)` (a precomputed
#' inverse) or `single_step(1 | id, pedigree = ped, markers = M)` (the engine
#' constructs `H^-1` from the pedigree + genotyped-subset markers via
#' `target = "single_step_construct"`). Adding `group` plus supplied `Gamma`
#' fits a supplied-`Gamma`, validation-scale `H^Gamma` path through
#' `target = "metafounder_single_step"`; `Gamma` is supplied, not estimated.
#' When `data` is an [hs_data()] container that bundles a pedigree and genotypes,
#' `single_step(1 | id)` resolves both from the bundle (the `animal(1 | id)`
#' precedent), so neither `pedigree =` nor `markers =` is required; explicit
#' arguments override the bundle. The remaining markers (`markers()`,
#' `marker_scan()`, `qtl_scan()`) are still inert syntax reservations that the
#' parser rejects with a planned-not-implemented message.
#'
#' @param formula A random-effect expression such as `1 | id`.
#' @param G,Ginv,H,Hinv Relationship or precision matrices for future genomic
#'   and single-step models.
#' @param M A marker or dosage matrix for future marker-effect and scan models.
#' @param map A marker map for future marker scans.
#' @param model Planned marker-effect mode.
#' @param position A chromosome-position table or variable for future QTL
#'   scans.
#' @param genotype_probs Genotype probabilities for future interval/QTL scans.
#' @param ... Reserved for future syntax.
#'
#' @return `NULL`, invisibly. Calls are interpreted by [hsquared()] when they
#'   appear inside model formulas.
#' @name genomic_markers
NULL

#' @rdname genomic_markers
#' @export
genomic <- function(formula, G = NULL, Ginv = NULL, markers = NULL, ...) {
  invisible(NULL)
}

#' @rdname genomic_markers
#' @param pedigree A pedigree data frame (`id`, `sire`, `dam`) for the
#'   single-step `H^-1` *construction* path
#'   (`single_step(1 | id, pedigree = ped, markers = M)`), in place of a
#'   precomputed `Hinv`.
#' @param markers An individual-by-marker numeric dosage matrix with rows named
#'   by genotyped ID. Entries must be the biallelic allele count `0`/`1`/`2`
#'   (or an imputed dosage in the closed interval `[0, 2]`); a matrix in the
#'   centered `-1`/`0`/`1` convention must be recoded (add `1`) before it is
#'   passed in. For `genomic()`, the engine uses sample allele frequencies,
#'   VanRaden method 1, and ridge `0.01`. For `single_step()`, this is the
#'   genotyped-subset marker matrix used to construct the genomic relationship.
#' @param group Animal-to-metafounder group labels for the future `H^Gamma`
#'   single-step path; must be supplied together with `Gamma`.
#' @param Gamma A supplied metafounder relationship matrix for the future
#'   `H^Gamma` single-step path; supplied, not estimated, and must be supplied
#'   together with `group`.
#' @param tau,omega,blend_weight,ridge Single-step construction tuning knobs
#'   (Aguilar et al. 2010); defaults `tau = omega = 1`, `blend_weight = ridge = 0`.
#' @export
single_step <- function(
  formula,
  H = NULL,
  Hinv = NULL,
  pedigree = NULL,
  markers = NULL,
  group = NULL,
  Gamma = NULL,
  tau = 1,
  omega = 1,
  blend_weight = 0,
  ridge = 0,
  ...
) {
  invisible(NULL)
}

#' @rdname genomic_markers
#' @export
markers <- function(M, model = c("random", "fixed", "scan"), ...) {
  invisible(NULL)
}

#' @rdname genomic_markers
#' @export
marker_scan <- function(M, map = NULL, ...) {
  invisible(NULL)
}

#' @rdname genomic_markers
#' @export
qtl_scan <- function(position, genotype_probs = NULL, ...) {
  invisible(NULL)
}
