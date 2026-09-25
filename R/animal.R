#' Animal-model formula marker
#'
#' `r lifecycle::badge("experimental")`
#'
#' `animal()` is the additive-genetic term in an `hsquared()` formula. Write
#' `animal(1 | id, pedigree = ped)`, or `animal(1 | id)` when `data` is an
#' [hs_data()] object with a pedigree component. Fitting happens in
#' [hsquared()], not here: calling `animal()` on its own is a syntax marker
#' and returns `NULL`.
#'
#' @param formula A random-effect expression. The v0.1 parser accepts
#'   `1 | id`.
#' @param pedigree A pedigree data frame with individual, sire, and dam columns.
#'   Optional only when the enclosing [hsquared()] or [model_spec()] call uses
#'   `data = hs_data(..., pedigree = ...)`. In v0.1 this route rejects any row
#'   with the same known sire and dam (selfing); no argument here reaches the
#'   engine's existing `allow_selfing` flag. The current workaround for a
#'   selfing or hermaphroditic pedigree is `relmat(1 | id, K = A)` with a
#'   hand-built or `AGHmatrix`-built relationship matrix `A` (an experimental
#'   route; see `?hs_control`).
#' @param ... Reserved for future syntax such as `cov =`.
#'
#' @return `NULL`, invisibly. The call is interpreted by [hsquared()] when it
#'   appears inside a model formula.
#' @examples
#' # n = 4 is a syntax demo, not a number for a paper.
#' ped <- data.frame(
#'   id = c("sire", "dam", "off1", "off2"),
#'   sire = c(NA, NA, "sire", "sire"),
#'   dam = c(NA, NA, "dam", "dam")
#' )
#' dat <- data.frame(
#'   id = c("sire", "dam", "off1", "off2"),
#'   sex = c("m", "f", "m", "f"),
#'   weight = c(42, 38, 40, 37)
#' )
#' hsquared(
#'   weight ~ sex + animal(1 | id, pedigree = ped),
#'   data = dat,
#'   control = hs_control(engine = "validate")
#' )
#' model_spec(
#'   weight ~ sex + animal(1 | id, pedigree = ped),
#'   data = dat
#' )
#' @export
animal <- function(formula, pedigree = NULL, ...) {
  invisible(NULL)
}
