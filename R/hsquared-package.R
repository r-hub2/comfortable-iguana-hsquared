#' hsquared: R Interface for Julia-Backed Quantitative-Genetic Models
#'
#' The hsquared package is the R-facing interface for heritability,
#' breeding-value, G-matrix, and inheritance-structured quantitative-genetic
#' models backed by the HSquared.jl Julia engine.
#'
#' v0.1 fits the univariate Gaussian animal model
#' `y ~ fixed + animal(1 | id, pedigree = ped)` by REML (average-information)
#' through the HSquared.jl engine: the default `hsquared()` call fits when a
#' local Julia and `HSquared.jl` are available, and otherwise errors with
#' install guidance. A `cbind()` multivariate Gaussian response also routes on
#' that default path; multivariate is covered at validation scale and stays experimental.
#' Narrow genomic GREML is an explicit experimental target (`target = "genomic"`)
#' and is covered at validation scale. The covered common-environment two-effect, independent multi-effect,
#' direct-maternal, and k = 2 random-regression models remain opt-in. Repeatability,
#' maternal two-effect, and non-Gaussian (`poisson(log)`/`binomial(logit)`,
#' Laplace marginal likelihood or variational ELBO) models are opt-in and
#' experimental; single-step stays **opt-in partial**. Factor-analytic G is
#' engine-covered on the Julia engine (`V4-FA`) and planned, not fitted, on the
#' R formula (Julia engine-covered is not R covered).
#'
#' @section Current limitations:
#' This package is an **experimental 0.9.0 release**; it is not production or
#' 1.0.0, and CRAN availability is tracked separately. Default [hsquared()] fitting requires a local Julia installation
#' and HSquared.jl; use [hs_control()] with `engine = "validate"` to check the
#' model contract without fitting. Report point estimates only for `covered`
#' routes; uncertainty intervals are experimental and not coverage-calibrated.
#'
#' What you may report is listed on
#' [Can I fit and report this?](
#' https://itchyshin.github.io/hsquared/articles/current-limits.html)
#' -- not in [validation_status()]. That article lists every covered route
#' and its reporting scope. `validation_status()` is a developer evidence
#' table of validation atoms; it is not the user-facing list. The
#' [model status](https://itchyshin.github.io/hsquared/articles/model-status.html)
#' and [validation evidence](
#' https://itchyshin.github.io/hsquared/articles/validation-evidence.html)
#' articles are published on the package website rather than as installed
#' vignettes.
#'
#' @keywords internal
## usethis namespace: start
#' @importFrom lifecycle deprecated
## usethis namespace: end
"_PACKAGE"
