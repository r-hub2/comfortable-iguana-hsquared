# Require a DESCRIPTION Suggests package on the maintainer / CI path
# (NOT_CRAN=true). On CRAN (or any lane without NOT_CRAN=true) fall back to a
# quiet skip so Suggests stay optional for CRAN. Headline flip gates (MV-1) must
# not become silent skips when Suggests fail to install under CI.
#
# Call after skip_on_cran() for CI-gated comparator legs. Prefer this over a
# bare skip_if_not_installed() for any gate whose absence would hide a red build.

hs_require_suggests <- function(pkg, gate = NULL) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    return(invisible(TRUE))
  }
  if (!isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))) {
    testthat::skip_if_not_installed(pkg)
  }
  label <- if (is.null(gate) || !nzchar(gate)) {
    sprintf("Suggests package '%s'", pkg)
  } else {
    sprintf("%s requires Suggests package '%s'", gate, pkg)
  }
  stop(
    label,
    " when NOT_CRAN=true (CI/maintainer path). Install Suggests, or unset ",
    "NOT_CRAN to allow a CRAN-lane skip.",
    call. = FALSE
  )
}
