# CRAN-safe skip for any test that may reach JuliaCall::julia_setup().
# Pattern stolen from drmTMB `drm_skip_live_julia()` after Ligges hung ~10k s
# inside julia_setup() on a CRAN-lane expect_error(engine = "julia") path.
#
# Opt in with HSQUARED_JULIA_TESTS=true. Repository CI / maintainer live suites
# use NOT_CRAN=true (this skip is then a no-op aside from skip_on_cran).
#
# Call this at the top of every live Julia test (or shared skip helper) before
# hs_julia_setup() / engine = "julia" fits. Intentional pre-setup gates
# (hsquared_unsupported_syntax) may still run on CRAN without this helper.

hs_skip_live_julia <- function() {
  if (identical(Sys.getenv("HSQUARED_JULIA_TESTS"), "true")) {
    return(invisible(TRUE))
  }
  if (
    !isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false"))) &&
      !interactive()
  ) {
    testthat::skip(
      "Live Julia skipped on CRAN lane (set NOT_CRAN=true or HSQUARED_JULIA_TESTS=true)."
    )
  }
  testthat::skip_on_cran()
  invisible(TRUE)
}

# Require a live Julia bridge on the Tier-1 / maintainer path that sets
# HSQUARED_REQUIRE_BRIDGE=true. Without that variable, missing JuliaCall /
# julia / Project.toml remains a quiet skip (CRAN and ordinary local runs).
# With it set, a failed provision must stop() — never skip-as-pass — so a
# future Tier-1 job cannot go green with the A26 parity legs absent.
#
# Mirror of hs_require_suggests() (A26b). Call from A26 live helpers (and any
# future Tier-1-gated bridge file) instead of skip_if_not(hs_julia_bridge_available()).

hs_require_bridge <- function(gate = NULL, project = NULL) {
  hs_skip_live_julia()
  project <- if (is.null(project)) {
    hsquared:::hs_default_julia_project()
  } else {
    project
  }
  if (isTRUE(hsquared:::hs_julia_bridge_available(project))) {
    return(invisible(TRUE))
  }
  if (!isTRUE(as.logical(Sys.getenv("HSQUARED_REQUIRE_BRIDGE", "false")))) {
    testthat::skip(
      "JuliaCall, Julia, and a local HSquared.jl project are required."
    )
  }
  label <- if (is.null(gate) || !nzchar(gate)) {
    "Live Julia bridge"
  } else {
    gate
  }
  stop(
    label,
    " requires JuliaCall + julia on PATH + an HSquared.jl Project.toml, and ",
    "HSQUARED_REQUIRE_BRIDGE=true was set. Fix the provisioning; do not ",
    "unset the variable to get a green run.",
    call. = FALSE
  )
}
