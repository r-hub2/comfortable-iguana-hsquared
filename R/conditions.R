# Structured conditions for hsquared errors.

# Raise a structured "unsupported syntax" error.
#
# Genuine grammar/target rejections -- planned, opt-in-required, or
# not-implemented formula terms, unsupported response grammar, and unsupported
# `engine_control$target` branches -- are raised through this helper so they
# carry a stable, catchable condition class:
#
#   c("hsquared_unsupported_syntax", "hsquared_error", "error", "condition")
#
# Data validation, wrong-object, dimension, numeric-invariant, and internal
# payload guards stay plain `stop()`; they are a different family and must not
# carry this class.
#
# The `...` message pieces are concatenated with `paste0()`, matching the
# existing `stop("a", var, "b")` style. `call. = FALSE` is accepted and
# discarded so call sites can keep their existing ergonomics: the raised
# condition never carries a call (mirroring the previous `call. = FALSE`).
hs_abort_unsupported_syntax <- function(..., call. = FALSE) {
  cond <- errorCondition(
    paste0(...),
    class = c("hsquared_unsupported_syntax", "hsquared_error")
  )
  stop(cond)
}

# One-line hints appended by `hs_julia_fit()` (hsquared#214, #217). A hint may
# name only a lever the failing target actually accepts: `max_dense_cells` at
# the two routes that forward it to a guarded engine entry point,
# `ridge`/`blend_weight` at the single-step H^-1 construction path, and a
# lever-free scale note everywhere else.

# Appended only at the two call sites that actually forward `max_dense_cells`
# to a guarded engine entry point (`fit_variance_components`,
# `fit_repeatability_reml`). Naming the lever anywhere else would point the
# user at a control that target does not accept (hsquared#212's defect class).
hs_dense_route_hint <- "raise max_dense_cells in engine_control or use a sparse route"

# The repeatability target's own hint. Unlike the generic one above, this route
# HAS a named sparse escape, so the message says which control to set instead of
# leaving "use a sparse route" as an instruction the reader cannot act on. The
# dense ceiling scales as nobs^2 + nanimals^2, so any real repeated-measures
# pedigree hits it long before the model is the problem.
hs_repeatability_dense_route_hint <- paste0(
  "set engine_control = list(target = \"repeatability\", ",
  "scale_method = \"auto\") to fit the same model through the sparse route, ",
  "or raise max_dense_cells in engine_control"
)

# Appended at the remaining dense/validation-scale routes, which enforce no
# cell cap and expose no lever to raise.
hs_dense_scale_hint <- "this route is dense and validation-scale (n <= ~1000); reduce the problem size or use a sparse target"

hs_single_step_ridge_hint <- "raise ridge or blend_weight in engine_control"

# Evaluate `expr` (a Julia-bridge fit call) and translate ANY error it raises
# into a structured condition of class:
#
#   c("hsquared_julia_error", "hsquared_error", "error", "condition")
#
# instead of letting a raw, untranslated Julia trace cross the bridge
# (hsquared#214, #217). `hint`, when supplied, is appended on its own line so
# the R-facing lever (e.g. `ridge`, `max_dense_cells`) is named at the point of
# failure; the original message is always preserved. A non-error result of
# `expr` passes through unchanged. Mirrors `hs_abort_unsupported_syntax()`: the
# raised condition never carries a call.
hs_julia_fit <- function(expr, hint = NULL) {
  tryCatch(
    expr,
    error = function(e) {
      msg <- conditionMessage(e)
      if (!is.null(hint) && nzchar(hint)) {
        msg <- paste0(msg, "\n", hint)
      }
      cond <- errorCondition(
        msg,
        class = c("hsquared_julia_error", "hsquared_error")
      )
      stop(cond)
    }
  )
}

# Post-fit quantities (SEs, intervals, plot data) are computed inside Julia
# `try` blocks so a failure there never aborts the fit. The bridge records each
# failure as `sprint(showerror, err)` under the engine function's name
# (hsquared / HSquared.jl#351; before that the `catch` swallowed it, and a
# missing SE was indistinguishable from an undefined one).
# `hs_attach_bridge_errors()` stores the record as `attr(fit, "bridge_errors")`
# (a named character vector) and raises ONE consolidated warning naming what
# failed; the warning shows the first 200 characters of each message, the
# attribute keeps them whole.
hs_format_bridge_errors <- function(errors) {
  n <- length(errors)
  shown <- substr(unname(errors), 1L, 200L)
  paste0(
    "The fit succeeded, but ",
    n,
    " post-fit ",
    ngettext(n, "quantity", "quantities"),
    " could not be computed and ",
    ngettext(n, "is", "are"),
    " absent from the fit. The engine threw:\n",
    paste0("  - ", names(errors), "(): ", shown, collapse = "\n"),
    "\nFull messages: attr(fit, \"bridge_errors\")."
  )
}

hs_attach_bridge_errors <- function(fit, errors) {
  errors <- errors[nzchar(errors)]
  if (length(errors) == 0L) {
    return(fit)
  }
  attr(fit, "bridge_errors") <- errors
  warning(hs_format_bridge_errors(errors), call. = FALSE)
  fit
}

# Raise a structured "out of range" error for a numeric value outside its
# valid domain -- the normalized-Legendre basis input `t`, a user-supplied
# `at` outside the fitted covariate range, or a `markers` dosage matrix
# outside its required interval. Carries a stable, catchable condition class:
#
#   c("hsquared_out_of_range", "hsquared_error", "error", "condition")
hs_abort_out_of_range <- function(...) {
  cond <- errorCondition(
    paste0(...),
    class = c("hsquared_out_of_range", "hsquared_error")
  )
  stop(cond)
}

# Raise a structured "boundary refused" error: an R-side defense-in-depth
# guard for the non-Gaussian three-field route (hsquared#222). The Julia
# payload builder (`nongaussian_three_field_payload`, HSquared.jl#342) already
# refuses to hand back a `boundary = true` fit, so this guard fires only if a
# FUTURE bridge route constructs the `nongaussian_three_field_v09` envelope
# without going through that refusing builder. Carries a stable, catchable
# condition class:
#
#   c("hsquared_boundary_refused", "hsquared_error", "error", "condition")
hs_abort_boundary_refused <- function(...) {
  cond <- errorCondition(
    paste0(...),
    class = c("hsquared_boundary_refused", "hsquared_error")
  )
  stop(cond)
}

# Session-scoped flags for warn-once messages. Reset in tests.
hs_session_flags <- new.env(parent = emptyenv())

hs_reset_session_flags <- function() {
  rm(list = ls(envir = hs_session_flags), envir = hs_session_flags)
  invisible(NULL)
}

# Default-path cbind() routes to the multivariate fitter (MV-4). G10 covered
# t=2 unstructured at validation scale; the experimental label is retained.
# Warn once per session so easy syntax does not look like interval-calibrated
# or k>=3 / diagonal coverage.
hs_warn_cbind_experimental_once <- function() {
  if (isTRUE(hs_session_flags$cbind_experimental)) {
    return(invisible(FALSE))
  }
  hs_session_flags$cbind_experimental <- TRUE
  warning(
    "This cbind() model fitted; multivariate is covered at validation scale (experimental).\n",
    "Report point estimates for t=2 unstructured G0/R0 only; intervals are not coverage-calibrated.",
    call. = FALSE
  )
  invisible(TRUE)
}

# Deparse a user formula or data argument into a pasteable snippet.
hs_deparse_user_expr <- function(x) {
  paste(deparse(x, width.cutoff = 500L), collapse = "\n    ")
}

# Pasteable next call for an opt-in engine target, using the formula the
# user already wrote.
hs_format_next_call <- function(formula, data_name, target) {
  paste0(
    "  hsquared(\n",
    "    ",
    hs_deparse_user_expr(formula),
    ",\n",
    "    data = ",
    data_name,
    ",\n",
    "    control = hs_control(\n",
    "      engine = \"julia\",\n",
    "      engine_control = list(target = \"",
    target,
    "\")\n",
    "    )\n",
    "  )"
  )
}

hs_opt_in_term_label <- function(type) {
  switch(
    type,
    common_env = "common_env()",
    permanent = "permanent()",
    maternal_genetic = "maternal_genetic()",
    iid_effects = "(1 | group)",
    random_regression = "rr(...)",
    paste0(type, "()")
  )
}

# Honest one-paragraph status for the pasteable next call. Status words
# match capability-status: no covered flip, no new estimand.
hs_opt_in_route_note <- function(type) {
  switch(
    type,
    common_env = paste0(
      "That route is covered for point estimates only (validation-scale, ",
      "common-environment leg). Intervals are experimental."
    ),
    permanent = paste0(
      "That route is experimental. Do not report these numbers from ",
      "hsquared alone."
    ),
    maternal_genetic = paste0(
      "`target = \"two_effect\"` is experimental (independent maternal). ",
      "`target = \"direct_maternal\"` is covered for point estimates only ",
      "(validation-scale Willham triple: h2_direct, m2, r_am -- not a ",
      "scalar h2). Intervals are experimental."
    ),
    random_regression = paste0(
      "That route is covered for point estimates only (validation-scale, ",
      "k = 2). Higher order is experimental."
    ),
    paste0(
      "That route is experimental and opt-in. See formula_status()."
    )
  )
}

# Default-path abort that prints the call the user should paste.
# Used for common_env / permanent / maternal_genetic / rr(...).
hs_abort_opt_in_next_call <- function(
  type,
  formula,
  data_name,
  target = NULL
) {
  if (identical(type, "random_regression")) {
    target <- "random_regression"
  } else if (is.null(target)) {
    target <- hs_second_effect_target(type)
  }
  label <- hs_opt_in_term_label(type)
  extra <- ""
  if (identical(type, "maternal_genetic")) {
    extra <- paste0(
      "\n\nCovered alternative (Willham triple, not a scalar h2):\n\n",
      hs_format_next_call(formula, data_name, "direct_maternal")
    )
  }
  hs_abort_unsupported_syntax(
    "`",
    label,
    "` is not on the default path.\n\n",
    "Closest working call:\n\n",
    hs_format_next_call(formula, data_name, target),
    extra,
    "\n\n",
    hs_opt_in_route_note(type),
    "\nSee: current-limits article, formula_status()."
  )
}


# Extra clause when engine = "julia" has a single_step() formula but the wrong
# target (or the default fit_animal_model). Names the opt-in fence explicitly.
hs_single_step_allowed_targets_note <- function() {
  paste0(
    " Ordinary/default single-step stays held. Julia `V2-SSHINV` is ",
    "engine-covered; that is not an R-public ordinary-path promotion. ",
    "public_covered_count stays 7."
  )
}

# Default-path copy for single_step(). The supplied-Hinv and engine-built
# construction routes are distinct, so a generic target-only suggestion would
# leave users unsure which formula arguments they need.
hs_abort_single_step_default_path <- function() {
  hs_abort_unsupported_syntax(
    "`single_step()` is not on the default `engine = \"fit\"` path. ",
    "R single-step remains experimental and opt-in. Ordinary/default ",
    "single-step stays held. Julia `V2-SSHINV` is engine-covered ",
    "(HSquared.jl cf2a9bbf / #301); that is not an R-public ordinary-path ",
    "promotion. public_covered_count stays 7.\n\n",
    "For a supplied inverse, use `single_step(1 | id, Hinv = Hinv)` with ",
    "`control = hs_control(engine = \"julia\", engine_control = list(",
    "target = \"single_step\"))`.\n\n",
    "For engine-built H^-1, use ",
    "`single_step(1 | id, pedigree = ped, markers = M)` with ",
    "`control = hs_control(engine = \"julia\", engine_control = list(",
    "target = \"single_step_construct\"))` (or the matching `hs_data()` ",
    "bundle). Neither route is a covered default R formula path.\n\n",
    "See: current-limits article, formula_status().",
    call. = FALSE
  )
}

# Default-path copy for genomic(). Covered GREML stays explicit
# target = "genomic"; ordinary-call / default-route activation remains held
# (design-44 G5; sealed BOUNDARY_HOLDOUT_FAIL). Do not silent-route.
hs_abort_genomic_default_path <- function(formula, data_name) {
  hs_abort_unsupported_syntax(
    "`genomic()` is not on the default `engine = \"fit\"` path.\n\n",
    "Ordinary-call / default-route genomic activation remains held ",
    "(design-44 G5; sealed BOUNDARY_HOLDOUT_FAIL). Covered genomic GREML ",
    "stays explicit opt-in only; this is not a silent ordinary-route ",
    "promotion and public_covered_count stays 7.\n\n",
    "Closest working call (validation-scale GREML; report ",
    "`genomic_variance_ratio`, not a bare pedigree h2):\n\n",
    hs_format_next_call(formula, data_name, "genomic"),
    "\n\n",
    "For marker-effect SNP-BLUP use the same `genomic(1 | id, markers = M)` ",
    "formula with `target = \"snp_blup\"` instead. Intervals and SE accessors ",
    "remain unavailable on the genomic GREML route.\n\n",
    "See: current-limits article, formula_status().",
    call. = FALSE
  )
}

# One rule for the public hsquared() door: ML is not a live path.
# REML = TRUE is the nearest working call. Internal spec builders and
# model_spec() may still construct an ML-labelled spec; engine = "julia"
# keeps its supplied-variance exemptions in R/hsquared.R.
hs_abort_reml_false <- function() {
  hs_abort_unsupported_syntax(
    "ML estimation (`REML = FALSE`) is not implemented.\n\n",
    "Closest working call: `REML = TRUE` (the default). ",
    "The live path estimates variance components by REML ",
    "(average-information REML).\n\n",
    "See: current-limits article, formula_status()."
  )
}

hs_public_hsquared_control <- function() {
  calls <- sys.calls()
  frames <- sys.frames()
  n <- min(length(calls), length(frames))
  if (n < 1L) {
    return(NULL)
  }
  for (i in seq_len(n)) {
    fn <- calls[[i]][[1L]]
    if (!identical(fn, quote(hsquared))) {
      next
    }
    env <- frames[[i]]
    if (!exists("control", envir = env, inherits = FALSE)) {
      next
    }
    cand <- get("control", envir = env, inherits = FALSE)
    if (inherits(cand, "hs_control")) {
      return(cand)
    }
  }
  NULL
}

hs_abort_reml_false_on_public_path <- function(REML) {
  if (isTRUE(REML)) {
    return(invisible(FALSE))
  }
  control <- hs_public_hsquared_control()
  if (is.null(control)) {
    return(invisible(FALSE))
  }
  if (
    identical(control$engine, "validate") || identical(control$engine, "fit")
  ) {
    hs_abort_reml_false()
  }
  invisible(FALSE)
}

# Repeated records per individual with NO permanent-environment term
# (HSquared.jl#352 item 3).
#
# With more than one record per animal and only the additive effect in the
# model, sigma^2_a is not narrow-sense additive variance: it absorbs the
# permanent-environment variance too, because nothing else in the model can
# explain why repeated records on the same individual resemble each other. The
# estimate is silently inflated and h^2 with it -- on the great tit data this
# returned a single animal = 1.110 where ASReml's vm(animal) + ide(animal) was
# 0.595 + 0.525.
#
# The upstream issue asked the R side to REFUSE this design. It warns instead:
# refusing would reject models that are legitimately specified this way (a
# deliberately repeatability-free analysis, or a design where the second
# component is known to be unidentifiable), and at the time the issue was filed
# there was no reachable alternative to point at. There is now, so the warning
# names it.
hs_warn_unmodelled_repeated_records <- function(spec) {
  animal <- spec$random$animal
  if (is.null(animal) || !is.null(spec$random$permanent)) {
    return(invisible(FALSE))
  }
  # Any other id-level random effect can absorb the between-record correlation
  # too, so only warn when the animal effect is genuinely alone on these ids.
  others <- spec$random[!names(spec$random) %in% "animal"]
  if (length(others) > 0L) {
    return(invisible(FALSE))
  }
  values <- animal$values
  if (is.null(values)) {
    return(invisible(FALSE))
  }
  values <- as.character(values)
  n_records <- length(values)
  n_ids <- length(unique(values))
  if (n_records <= n_ids) {
    return(invisible(FALSE))
  }
  preamble <- paste0(
    "The data have repeated records per individual (",
    n_records,
    " records for ",
    n_ids,
    " individuals), but the model has no ",
    "`permanent(1 | ...)` term. The additive variance will ABSORB the ",
    "permanent-environment variance, so `animal` and the heritability derived ",
    "from it are inflated -- they are not narrow-sense quantities here.\n"
  )
  # A multivariate `cbind()` response reaches an animal-only fitter: the spec
  # fence admits a single random-intercept animal effect and nothing else, so
  # `permanent()` is NOT accepted on this call (hsquared#237). Naming
  # `permanent()` *on the cbind formula* would be hsquared#212's defect class.
  # Point instead at separate live routes (univariate repeatability, or
  # cbind without PE) and keep the engine-debt note explicit.
  if (isTRUE(spec$response$multivariate)) {
    warning(
      preamble,
      "A multivariate `cbind()` response cannot carry a `permanent()` term ",
      "on the current route (hsquared#237): the engine has no multi-trait PE ",
      "fitter yet (G0 x A + P0 x I). Closest live paths: (1) fit each trait ",
      "univariately with `animal(...) + permanent(1 | <id>)` and ",
      "`control = hs_control(engine = \"julia\", engine_control = list(",
      "target = \"repeatability\"))`; (2) keep `cbind(...)` without ",
      "`permanent()` and treat animal/G0 as absorbing PE (not narrow-sense). ",
      "Full MV+PE remains deferred; public_covered_count stays 7. See ",
      "docs/design/57-mv-pe-cbind-permanent-237.md.\n",
      "Suppress with suppressWarnings() if this is intended.",
      call. = FALSE
    )
    return(invisible(TRUE))
  }
  warning(
    preamble,
    "To separate them, add a permanent-environment effect:\n",
    "  hsquared(\n",
    "    <response> ~ <fixed> + animal(1 | <id>, pedigree = ped) + ",
    "permanent(1 | <id>),\n",
    "    data = <data>,\n",
    "    control = hs_control(engine = \"julia\", engine_control = list(\n",
    "      target = \"repeatability\", scale_method = \"auto\"))\n",
    "  )\n",
    "Suppress with suppressWarnings() if the single-effect model is intended.",
    call. = FALSE
  )
  invisible(TRUE)
}
