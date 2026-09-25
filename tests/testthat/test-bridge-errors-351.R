# Post-fit engine failures must surface, never vanish (hsquared / HSquared.jl#351).
#
# The bridge computes SEs, intervals, and plot data inside Julia `try` blocks so
# a failure there never aborts the fit. Before this slice the `catch` swallowed
# the error: a collaborator's real fit came back with SE = NA for every variance
# component and no message, so nobody could tell "undefined here" from "the
# engine threw" (finding: Szymon Drobniak, 2026-09-17). Now each `catch` records
# `sprint(showerror, err)` under the engine function's name; after the fit the
# bridge attaches the record as `attr(fit, "bridge_errors")` (a named character
# vector) and raises ONE consolidated warning naming what failed.
#
# The pure-R tests need no Julia. The live tests force a throw by rewriting the
# engine call inside the Julia command string at the JuliaCall boundary, which
# is the only lever that reproduces the defect against unmodified engine code.

# --- pure R -------------------------------------------------------------------

test_that("hs_format_bridge_errors() names every failed function in one message", {
  errors <- c(
    variance_component_standard_errors = "SingularException(2)",
    breeding_values_plot_data = "OutOfMemoryError()"
  )
  msg <- hsquared:::hs_format_bridge_errors(errors)
  expect_type(msg, "character")
  expect_length(msg, 1L)
  expect_match(msg, "variance_component_standard_errors()", fixed = TRUE)
  expect_match(msg, "SingularException(2)", fixed = TRUE)
  expect_match(msg, "breeding_values_plot_data()", fixed = TRUE)
  expect_match(msg, "OutOfMemoryError()", fixed = TRUE)
  expect_match(msg, "attr(fit, \"bridge_errors\")", fixed = TRUE)
})

test_that("hs_format_bridge_errors() shows at most 200 characters per message", {
  long <- strrep("x", 500L)
  msg <- hsquared:::hs_format_bridge_errors(c(heritability_interval = long))
  expect_match(msg, strrep("x", 200L), fixed = TRUE)
  expect_false(grepl(strrep("x", 201L), msg, fixed = TRUE))
})

test_that("hs_attach_bridge_errors() is a no-op when nothing failed", {
  fit <- hsquared:::hs_new_fit(
    spec = list(method = "REML", family = list(family = "gaussian")),
    payload = list(y = seq_len(4)),
    result = list(heritability = data.frame(term = "animal", estimate = 0.4))
  )
  expect_no_warning(out <- hsquared:::hs_attach_bridge_errors(fit, character()))
  expect_identical(out, fit)
  expect_null(attr(out, "bridge_errors"))
})

test_that("hs_attach_bridge_errors() attaches the record and warns exactly once", {
  fit <- hsquared:::hs_new_fit(
    spec = list(method = "REML", family = list(family = "gaussian")),
    payload = list(y = seq_len(4)),
    result = list(heritability = data.frame(term = "animal", estimate = 0.4))
  )
  errors <- c(
    variance_component_standard_errors = "SingularException(2)",
    heritability_standard_error = "SingularException(2)"
  )
  warnings <- character()
  out <- withCallingHandlers(
    hsquared:::hs_attach_bridge_errors(fit, errors),
    warning = function(cnd) {
      warnings <<- c(warnings, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    }
  )
  expect_length(warnings, 1L)
  expect_match(warnings, "variance_component_standard_errors()", fixed = TRUE)
  expect_match(warnings, "heritability_standard_error()", fixed = TRUE)
  expect_s3_class(out, "hsquared_fit")
  expect_identical(attr(out, "bridge_errors"), errors)
})

test_that("hs_julia_try_slot() records the error text instead of swallowing it", {
  cmd <- hsquared:::hs_julia_try_slot(
    "hsq_vcse",
    "HSquared.variance_component_standard_errors(hsq_fit)",
    "variance_component_standard_errors"
  )
  expect_match(cmd, "^hsq_vcse = try; ")
  expect_match(cmd, "catch err;", fixed = TRUE)
  expect_match(
    cmd,
    "hsq_bridge_errors[\"variance_component_standard_errors\"] = sprint(showerror, err); nothing; end;",
    fixed = TRUE
  )
  expect_false(grepl("catch; nothing", cmd, fixed = TRUE))
})

# --- live (skip-guarded on the local HSquared.jl bridge) ----------------------

hs_bridge_errors_351_skip <- function() {
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for this live test."
  )
}

# Rewrite `from` -> `to` inside every Julia command the bridge sends, for the
# duration of the calling test. `from` is the engine call text as the bridge
# writes it; `to` is a Julia expression that throws.
hs_local_julia_fault <- function(from, to, env = parent.frame()) {
  real <- JuliaCall::julia_command
  testthat::local_mocked_bindings(
    julia_command = function(cmd, ...) {
      real(gsub(from, to, cmd, fixed = TRUE), ...)
    },
    .package = "JuliaCall",
    .env = env
  )
}

hs_collect_bridge_warnings <- function(expr) {
  warnings <- character()
  value <- withCallingHandlers(
    expr,
    warning = function(cnd) {
      warnings <<- c(warnings, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = warnings[grepl("bridge_errors", warnings, fixed = TRUE)])
}

test_that("a throwing SE call surfaces one warning and attr(fit, 'bridge_errors') [live]", {
  hs_bridge_errors_351_skip()
  hs_local_julia_fault(
    "HSquared.variance_component_standard_errors(hsq_fit)",
    "error(\"forced bridge fault 351\")"
  )
  fx <- hsquared:::hs_mrode_supplied_variance_validation_fixture()

  got <- hs_collect_bridge_warnings(
    hsquared(fx$formula, data = fx$data, family = stats::gaussian(), REML = TRUE)
  )
  fit <- got$value

  # (a) the fit still returns, on the default ai_reml route
  expect_s3_class(fit, "hsquared_fit")
  expect_equal(fit$spec$target, "ai_reml")
  # (b) exactly one warning, naming the function that threw
  expect_length(got$warnings, 1L)
  expect_match(got$warnings, "variance_component_standard_errors()", fixed = TRUE)
  expect_match(got$warnings, "forced bridge fault 351", fixed = TRUE)
  # (c) the message is recoverable from the fit
  errors <- attr(fit, "bridge_errors")
  expect_type(errors, "character")
  expect_named(errors, "variance_component_standard_errors")
  expect_match(errors[["variance_component_standard_errors"]], "forced bridge fault 351", fixed = TRUE)
  # the failed quantity is absent, not a silent NA; its siblings are untouched
  expect_null(fit$result$variance_component_se)
  expect_error(variance_component_standard_errors(fit), "does not contain")
  expect_true(is.finite(heritability_standard_error(fit)$se))
  expect_s3_class(heritability_interval(fit), "data.frame")
})

test_that("a throwing plot-data call (OutOfMemoryError) is reported, not dropped [live]", {
  hs_bridge_errors_351_skip()
  hs_local_julia_fault(
    "HSquared.breeding_values_plot_data(hsq_fit)",
    "throw(OutOfMemoryError())"
  )
  fx <- hsquared:::hs_mrode_supplied_variance_validation_fixture()

  got <- hs_collect_bridge_warnings(
    hsquared(fx$formula, data = fx$data, family = stats::gaussian(), REML = TRUE)
  )
  fit <- got$value

  expect_s3_class(fit, "hsquared_fit")
  expect_length(got$warnings, 1L)
  expect_match(got$warnings, "breeding_values_plot_data()", fixed = TRUE)
  expect_match(got$warnings, "OutOfMemoryError", fixed = TRUE)
  expect_named(attr(fit, "bridge_errors"), "breeding_values_plot_data")
  expect_null(fit$result$breeding_values_plot_data)
  expect_false(is.null(fit$result$variance_components_plot_data))
  # the SEs were unaffected by the plot-data failure
  expect_s3_class(variance_component_standard_errors(fit), "data.frame")
})

test_that("two post-fit failures are consolidated into ONE warning [live]", {
  hs_bridge_errors_351_skip()
  hs_local_julia_fault(
    "HSquared.variance_component_standard_errors(hsq_fit)",
    "error(\"fault A\")"
  )
  hs_local_julia_fault(
    "HSquared.heritability_standard_error(hsq_fit)",
    "error(\"fault B\")"
  )
  fx <- hsquared:::hs_mrode_supplied_variance_validation_fixture()

  got <- hs_collect_bridge_warnings(
    hsquared(fx$formula, data = fx$data, family = stats::gaussian(), REML = TRUE)
  )

  expect_length(got$warnings, 1L)
  expect_match(got$warnings, "fault A", fixed = TRUE)
  expect_match(got$warnings, "fault B", fixed = TRUE)
  expect_setequal(
    names(attr(got$value, "bridge_errors")),
    c("variance_component_standard_errors", "heritability_standard_error")
  )
})

test_that("a clean fit carries no bridge_errors attribute and no warning [live]", {
  hs_bridge_errors_351_skip()
  fx <- hsquared:::hs_mrode_supplied_variance_validation_fixture()
  got <- hs_collect_bridge_warnings(
    hsquared(fx$formula, data = fx$data, family = stats::gaussian(), REML = TRUE)
  )
  expect_length(got$warnings, 0L)
  expect_null(attr(got$value, "bridge_errors"))
  expect_s3_class(variance_component_standard_errors(got$value), "data.frame")
})

test_that("a throwing multivariate covariance SE call surfaces on the cbind route [live]", {
  hs_bridge_errors_351_skip()
  hs_local_julia_fault(
    "HSquared.multivariate_covariance_standard_errors(",
    "error(\"forced multivariate fault 351\"); ("
  )
  ped <- data.frame(
    id = c("s1", "d1", "s2", "d2", "a", "b", "c", "d"),
    sire = c(NA, NA, NA, NA, "s1", "s1", "s2", "s2"),
    dam = c(NA, NA, NA, NA, "d1", "d1", "d2", "d2")
  )
  dat <- data.frame(
    y1 = c(1.0, 1.8, 1.2, 2.0, 3.0, 3.4, 2.8, 3.2),
    y2 = c(2.1, 1.5, 2.2, 1.7, 3.1, NA, 3.0, 2.8),
    id = ped$id
  )

  got <- hs_collect_bridge_warnings(
    hsquared(
      cbind(y1, y2) ~ animal(1 | id, pedigree = ped),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE,
      control = hs_control(
        engine = "julia",
        engine_control = list(target = "multivariate", iterations = 400L)
      )
    )
  )
  fit <- got$value

  expect_s3_class(fit, "hsquared_fit")
  expect_equal(fit$spec$target, "multivariate")
  expect_length(got$warnings, 1L)
  expect_match(got$warnings, "multivariate_covariance_standard_errors()", fixed = TRUE)
  expect_match(got$warnings, "forced multivariate fault 351", fixed = TRUE)
  expect_named(attr(fit, "bridge_errors"), "multivariate_covariance_standard_errors")
  expect_null(fit$result$covariance_standard_errors)
  expect_error(covariance_standard_errors(fit), "does not contain")
})
