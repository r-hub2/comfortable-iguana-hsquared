# Opt-in, experimental common-environment model (animal additive genetic + an
# IID environmental effect). Surfaces the twin's fit_two_effect_reml behind the
# explicit `engine = "julia", target = "two_effect"` path. REML only.

test_that("the parser accepts common_env() alongside animal()", {
  ped <- data.frame(
    id = c("a", "b", "c", "d"),
    sire = c(NA, NA, "a", "a"),
    dam = c(NA, NA, "b", "b")
  )
  dat <- data.frame(
    y = c(1, 2, 3, 4),
    id = c("a", "b", "c", "d"),
    litter = c("l1", "l1", "l2", "l2")
  )

  spec <- hsquared:::hs_build_model_spec(
    y ~ animal(1 | id, pedigree = ped) + common_env(1 | litter),
    data = dat,
    family = stats::gaussian(),
    REML = TRUE
  )

  expect_named(spec$random, c("animal", "common_env"))
  expect_equal(spec$random$common_env$type, "common_env")
  expect_equal(spec$random$common_env$group, "litter")
  expect_equal(spec$random$common_env$design, "intercept")
  expect_equal(spec$random$common_env$relationship, "identity")
  expect_equal(spec$random$common_env$levels, c("l1", "l2"))
  expect_match(spec$bridge$target, "two_effect", fixed = TRUE)
})

test_that("at most one second random effect is allowed", {
  ped <- data.frame(id = c("a", "b"), sire = c(NA, NA), dam = c(NA, NA))
  dat <- data.frame(
    y = c(1, 2, 1.5, 2.5),
    id = c("a", "b", "a", "b"),
    litter = c("l1", "l1", "l2", "l2")
  )
  expect_error(
    hsquared:::hs_build_model_spec(
      y ~ animal(1 | id, pedigree = ped) +
        permanent(1 | id) +
        common_env(1 | litter),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "at most one additional random effect",
    fixed = TRUE
  )
})

test_that("common_env() must be intercept-only and a known column", {
  ped <- data.frame(id = c("a", "b"), sire = c(NA, NA), dam = c(NA, NA))
  dat <- data.frame(y = c(1, 2), id = c("a", "b"), litter = c("l1", "l2"))
  expect_error(
    hsquared:::hs_build_model_spec(
      y ~ animal(1 | id, pedigree = ped) + common_env(age | litter),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "random-intercept",
    fixed = TRUE
  )
  expect_error(
    hsquared:::hs_build_model_spec(
      y ~ animal(1 | id, pedigree = ped) + common_env(1 | nope),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "was not found in `data`",
    fixed = TRUE
  )
})

test_that("two_effect is a valid opt-in julia target", {
  expect_equal(
    hsquared:::hs_validate_julia_target("two_effect"),
    "two_effect"
  )
})

test_that("the default engine = \"fit\" rejects a common_env() formula", {
  ped <- data.frame(id = c("a", "b"), sire = c(NA, NA), dam = c(NA, NA))
  dat <- data.frame(
    y = c(1, 2, 3, 4),
    id = c("a", "b", "a", "b"),
    litter = c("l1", "l1", "l2", "l2")
  )
  expect_error(
    hsquared(
      y ~ animal(1 | id, pedigree = ped) + common_env(1 | litter),
      data = dat,
      family = stats::gaussian()
    ),
    "Closest working call",
    fixed = TRUE
  )
  err <- tryCatch(
    hsquared(
      y ~ animal(1 | id, pedigree = ped) + common_env(1 | litter),
      data = dat,
      family = stats::gaussian()
    ),
    error = function(e) conditionMessage(e)
  )
  expect_true(grepl(
    "`common_env()` is not on the default path",
    err,
    fixed = TRUE
  ))
  expect_true(grepl("target = \"two_effect\"", err, fixed = TRUE))
  expect_true(grepl(
    "y ~ animal(1 | id, pedigree = ped) + common_env(1 | litter)",
    err,
    fixed = TRUE
  ))
  expect_true(grepl("data = dat", err, fixed = TRUE))
  expect_true(grepl("covered for point estimates only", err, fixed = TRUE))
})

test_that("target = \"two_effect\" requires a common_env() term", {
  ped <- data.frame(id = c("a", "b"), sire = c(NA, NA), dam = c(NA, NA))
  dat <- data.frame(y = c(1, 2), id = c("a", "b"))
  expect_error(
    hsquared(
      y ~ animal(1 | id, pedigree = ped),
      data = dat,
      family = stats::gaussian(),
      control = hs_control(
        engine = "julia",
        engine_control = list(target = "two_effect")
      )
    ),
    "requires a `common_env(1 | group)` term",
    fixed = TRUE
  )
})

test_that("the two-effect bridge requires an internal payload", {
  expect_error(
    hsquared:::hs_fit_julia_two_effect_payload(list()),
    "`payload` must be an internal `hs_bridge_payload`.",
    fixed = TRUE
  )
})

test_that("common_env_effects extractor requires a fitted object", {
  expect_error(
    common_env_effects(1),
    "requires an `hsquared_fit`",
    fixed = TRUE
  )
})

test_that("common_env_proportion / interval extractors require a fitted object", {
  expect_error(
    common_env_proportion(1),
    "requires an `hsquared_fit`",
    fixed = TRUE
  )
  expect_error(
    common_env_proportion_interval(1),
    "requires an `hsquared_fit`",
    fixed = TRUE
  )
})

test_that("common_env_proportion() returns c2 with a Falconer interpretation", {
  result <- list(
    common_env_proportion = data.frame(
      term = "common_env",
      estimate = 0.23
    )
  )
  fit <- hsquared:::hs_new_fit(
    spec = list(method = "REML", family = list(family = "gaussian")),
    payload = list(y = seq_len(10)),
    result = result
  )
  ce <- common_env_proportion(fit)
  expect_equal(ce$estimate, 0.23)
  expect_match(attr(ce, "interpretation"), "variance ratio", fixed = TRUE)
  expect_match(attr(ce, "interpretation"), "NOT a heritability", fixed = TRUE)
})

test_that("common_env_proportion_interval() returns the ratio2 CI when present", {
  ci <- data.frame(
    estimate = 0.23,
    lower = 0.10,
    upper = 0.44,
    level = 0.95,
    se = 0.08,
    lower_clamped = FALSE,
    upper_clamped = FALSE,
    boundary = FALSE,
    stringsAsFactors = FALSE
  )
  fit <- hsquared:::hs_new_fit(
    spec = list(method = "REML", family = list(family = "gaussian")),
    payload = list(y = seq_len(10)),
    result = list(common_env_proportion_interval = ci)
  )
  out <- common_env_proportion_interval(fit)
  expect_equal(out$estimate, 0.23)
  expect_true(out$lower <= out$estimate && out$estimate <= out$upper)
  expect_false(out$boundary)
  expect_match(attr(out, "interpretation"), "variance ratio", fixed = TRUE)
})

test_that("common_env_proportion_interval() errors clearly without the field", {
  fit_no_ci <- hsquared:::hs_new_fit(
    spec = list(method = "REML", family = list(family = "gaussian")),
    payload = list(y = seq_len(10)),
    result = list(common_env_proportion = data.frame(term = "common_env", estimate = 0.2))
  )
  expect_error(
    common_env_proportion_interval(fit_no_ci),
    "common-environment variance-ratio interval"
  )
})

test_that("hs_normalize_two_effect_ratio_interval() builds the two shapes", {
  # ratio1 (h2) carries a method column; ratio2 (c2) omits it.
  r1 <- hsquared:::hs_normalize_two_effect_ratio_interval(
    estimate = 0.42, lower = 0.21, upper = 0.66, se = 0.11,
    lower_clamped = FALSE, upper_clamped = FALSE, boundary = FALSE,
    level = 0.95, with_method = TRUE
  )
  expect_equal(nrow(r1), 1L)
  expect_equal(r1$estimate, 0.42)
  expect_equal(r1$method, "delta")
  expect_true(r1$lower <= r1$estimate && r1$estimate <= r1$upper)

  r2 <- hsquared:::hs_normalize_two_effect_ratio_interval(
    estimate = 0.23, lower = 0.10, upper = 0.44, se = 0.08,
    lower_clamped = FALSE, upper_clamped = FALSE, boundary = FALSE,
    level = 0.95, with_method = FALSE
  )
  expect_false("method" %in% names(r2))
  expect_true(all(c("lower_clamped", "upper_clamped", "boundary") %in% names(r2)))
})

test_that("hs_normalize_two_effect_ratio_interval() flags a boundary component", {
  # Engine returns NaN bounds at a sigma -> 0 boundary; normalize to NA, flagged.
  bd <- hsquared:::hs_normalize_two_effect_ratio_interval(
    estimate = 0.0, lower = NaN, upper = NaN, se = NaN,
    lower_clamped = FALSE, upper_clamped = FALSE, boundary = TRUE,
    level = 0.95, with_method = FALSE
  )
  expect_true(is.na(bd$lower))
  expect_true(is.na(bd$upper))
  expect_true(is.na(bd$se))
  expect_true(bd$boundary)
})

test_that("hs_attach_two_effect_intervals routes ratio2 to the common-env field", {
  raw_ci <- list(
    level = 0.95,
    r1_estimate = 0.42, r1_lower = 0.21, r1_upper = 0.66, r1_se = 0.11,
    r1_lower_clamped = FALSE, r1_upper_clamped = FALSE, r1_boundary = FALSE,
    r2_estimate = 0.23, r2_lower = 0.10, r2_upper = 0.44, r2_se = 0.08,
    r2_lower_clamped = FALSE, r2_upper_clamped = FALSE, r2_boundary = FALSE
  )
  payload <- list(effect2 = list(type = "common_env"))
  result <- hsquared:::hs_attach_two_effect_intervals(list(), raw_ci, payload)
  # ratio1 -> heritability_interval (so heritability_interval() resolves);
  # ratio2 -> common_env_proportion_interval.
  expect_true(!is.null(result$heritability_interval))
  expect_equal(result$heritability_interval$estimate, 0.42)
  expect_equal(result$heritability_interval$method, "delta")
  expect_true(!is.null(result$common_env_proportion_interval))
  expect_equal(result$common_env_proportion_interval$estimate, 0.23)
  expect_null(result$maternal_proportion_interval)
})

test_that("heritability_interval() resolves on a two-effect fit", {
  # After the bridge attaches ratio1 to `heritability_interval`, the existing
  # extractor works on a two-effect fit with no special-casing.
  fit <- hsquared:::hs_new_fit(
    spec = list(method = "REML", family = list(family = "gaussian"),
                target = "two_effect"),
    payload = list(y = seq_len(10)),
    result = list(
      heritability_interval = data.frame(
        estimate = 0.42, lower = 0.21, upper = 0.66, level = 0.95,
        se = 0.11, lower_clamped = FALSE, upper_clamped = FALSE,
        boundary = FALSE, method = "delta", stringsAsFactors = FALSE
      )
    )
  )
  hi <- heritability_interval(fit)
  expect_equal(hi$estimate, 0.42)
  expect_equal(hi$method, "delta")
})

test_that("hsquared fits the opt-in common-environment model", {
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for live two-effect."
  )

  # Founders a-d plus offspring e-h, grouped into litters; the common-environment
  # effect (litter) is shared by litter-mates and independent across litters.
  ped <- data.frame(
    id = c("a", "b", "c", "d", "e", "f", "g", "h"),
    sire = c(NA, NA, NA, NA, "a", "a", "c", "c"),
    dam = c(NA, NA, NA, NA, "b", "b", "d", "d")
  )
  set.seed(7)
  ids <- ped$id
  litter <- c("l1", "l1", "l2", "l2", "l3", "l3", "l4", "l4")
  ce <- stats::setNames(stats::rnorm(4, 0, 0.8), c("l1", "l2", "l3", "l4"))
  dat <- data.frame(
    y = 3 + ce[litter] + stats::rnorm(8, 0, 0.7),
    id = ids,
    litter = litter
  )

  fit <- hsquared(
    y ~ animal(1 | id, pedigree = ped) + common_env(1 | litter),
    data = dat,
    family = stats::gaussian(),
    control = hs_control(
      engine = "julia",
      engine_control = list(target = "two_effect")
    )
  )

  expect_s3_class(fit, "hsquared_fit")
  expect_equal(fit$spec$target, "two_effect")
  vc <- variance_components(fit)
  expect_equal(vc$component, c("animal", "common_env", "residual"))
  expect_true(all(is.finite(vc$estimate)) && all(vc$estimate > 0))
  h2 <- heritability(fit)$estimate
  expect_true(is.finite(h2) && h2 >= 0 && h2 < 1)
  expect_equal(nrow(common_env_effects(fit)), 4L)
  expect_equal(
    fit_diagnostics(fit)$value[
      fit_diagnostics(fit)$metric == "variance_components_source"
    ],
    "estimated_two_effect_reml"
  )

  # R1a: reachable c2 accessor with the Falconer fence attached.
  ce <- common_env_proportion(fit)
  expect_equal(nrow(ce), 1L)
  c2 <- ce$estimate
  expect_true(is.finite(c2) && c2 >= 0 && c2 < 1)
  expect_match(attr(ce, "interpretation"), "variance ratio", fixed = TRUE)

  # R1b: the ratio interval is wired opportunistically. When the engine returns
  # it (a non-boundary, positive-definite optimum on this small design), both
  # heritability_interval() (ratio1) and common_env_proportion_interval()
  # (ratio2) resolve; otherwise the fields are absent (documented boundary
  # behavior), so guard the assertions on presence rather than asserting the
  # engine always produces a CI on 8 records.
  if (!is.null(fit$result$heritability_interval)) {
    hi <- heritability_interval(fit)
    expect_equal(hi$estimate, h2, tolerance = 1e-6)
    expect_equal(hi$method, "delta")
    expect_true(hi$boundary || (hi$lower <= hi$estimate && hi$estimate <= hi$upper))
  }
  if (!is.null(fit$result$common_env_proportion_interval)) {
    ci <- common_env_proportion_interval(fit)
    expect_equal(ci$estimate, c2, tolerance = 1e-6)
    expect_true(ci$boundary || (ci$lower <= ci$estimate && ci$estimate <= ci$upper))
    expect_match(attr(ci, "interpretation"), "variance ratio", fixed = TRUE)
  }
})
