# Opt-in, experimental maternal-genetic two-effect model: direct additive
# genetic effect + maternal additive genetic effect (both pedigree A, the
# maternal effect expressed through the dam). Surfaces fit_two_effect_reml with
# a pedigree relationship for the second effect. REML only.

test_that("the parser accepts maternal_genetic() as a pedigree second effect", {
  ped <- data.frame(
    id = c("a", "b", "c", "d", "e"),
    sire = c(NA, NA, NA, "a", "a"),
    dam = c(NA, NA, NA, "b", "c")
  )
  dat <- data.frame(
    y = c(1, 2, 3, 4),
    id = c("b", "c", "d", "e"),
    mum = c("a", "a", "b", "c")
  )

  spec <- hsquared:::hs_build_model_spec(
    y ~ animal(1 | id, pedigree = ped) + maternal_genetic(1 | mum),
    data = dat,
    family = stats::gaussian(),
    REML = TRUE
  )

  expect_named(spec$random, c("animal", "maternal_genetic"))
  expect_equal(spec$random$maternal_genetic$type, "maternal_genetic")
  expect_equal(spec$random$maternal_genetic$group, "mum")
  expect_equal(spec$random$maternal_genetic$relationship, "pedigree")
  expect_match(spec$bridge$target, "two_effect", fixed = TRUE)
})

test_that("maternal_genetic() dams must be in the pedigree", {
  ped <- data.frame(id = c("a", "b"), sire = c(NA, NA), dam = c(NA, NA))
  dat <- data.frame(y = c(1, 2), id = c("a", "b"), mum = c("a", "ghost"))
  expect_error(
    hsquared:::hs_build_model_spec(
      y ~ animal(1 | id, pedigree = ped) + maternal_genetic(1 | mum),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "not in the pedigree",
    fixed = TRUE
  )
})

hs_expect_maternal_parser_next_call <- function(err) {
  expect_true(grepl("column of `data`", err, fixed = TRUE))
  expect_true(grepl("mothers of the records", err, fixed = TRUE))
  expect_true(grepl("animal pedigree", err, fixed = TRUE))
  expect_true(grepl(
    "does not take dams from the pedigree table",
    err,
    fixed = TRUE
  ))
  expect_true(grepl("Closest working call", err, fixed = TRUE))
  expect_true(grepl("target = \"direct_maternal\"", err, fixed = TRUE))
  expect_true(grepl("target = \"two_effect\"", err, fixed = TRUE))
  expect_true(grepl("Willham triple", err, fixed = TRUE))
  expect_true(grepl("not a scalar h2", err, fixed = TRUE))
  expect_true(grepl("neither is auto-routed", err, fixed = TRUE))
  expect_lt(
    regexpr("direct_maternal", err, fixed = TRUE)[[1L]],
    regexpr("two_effect", err, fixed = TRUE)[[1L]]
  )
}

test_that("extra pedigree= on maternal_genetic() pastes the covered next call", {
  ped <- data.frame(id = c("a", "b"), sire = c(NA, NA), dam = c(NA, NA))
  dat <- data.frame(y = c(1, 2), id = c("a", "b"))
  err <- tryCatch(
    hsquared:::hs_build_model_spec(
      y ~ animal(1 | id, pedigree = ped) +
        maternal_genetic(1 | dam, pedigree = ped),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    error = function(e) conditionMessage(e)
  )
  expect_true(grepl("takes no extra arguments", err, fixed = TRUE))
  expect_true(grepl("pedigree=", err, fixed = TRUE))
  expect_false(grepl("was not found in `data`", err, fixed = TRUE))
  hs_expect_maternal_parser_next_call(err)
})

test_that("missing dam column on maternal_genetic() pastes the covered next call", {
  ped <- data.frame(id = c("a", "b"), sire = c(NA, NA), dam = c(NA, NA))
  dat <- data.frame(y = c(1, 2), id = c("a", "b"))
  err <- tryCatch(
    hsquared:::hs_build_model_spec(
      y ~ animal(1 | id, pedigree = ped) + maternal_genetic(1 | dam),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    error = function(e) conditionMessage(e)
  )
  expect_true(grepl(
    "`maternal_genetic()` grouping variable `dam` was not found in `data`",
    err,
    fixed = TRUE
  ))
  expect_false(grepl("takes no extra arguments", err, fixed = TRUE))
  hs_expect_maternal_parser_next_call(err)
})

test_that("maternal_genetic() must be intercept-only", {
  ped <- data.frame(id = c("a", "b"), sire = c(NA, NA), dam = c(NA, NA))
  dat <- data.frame(y = c(1, 2), id = c("a", "b"), mum = c("a", "b"))
  expect_error(
    hsquared:::hs_build_model_spec(
      y ~ animal(1 | id, pedigree = ped) + maternal_genetic(age | mum),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "random-intercept",
    fixed = TRUE
  )
})

test_that("the default engine = \"fit\" rejects a maternal_genetic() formula", {
  ped <- data.frame(id = c("a", "b"), sire = c(NA, NA), dam = c(NA, NA))
  dat <- data.frame(y = c(1, 2), id = c("a", "b"), mum = c("a", "b"))
  err <- tryCatch(
    hsquared(
      y ~ animal(1 | id, pedigree = ped) + maternal_genetic(1 | mum),
      data = dat,
      family = stats::gaussian()
    ),
    error = function(e) conditionMessage(e)
  )
  expect_true(grepl(
    "`maternal_genetic()` is not on the default path",
    err,
    fixed = TRUE
  ))
  expect_true(grepl("Closest working call", err, fixed = TRUE))
  expect_true(grepl("target = \"two_effect\"", err, fixed = TRUE))
  expect_true(grepl("target = \"direct_maternal\"", err, fixed = TRUE))
  expect_true(grepl("Willham triple", err, fixed = TRUE))
  expect_true(grepl("not a scalar h2", err, fixed = TRUE))
  expect_true(grepl("neither is auto-routed", err, fixed = TRUE))
  # Covered sibling is named first; two_effect stays the experimental alternative.
  expect_lt(
    regexpr("direct_maternal", err, fixed = TRUE)[[1L]],
    regexpr("two_effect", err, fixed = TRUE)[[1L]]
  )
})

test_that("hs_second_effect_target does not auto-route maternal_genetic", {
  expect_identical(
    hsquared:::hs_second_effect_target("maternal_genetic"),
    "two_effect"
  )
})

test_that("a maternal_genetic() formula needs the two_effect target", {
  ped <- data.frame(id = c("a", "b"), sire = c(NA, NA), dam = c(NA, NA))
  dat <- data.frame(y = c(1, 2), id = c("a", "b"), mum = c("a", "b"))
  err <- tryCatch(
    hsquared(
      y ~ animal(1 | id, pedigree = ped) + maternal_genetic(1 | mum),
      data = dat,
      family = stats::gaussian(),
      control = hs_control(
        engine = "julia",
        engine_control = list(target = "ai_reml")
      )
    ),
    error = function(e) conditionMessage(e)
  )
  expect_true(grepl("needs `target = \"two_effect\"`", err, fixed = TRUE))
  expect_true(grepl("direct_maternal", err, fixed = TRUE))
  expect_true(grepl("is covered", err, fixed = TRUE))
  expect_true(grepl("is experimental", err, fixed = TRUE))
})

test_that("maternal_proportion / interval extractors require a fitted object", {
  expect_error(
    maternal_proportion(1),
    "requires an `hsquared_fit`",
    fixed = TRUE
  )
  expect_error(
    maternal_proportion_interval(1),
    "requires an `hsquared_fit`",
    fixed = TRUE
  )
})

test_that("maternal_proportion() returns m2 with a Falconer interpretation", {
  fit <- hsquared:::hs_new_fit(
    spec = list(method = "REML", family = list(family = "gaussian")),
    payload = list(y = seq_len(10)),
    result = list(
      maternal_proportion = data.frame(
        term = "maternal_genetic",
        estimate = 0.18
      )
    )
  )
  m <- maternal_proportion(fit)
  expect_equal(m$estimate, 0.18)
  expect_match(attr(m, "interpretation"), "maternal", fixed = TRUE)
  expect_match(attr(m, "interpretation"), "NOT a heritability", fixed = TRUE)
})

test_that("hs_attach_two_effect_intervals routes ratio2 to the maternal field", {
  raw_ci <- list(
    level = 0.95,
    r1_estimate = 0.42,
    r1_lower = 0.21,
    r1_upper = 0.66,
    r1_se = 0.11,
    r1_lower_clamped = FALSE,
    r1_upper_clamped = FALSE,
    r1_boundary = FALSE,
    r2_estimate = 0.18,
    r2_lower = 0.05,
    r2_upper = 0.48,
    r2_se = 0.09,
    r2_lower_clamped = FALSE,
    r2_upper_clamped = FALSE,
    r2_boundary = FALSE
  )
  payload <- list(effect2 = list(type = "maternal_genetic"))
  result <- hsquared:::hs_attach_two_effect_intervals(list(), raw_ci, payload)
  expect_true(!is.null(result$maternal_proportion_interval))
  expect_equal(result$maternal_proportion_interval$estimate, 0.18)
  expect_null(result$common_env_proportion_interval)
  # ratio1 still feeds heritability_interval regardless of the second-effect type
  expect_equal(result$heritability_interval$estimate, 0.42)
})

test_that("maternal_proportion_interval() returns the ratio2 CI with the fence", {
  ci <- data.frame(
    estimate = 0.18,
    lower = 0.05,
    upper = 0.48,
    level = 0.95,
    se = 0.09,
    lower_clamped = FALSE,
    upper_clamped = FALSE,
    boundary = FALSE,
    stringsAsFactors = FALSE
  )
  fit <- hsquared:::hs_new_fit(
    spec = list(method = "REML", family = list(family = "gaussian")),
    payload = list(y = seq_len(10)),
    result = list(maternal_proportion_interval = ci)
  )
  out <- maternal_proportion_interval(fit)
  expect_equal(out$estimate, 0.18)
  expect_match(attr(out, "interpretation"), "maternal", fixed = TRUE)
})

test_that("hsquared fits the opt-in maternal-genetic model", {
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for live maternal fit."
  )

  # Founders a-f; offspring g-n each with a dam among the founders/earlier
  # animals. The maternal effect is expressed through the dam (pedigree A2).
  ped <- data.frame(
    id = c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l"),
    sire = c(NA, NA, NA, NA, NA, NA, "a", "a", "c", "c", "e", "e"),
    dam = c(NA, NA, NA, NA, NA, NA, "b", "b", "d", "d", "f", "f")
  )
  set.seed(3)
  rec_id <- c("g", "h", "i", "j", "k", "l")
  rec_mum <- c("b", "b", "d", "d", "f", "f")
  mat <- stats::setNames(
    stats::rnorm(6, 0, 0.8),
    c("a", "b", "c", "d", "e", "f")
  )
  dat <- data.frame(
    y = 4 + mat[rec_mum] + stats::rnorm(6, 0, 0.7),
    id = rec_id,
    mum = rec_mum
  )

  fit <- hsquared(
    y ~ animal(1 | id, pedigree = ped) + maternal_genetic(1 | mum),
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
  expect_equal(vc$component, c("animal", "maternal_genetic", "residual"))
  expect_true(all(is.finite(vc$estimate)) && all(vc$estimate > 0))
  expect_equal(
    fit_diagnostics(fit)$value[
      fit_diagnostics(fit)$metric == "variance_components_source"
    ],
    "estimated_two_effect_reml"
  )

  # R1a: reachable m2 accessor with the Falconer fence.
  m <- maternal_proportion(fit)
  expect_equal(nrow(m), 1L)
  m2 <- m$estimate
  expect_true(is.finite(m2) && m2 >= 0 && m2 < 1)
  expect_match(attr(m, "interpretation"), "maternal", fixed = TRUE)

  # R1b: opportunistic ratio interval (guard on presence; a boundary/flat
  # optimum legitimately omits the CI on this small design).
  if (!is.null(fit$result$maternal_proportion_interval)) {
    mi <- maternal_proportion_interval(fit)
    expect_equal(mi$estimate, m2, tolerance = 1e-6)
    expect_true(
      mi$boundary || (mi$lower <= mi$estimate && mi$estimate <= mi$upper)
    )
  }
})
