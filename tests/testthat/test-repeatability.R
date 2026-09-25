# Opt-in, experimental repeatability (permanent-environment) model.
# Mirrors the sparse_reml/ai_reml opt-in pattern: the default engine = "fit"
# stays single-effect; the repeatability model is reachable only via the
# explicit `engine = "julia", target = "repeatability"` path. REML only.

test_that("the parser accepts permanent() alongside animal() as the PE effect", {
  ped <- data.frame(
    id = c("a", "b", "c"),
    sire = c(NA, NA, "a"),
    dam = c(NA, NA, "b")
  )
  # repeated records: a and b each measured twice (needed for identifiability)
  dat <- data.frame(
    y = c(1, 2, 3, 1.5, 2.5),
    id = c("a", "b", "c", "a", "b")
  )

  spec <- hsquared:::hs_build_model_spec(
    y ~ animal(1 | id, pedigree = ped) + permanent(1 | id),
    data = dat,
    family = stats::gaussian(),
    REML = TRUE
  )

  expect_named(spec$random, c("animal", "permanent"))
  expect_equal(spec$random$permanent$type, "permanent")
  expect_equal(spec$random$permanent$group, "id")
  expect_equal(spec$random$permanent$design, "intercept")
  expect_equal(spec$random$permanent$relationship, "identity")
  # the PE effect shares the animal incidence (same per-record id values)
  expect_equal(spec$random$permanent$values, spec$random$animal$values)
  expect_match(spec$bridge$target, "repeatability", fixed = TRUE)
})

test_that("a single animal() formula still builds a single-effect spec", {
  ped <- data.frame(id = c("a", "b"), sire = c(NA, NA), dam = c(NA, NA))
  dat <- data.frame(y = c(1, 2), id = c("a", "b"))

  spec <- hsquared:::hs_build_model_spec(
    y ~ animal(1 | id, pedigree = ped),
    data = dat,
    family = stats::gaussian(),
    REML = TRUE
  )
  expect_named(spec$random, "animal")
  expect_false("repeatability" %in% spec$bridge$target)
})

test_that("permanent() must share the animal grouping and be intercept-only", {
  ped <- data.frame(id = c("a", "b"), sire = c(NA, NA), dam = c(NA, NA))
  dat <- data.frame(y = c(1, 2), id = c("a", "b"), cage = c("c1", "c2"))

  expect_error(
    hsquared:::hs_build_model_spec(
      y ~ animal(1 | id, pedigree = ped) + permanent(1 | cage),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "same grouping variable",
    fixed = TRUE
  )

  expect_error(
    hsquared:::hs_build_model_spec(
      y ~ animal(1 | id, pedigree = ped) + permanent(cage | id),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "random-intercept",
    fixed = TRUE
  )
})

test_that("only one second random effect is supported", {
  ped <- data.frame(id = c("a", "b"), sire = c(NA, NA), dam = c(NA, NA))
  dat <- data.frame(y = c(1, 2), id = c("a", "b"))

  expect_error(
    hsquared:::hs_build_model_spec(
      y ~ animal(1 | id, pedigree = ped) +
        permanent(1 | id) +
        permanent(1 | id),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "at most one additional random effect",
    fixed = TRUE
  )
})

test_that("other planned QG markers are still rejected", {
  ped <- data.frame(id = c("a", "b"), sire = c(NA, NA), dam = c(NA, NA))
  dat <- data.frame(y = c(1, 2), id = c("a", "b"))

  expect_error(
    hsquared:::hs_build_model_spec(
      y ~ animal(1 | id, pedigree = ped) + dominance(1 | id),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "`dominance()` is planned, not implemented.",
    fixed = TRUE
  )
})

test_that("repeatability is a valid opt-in julia target", {
  expect_equal(
    hsquared:::hs_validate_julia_target("repeatability"),
    "repeatability"
  )
})

test_that("the default engine = \"fit\" rejects a permanent() formula", {
  ped <- data.frame(id = c("a", "b"), sire = c(NA, NA), dam = c(NA, NA))
  dat <- data.frame(y = c(1, 2, 1.5, 2.5), id = c("a", "b", "a", "b"))
  expect_error(
    hsquared(
      y ~ animal(1 | id, pedigree = ped) + permanent(1 | id),
      data = dat,
      family = stats::gaussian()
    ),
    "Closest working call",
    fixed = TRUE
  )
  err <- tryCatch(
    hsquared(
      y ~ animal(1 | id, pedigree = ped) + permanent(1 | id),
      data = dat,
      family = stats::gaussian()
    ),
    error = function(e) conditionMessage(e)
  )
  expect_true(grepl(
    "`permanent()` is not on the default path",
    err,
    fixed = TRUE
  ))
  expect_true(grepl("target = \"repeatability\"", err, fixed = TRUE))
  expect_true(grepl("That route is experimental", err, fixed = TRUE))
})

test_that("target = \"repeatability\" requires a permanent() term", {
  ped <- data.frame(id = c("a", "b"), sire = c(NA, NA), dam = c(NA, NA))
  dat <- data.frame(y = c(1, 2), id = c("a", "b"))
  expect_error(
    hsquared(
      y ~ animal(1 | id, pedigree = ped),
      data = dat,
      family = stats::gaussian(),
      control = hs_control(
        engine = "julia",
        engine_control = list(target = "repeatability")
      )
    ),
    "requires a `permanent(1 | id)` term",
    fixed = TRUE
  )
})

test_that("a permanent() formula needs the repeatability target", {
  ped <- data.frame(id = c("a", "b"), sire = c(NA, NA), dam = c(NA, NA))
  dat <- data.frame(y = c(1, 2, 1.5, 2.5), id = c("a", "b", "a", "b"))
  expect_error(
    hsquared(
      y ~ animal(1 | id, pedigree = ped) + permanent(1 | id),
      data = dat,
      family = stats::gaussian(),
      control = hs_control(
        engine = "julia",
        engine_control = list(target = "ai_reml")
      )
    ),
    "needs `target = \"repeatability\"`",
    fixed = TRUE
  )
})

test_that("the repeatability target rejects REML = FALSE (ML)", {
  ped <- data.frame(id = c("a", "b"), sire = c(NA, NA), dam = c(NA, NA))
  dat <- data.frame(y = c(1, 2, 1.5, 2.5), id = c("a", "b", "a", "b"))
  # REML-only: this is a pure request-validity error and fires before any
  # Julia-engine call, so it needs no live bridge.
  expect_error(
    hsquared(
      y ~ animal(1 | id, pedigree = ped) + permanent(1 | id),
      data = dat,
      family = stats::gaussian(),
      REML = FALSE,
      control = hs_control(
        engine = "julia",
        engine_control = list(target = "repeatability")
      )
    ),
    "ML estimation",
    fixed = TRUE
  )
})

test_that("repeatability extractors require a fitted object", {
  expect_error(repeatability(1), "requires an `hsquared_fit`", fixed = TRUE)
  expect_error(
    permanent_effects(1),
    "requires an `hsquared_fit`",
    fixed = TRUE
  )
})

test_that("the repeatability bridge requires an internal payload", {
  expect_error(
    hsquared:::hs_fit_julia_repeatability_payload(list()),
    "`payload` must be an internal `hs_bridge_payload`.",
    fixed = TRUE
  )
})

test_that("the repeatability initial validator enforces three named components", {
  expect_error(
    hsquared:::hs_validate_repeatability_initial(c(sigma_a2 = 1, sigma_e2 = 1)),
    "sigma_a2`, `sigma_pe2`, and `sigma_e2`",
    fixed = TRUE
  )
  expect_error(
    hsquared:::hs_validate_repeatability_initial(
      c(sigma_a2 = 1, sigma_pe2 = -1, sigma_e2 = 1)
    ),
    "finite and positive",
    fixed = TRUE
  )
})

test_that("hsquared fits the opt-in repeatability model on repeated records", {
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for live repeatability."
  )

  # A related pedigree (d, e are offspring of the founders, so A != I) with
  # repeated records (each animal measured three times). Relatedness plus
  # repeated records is what lets the additive (Va) and permanent-environment
  # (Vpe) variances be separated.
  ped <- data.frame(
    id = c("a", "b", "c", "d", "e"),
    sire = c(NA, NA, NA, "a", "a"),
    dam = c(NA, NA, NA, "b", "c")
  )
  set.seed(11)
  ids <- rep(c("a", "b", "c", "d", "e"), each = 3)
  pe <- stats::setNames(stats::rnorm(5, 0, 1.0), c("a", "b", "c", "d", "e"))
  dat <- data.frame(
    y = 2 + pe[ids] + stats::rnorm(15, 0, 0.7),
    id = ids
  )

  # PINNED (hsquared#235): on this tiny fixture the engine's REML information
  # is not positive definite, so repeatability_interval() refuses and the
  # bridge surfaces it as one warning instead of a silent NA. Asserting it
  # means the surfacing cannot regress to silence unnoticed. If the engine
  # later makes this interval computable, this expectation FAILS -- that is
  # the intent: update it deliberately, do not muffle it.
  expect_warning(
    fit <- hsquared(
      y ~ animal(1 | id, pedigree = ped) + permanent(1 | id),
      data = dat,
      family = stats::gaussian(),
      control = hs_control(
        engine = "julia",
        engine_control = list(target = "repeatability")
      )
    ),
    "repeatability_interval\\(\\).*REML information is not positive definite"
  )

  expect_s3_class(fit, "hsquared_fit")
  expect_equal(fit$spec$target, "repeatability")
  vc <- variance_components(fit)
  expect_equal(vc$component, c("animal", "permanent", "residual"))
  expect_true(all(is.finite(vc$estimate)) && all(vc$estimate > 0))
  r <- repeatability(fit)$estimate
  expect_true(is.finite(r) && r > 0 && r < 1)
  h2 <- heritability(fit)$estimate
  expect_true(is.finite(h2) && h2 >= 0 && h2 < 1)
  expect_equal(nrow(permanent_effects(fit)), 5L)
  expect_equal(
    fit_diagnostics(fit)$value[
      fit_diagnostics(fit)$metric == "variance_components_source"
    ],
    "estimated_repeatability_reml"
  )
})
