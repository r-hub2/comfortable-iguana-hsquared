# CI-runnable honesty tests (no Julia required). These pin the user-facing
# wording and the reported Julia targets for three fixes:
#   #6 engine-setup error names every recovery lever an applied user needs;
#   #2 a genuinely-unsupported marker points at `formula_status()` and does not
#      over-claim that only `animal(...)` parses;
#   #1 a multivariate spec reports the multivariate Julia fit target, matching
#      the `spec$bridge$target` source of truth rather than the univariate one.
# None of these touch the live bridge, so they must run on CI without Julia.

test_that("the engine-setup error names every install and recovery lever", {
  ped <- data.frame(
    id = c("a", "b", "c"),
    sire = c(NA, NA, "a"),
    dam = c(NA, NA, "b")
  )
  dat <- data.frame(
    y = c(1, 2, 3),
    sex = c("m", "f", "m"),
    id = c("a", "b", "c")
  )

  # `engine = "fit"` with a bogus `julia_project` cannot find the engine, so it
  # must error with actionable setup guidance rather than fitting nothing. The
  # message has to spell out, by exact literal, each lever an applied user can
  # pull: what is missing, the env-var and `engine_control` escape hatches, how
  # to obtain the source, and the no-Julia `validate` fallback.
  err <- expect_error(
    hsquared(
      y ~ sex + animal(1 | id, pedigree = ped),
      data = dat,
      control = hs_control(
        engine = "fit",
        engine_control = list(julia_project = tempfile())
      )
    )
  )

  expect_match(
    conditionMessage(err),
    "requires the HSquared.jl Julia",
    fixed = TRUE
  )
  expect_match(
    conditionMessage(err),
    "HSQUARED_JULIA_PROJECT",
    fixed = TRUE
  )
  expect_match(conditionMessage(err), "julia_project", fixed = TRUE)
  expect_match(conditionMessage(err), "engine_control", fixed = TRUE)
  expect_match(conditionMessage(err), "git clone", fixed = TRUE)
  expect_match(
    conditionMessage(err),
    "https://github.com/itchyshin/HSquared.jl",
    fixed = TRUE
  )
  expect_match(conditionMessage(err), "validate", fixed = TRUE)
})

test_that("an unsupported marker points at formula_status without over-claiming", {
  ped <- data.frame(
    id = c("a", "b", "c", "d"),
    sire = c(NA, NA, "a", "a"),
    dam = c(NA, NA, "b", "c")
  )
  dat <- data.frame(
    y = c(1, 2, 3),
    sex = c("f", "m", "f"),
    id = c("a", "c", "d")
  )

  # `epistasis()` is a genuinely planned-but-unimplemented marker. Even on the
  # no-fit `validate` path it must error during spec construction, before the
  # validate stop, and route the user to the live status table.
  err <- expect_error(
    hsquared(
      y ~ sex + epistasis(1 | id),
      data = dat,
      control = hs_control(engine = "validate")
    )
  )

  expect_true(grepl("formula_status", conditionMessage(err), fixed = TRUE))

  # The message must not over-claim that only the additive-genetic term parses:
  # several other terms parse (e.g. `permanent`, `genomic`), and `formula_status()`
  # is the honest source. Pin the absence of the over-narrow literal.
  expect_false(grepl(
    "only `animal(1 | id, pedigree = ped)`",
    conditionMessage(err),
    fixed = TRUE
  ))
})

test_that("a multivariate spec reports the multivariate Julia fit target", {
  ped <- data.frame(
    id = c("sire", "dam", "calf1", "calf2"),
    sire = c(NA, NA, "sire", "sire"),
    dam = c(NA, NA, "dam", "dam")
  )
  dat <- data.frame(
    y1 = c(1, 2, 3, 4),
    y2 = c(1.5, NA, 3.5, 4.5),
    sex = c("m", "f", "f", "m"),
    id = ped$id
  )
  mv_formula <- cbind(y1, y2) ~ sex + animal(1 | id, pedigree = ped)

  # Source of truth: the model spec records the multivariate REML target.
  spec <- hsquared:::hs_build_model_spec(
    mv_formula,
    data = dat,
    family = stats::gaussian(),
    REML = TRUE
  )
  expect_match(spec$bridge$target, "fit_multivariate_reml", fixed = TRUE)

  # The no-fit `validate` path reports the Julia fit target in its message. For a
  # `cbind(...)` response it must name the multivariate target, not the
  # univariate `fit_animal_model` one.
  expect_message(
    validated_spec <- hsquared(
      mv_formula,
      data = dat,
      control = hs_control(engine = "validate")
    ),
    "fit_multivariate_reml",
    fixed = TRUE
  )
  expect_match(
    validated_spec$bridge$target,
    "fit_multivariate_reml",
    fixed = TRUE
  )

  # The public `model_spec()` inspector surfaces the same fit target, and its
  # `summary()` "fit" row must agree with it.
  spec_obj <- model_spec(mv_formula, data = dat)
  expect_match(spec_obj$julia$fit_target, "fit_multivariate_reml", fixed = TRUE)

  fit_row <- summary(spec_obj)$julia
  reported_fit_target <- fit_row$target[fit_row$stage == "fit"]
  expect_match(reported_fit_target, "fit_multivariate_reml", fixed = TRUE)

  # The reported target must equal the `spec$bridge$target` source of truth (the
  # inspector only namespaces it with the engine module). Equality, not a mere
  # substring, would catch a silent regression to the univariate target.
  expect_identical(
    spec_obj$julia$fit_target,
    paste0("HSquared.", spec$bridge$target)
  )
  expect_identical(reported_fit_target, spec_obj$julia$fit_target)
})
