# Genuine unsupported-syntax / grammar / unsupported-target rejections carry the
# structured `hsquared_unsupported_syntax` condition class (raised via
# `hs_abort_unsupported_syntax()`); data-validation / wrong-object guards stay
# plain `stop()` and must NOT carry it. These tests pin a representative sample
# across the reclassed sites in R/model-spec.R, R/hsquared.R, R/julia-bridge.R,
# and R/model-spec-inspect.R.

hs_usc_ped <- function() {
  data.frame(
    id = c("a", "b", "c", "d"),
    sire = c(NA, NA, "a", "a"),
    dam = c(NA, NA, "b", "b")
  )
}

hs_usc_dat <- function() {
  data.frame(
    y = c(1, 2, 3, 4),
    id = c("a", "b", "c", "d"),
    sex = c("f", "m", "f", "m"),
    age = c(1, 2, 3, 4),
    litter = c("l1", "l1", "l2", "l2"),
    dam = c("a", "a", "b", "b")
  )
}

test_that("the structured condition carries the documented class stack", {
  ped <- hs_usc_ped()
  dat <- hs_usc_dat()

  cnd <- expect_error(
    hsquared:::hs_build_model_spec(
      y ~ animal(1 | id, pedigree = ped, cov = us()),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    class = "hsquared_unsupported_syntax"
  )
  expect_s3_class(cnd, "hsquared_error")
  expect_s3_class(cnd, "error")
  expect_true(inherits(
    cnd,
    c(
      "hsquared_unsupported_syntax",
      "hsquared_error",
      "error",
      "condition"
    ),
    which = FALSE
  ))
  # `call. = FALSE` ergonomics preserved: no call attached to the condition.
  expect_null(conditionCall(cnd))
})

test_that("model-spec grammar/planned-term rejections carry the class", {
  ped <- hs_usc_ped()
  dat <- hs_usc_dat()
  Ginv <- diag(4)
  dimnames(Ginv) <- list(dat$id, dat$id)

  # planned `animal(cov = ...)` term
  expect_error(
    hsquared:::hs_build_model_spec(
      y ~ animal(1 | id, pedigree = ped, cov = us()),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    class = "hsquared_unsupported_syntax"
  )

  # planned offset() term
  expect_error(
    hsquared:::hs_build_model_spec(
      y ~ offset(age) + animal(1 | id, pedigree = ped),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    class = "hsquared_unsupported_syntax"
  )

  # unsupported family/target on this path
  expect_error(
    hsquared:::hs_build_model_spec(
      y ~ sex + animal(1 | id, pedigree = ped),
      data = dat,
      family = stats::gaussian("log"),
      REML = TRUE
    ),
    class = "hsquared_unsupported_syntax"
  )

  # more than one primary effect (formula-grammar cardinality)
  expect_error(
    hsquared:::hs_build_model_spec(
      y ~ animal(1 | id, pedigree = ped) + genomic(1 | id, Ginv = Ginv),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    class = "hsquared_unsupported_syntax"
  )

  # planned marker() term
  expect_error(
    hsquared:::hs_build_model_spec(
      y ~ sex + markers(M, model = "random"),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    class = "hsquared_unsupported_syntax"
  )

  # recognized effect/marker call nested inside an interaction or function call
  # (formula-grammar placement rejection)
  expect_error(
    hsquared:::hs_build_model_spec(
      y ~ sex:dominance(1 | id) + animal(1 | id, pedigree = ped),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    class = "hsquared_unsupported_syntax"
  )

  # random-slope i.i.d. effect (planned, not implemented)
  expect_error(
    hsquared:::hs_build_model_spec(
      y ~ animal(1 | id, pedigree = ped) + (age | litter),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    class = "hsquared_unsupported_syntax"
  )

  # common-environment slope (planned, not implemented)
  expect_error(
    hsquared:::hs_build_model_spec(
      y ~ animal(1 | id, pedigree = ped) + common_env(age | litter),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    class = "hsquared_unsupported_syntax"
  )
})

test_that("hsquared() engine/target gating rejections carry the class", {
  ped <- hs_usc_ped()
  dat <- hs_usc_dat()

  # ML (REML = FALSE) on the default fit path (raised before any Julia check)
  expect_error(
    hsquared(
      y ~ animal(1 | id, pedigree = ped),
      data = dat,
      REML = FALSE
    ),
    class = "hsquared_unsupported_syntax"
  )

  # univariate response with target = "multivariate" (target/response mismatch)
  expect_error(
    hsquared(
      y ~ animal(1 | id, pedigree = ped),
      data = dat,
      control = hs_control(
        engine = "julia",
        engine_control = list(target = "multivariate")
      )
    ),
    class = "hsquared_unsupported_syntax"
  )
})

test_that("the julia-bridge target enumerator rejection carries the class", {
  ped <- hs_usc_ped()
  dat <- hs_usc_dat()

  expect_error(
    hsquared(
      y ~ animal(1 | id, pedigree = ped),
      data = dat,
      control = hs_control(
        engine = "julia",
        engine_control = list(target = "not_a_real_target")
      )
    ),
    class = "hsquared_unsupported_syntax"
  )
})

test_that("julia-bridge genetic_structure/rank engine_control rejections carry the class", {
  ped <- hs_usc_ped()
  dat <- hs_usc_dat()

  # invalid genetic_structure enumerator value (allowed-values rejection)
  expect_error(
    hsquared(
      y ~ animal(1 | id, pedigree = ped),
      data = dat,
      control = hs_control(
        engine = "julia",
        engine_control = list(genetic_structure = "not_a_structure")
      )
    ),
    class = "hsquared_unsupported_syntax"
  )

  # genetic_structure supplied without the multivariate target (planned gate)
  expect_error(
    hsquared(
      y ~ animal(1 | id, pedigree = ped),
      data = dat,
      control = hs_control(
        engine = "julia",
        engine_control = list(genetic_structure = "diagonal")
      )
    ),
    class = "hsquared_unsupported_syntax"
  )

  # structured lowrank / factor_analytic G0 (planned, not implemented)
  expect_error(
    hsquared(
      y ~ animal(1 | id, pedigree = ped),
      data = dat,
      control = hs_control(
        engine = "julia",
        engine_control = list(
          target = "multivariate",
          genetic_structure = "lowrank"
        )
      )
    ),
    class = "hsquared_unsupported_syntax"
  )

  # `rank` reserved for future structured-covariance controls
  expect_error(
    hsquared(
      y ~ animal(1 | id, pedigree = ped),
      data = dat,
      control = hs_control(
        engine = "julia",
        engine_control = list(
          target = "multivariate",
          genetic_structure = "diagonal",
          rank = 2
        )
      )
    ),
    class = "hsquared_unsupported_syntax"
  )
})

test_that("model_spec() opt-in-target preview rejection carries the class", {
  dat <- data.frame(y = c(1, 2, 3, 4), id = c("a", "b", "c", "d"))
  Ginv <- diag(4)
  dimnames(Ginv) <- list(dat$id, dat$id)

  expect_error(
    model_spec(
      y ~ genomic(1 | id, Ginv = Ginv),
      data = dat
    ),
    class = "hsquared_unsupported_syntax"
  )
})

test_that("data-validation guards do NOT carry the unsupported-syntax class", {
  ped <- hs_usc_ped()
  dat <- hs_usc_dat()

  # wrong-object guard: `control` not built by hs_control()
  cnd_control <- expect_error(
    hsquared(
      y ~ animal(1 | id, pedigree = ped),
      data = dat,
      control = "not a control"
    )
  )
  expect_false(inherits(cnd_control, "hsquared_unsupported_syntax"))

  # data-content guard that even mentions "not implemented" (missing response
  # values) still stays a plain error, NOT an unsupported-syntax condition.
  dat_na <- dat
  dat_na$y[[1L]] <- NA_real_
  cnd_na <- expect_error(
    hsquared:::hs_build_model_spec(
      y ~ animal(1 | id, pedigree = ped),
      data = dat_na,
      family = stats::gaussian(),
      REML = TRUE
    )
  )
  expect_false(inherits(cnd_na, "hsquared_unsupported_syntax"))
})
