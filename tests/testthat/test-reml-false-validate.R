# R3 / Pat P4: REML = FALSE must not succeed on engine = "validate"
# and still describe a REML contract. One rule with the default fit
# path. model_spec() / internal builders stay able to label a spec ML.
# No covered flip.

hs_r3_ped <- function() {
  data.frame(
    id = c("a", "b", "c"),
    sire = c(NA, NA, "a"),
    dam = c(NA, NA, "b")
  )
}

hs_r3_dat <- function() {
  data.frame(
    y = c(1, 2, 3),
    sex = c("m", "f", "m"),
    id = c("a", "b", "c")
  )
}

test_that("engine = validate rejects REML = FALSE with the REML = TRUE path", {
  ped <- hs_r3_ped()
  dat <- hs_r3_dat()

  expect_error(
    hsquared(
      y ~ sex + animal(1 | id, pedigree = ped),
      data = dat,
      family = stats::gaussian(),
      REML = FALSE,
      control = hs_control(engine = "validate")
    ),
    "Closest working call: `REML = TRUE`",
    fixed = TRUE
  )
})

test_that("engine = validate REML = FALSE is unsupported syntax, not a validated spec", {
  ped <- hs_r3_ped()
  dat <- hs_r3_dat()

  err <- expect_error(
    hsquared(
      y ~ sex + animal(1 | id, pedigree = ped),
      data = dat,
      REML = FALSE,
      control = hs_control(engine = "validate")
    ),
    class = "hsquared_unsupported_syntax"
  )
  msg <- conditionMessage(err)
  expect_match(msg, "Closest working call: `REML = TRUE`", fixed = TRUE)
  # Test of the test: the old success path confirmed the contract.
  expect_false(grepl("Validated the v0.1", msg, fixed = TRUE))
})

test_that("default fit path REML = FALSE names the same nearest path", {
  ped <- hs_r3_ped()
  dat <- hs_r3_dat()

  expect_error(
    hsquared(
      y ~ sex + animal(1 | id, pedigree = ped),
      data = dat,
      family = stats::gaussian(),
      REML = FALSE
    ),
    "Closest working call: `REML = TRUE`",
    fixed = TRUE
  )
})

test_that("validate still accepts REML = TRUE", {
  ped <- hs_r3_ped()
  dat <- hs_r3_dat()

  expect_message(
    spec <- hsquared(
      y ~ sex + animal(1 | id, pedigree = ped),
      data = dat,
      REML = TRUE,
      control = hs_control(engine = "validate")
    ),
    "Validated the v0.1 animal-model contract",
    fixed = TRUE
  )
  expect_type(spec, "list")
  expect_equal(spec$method, "REML")
})

test_that("model_spec and internal builders may still label an ML spec", {
  ped <- hs_r3_ped()
  dat <- hs_r3_dat()

  spec <- hsquared:::hs_build_model_spec(
    y ~ sex + animal(1 | id, pedigree = ped),
    data = dat,
    family = stats::gaussian(),
    REML = FALSE
  )
  expect_equal(spec$method, "ML")

  preview <- model_spec(
    y ~ sex + animal(1 | id, pedigree = ped),
    data = dat,
    family = stats::gaussian(),
    REML = FALSE
  )
  expect_s3_class(preview, "hs_model_spec")
  expect_equal(preview$method, "ML")
})
