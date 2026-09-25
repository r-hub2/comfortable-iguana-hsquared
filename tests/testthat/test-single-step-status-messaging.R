test_that("formula_status distinguishes the two single-step on-ramps", {
  status <- formula_status()

  supplied <- status[
    status$term == "single_step(1 | id, Hinv = Hinv)",
    ,
    drop = FALSE
  ]
  expect_length(supplied$current_behavior, 1L)
  expect_match(supplied$current_behavior, "supplied-inverse", fixed = TRUE)
  expect_match(
    supplied$current_behavior,
    "does not construct H^-1 from pedigree and markers",
    fixed = TRUE
  )
  expect_match(
    supplied$current_behavior,
    "target = \"single_step_construct\"",
    fixed = TRUE
  )

  constructed <- status[
    status$term == "single_step(1 | id, pedigree = ped, markers = M)",
    ,
    drop = FALSE
  ]
  expect_length(constructed$current_behavior, 1L)
  expect_match(
    constructed$current_behavior,
    "separate from the supplied-`Hinv` route",
    fixed = TRUE
  )
  expect_match(
    constructed$current_behavior,
    "not a default-route promotion",
    fixed = TRUE
  )
  expect_match(
    constructed$current_behavior,
    "public_covered_count stays 7",
    fixed = TRUE
  )

  bundle <- status[
    grepl("^single_step\\(1 \\| id\\) with data =", status$term),
    ,
    drop = FALSE
  ]
  expect_length(bundle$current_behavior, 1L)
  expect_match(
    bundle$current_behavior,
    "plain data.frame still needs an explicit pedigree and marker matrix",
    fixed = TRUE
  )
})

test_that("single-step construction errors keep the two on-ramps distinct", {
  dat <- data.frame(y = c(1, 2), id = c("a", "b"))

  expect_error(
    hsquared:::hs_build_model_spec(
      y ~ single_step(1 | id),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "or `pedigree = ped` + `markers = M`",
    fixed = TRUE
  )

  expect_error(
    hsquared:::hs_build_model_spec(
      y ~ single_step(
        1 | id,
        pedigree = data.frame(
          id = c("a", "b"),
          sire = c(NA, NA),
          dam = c(NA, NA)
        )
      ),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "Supply a precomputed `Hinv` instead via",
    fixed = TRUE
  )
})
