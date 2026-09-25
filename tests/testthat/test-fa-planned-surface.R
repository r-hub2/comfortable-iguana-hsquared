test_that("the R FA surface remains planned while the engine is covered", {
  status <- validation_status()

  # The public count is a claim-surface fence, not a reason to promote the
  # engine's V4-FA evidence into the R package.
  capability_path <- testthat::test_path(
    "..",
    "..",
    "docs",
    "design",
    "capability-status.md"
  )
  testthat::skip_if_not(
    file.exists(capability_path),
    "repository-only capability-status.md is unavailable in the build tarball"
  )
  capability_doc <- paste(
    readLines(capability_path, warn = FALSE),
    collapse = "\n"
  )
  expect_match(
    capability_doc,
    "factor-analytic G matrices | planned",
    fixed = TRUE
  )
  expect_true(grepl("public_covered_count[^\\n]*7", capability_doc))

  formulas <- formula_status()
  fa_term <- "animal(trait | id, pedigree = ped, cov = fa(K = 2))"
  fa_note <- formulas$current_behavior[formulas$term == fa_term]
  expect_length(fa_note, 1L)
  expect_match(fa_note, "Julia V4-FA is engine-covered", fixed = TRUE)
  expect_match(fa_note, "Not an R-public FA claim", fixed = TRUE)
  expect_match(fa_note, "public_covered_count stays 7", fixed = TRUE)

  # Live claim surfaces must track experimental 0.9.0 (not a stale 0.8.0 pin).
  expect_match(
    capability_doc,
    "Experimental stays **0.9.0**. Not an R-public FA flip.",
    fixed = TRUE
  )
  expect_match(
    capability_doc,
    "Experimental stays **0.9.0**. Not an R-public SS flip.",
    fixed = TRUE
  )
  expect_false(grepl(
    "Experimental stays \\*\\*0\\.8\\.0\\*\\*\\. Not an R-public FA flip",
    capability_doc
  ))
  expect_false(grepl(
    "Experimental stays \\*\\*0\\.8\\.0\\*\\*\\. Not an R-public SS flip",
    capability_doc
  ))
  expect_match(capability_doc, "Post-#366 T2", fixed = TRUE)
})

test_that("factor-analytic controls reject with the R boundary", {
  expect_error(
    hsquared:::hs_validate_genetic_structure_control(
      hs_control(
        engine = "julia",
        engine_control = list(
          target = "multivariate",
          genetic_structure = "factor_analytic"
        )
      ),
      "multivariate"
    ),
    "not an R-public factor-analytic fit",
    fixed = TRUE
  )
  expect_error(
    hsquared:::hs_validate_genetic_structure_control(
      hs_control(
        engine = "julia",
        engine_control = list(
          target = "multivariate",
          genetic_structure = "factor_analytic"
        )
      ),
      "multivariate"
    ),
    "public_covered_count stays 7",
    fixed = TRUE
  )
  expect_error(
    hsquared:::hs_validate_genetic_structure_control(
      hs_control(
        engine = "julia",
        engine_control = list(
          target = "multivariate",
          genetic_structure = "lowrank"
        )
      ),
      "multivariate"
    ),
    "public_covered_count stays 7",
    fixed = TRUE
  )
})

test_that("FA formula grammar rejects before model construction", {
  ped <- data.frame(
    id = c("a", "b"),
    sire = c(NA, NA),
    dam = c(NA, NA)
  )
  dat <- data.frame(y = c(1, 2), id = ped$id)

  expect_error(
    hsquared:::hs_build_model_spec(
      y ~ animal(1 | id, pedigree = ped, cov = fa(K = 2)),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "not an R-public fit",
    fixed = TRUE
  )
})
