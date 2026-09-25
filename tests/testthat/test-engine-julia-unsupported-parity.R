# Parity-09 unsupported-syntax regressions (engine = "julia").
# Closes SILENT-GAPs SG1–SG3 from docs/design/55-r-julia-engine-julia-parity-09.md
# and docs/design/56-r-julia-parity-09-error-audit.md. No live Julia required.

hs_parity09_tiny_ped <- function() {
  data.frame(
    id = c("a", "b", "c"),
    sire = c(NA, NA, "a"),
    dam = c(NA, NA, "b"),
    stringsAsFactors = FALSE
  )
}

test_that("SG1: cbind + permanent names repeatability and MV-without-PE tips (D3)", {
  ped <- hs_parity09_tiny_ped()
  dat <- data.frame(
    y1 = c(1, 2, 3, 1.5, 2.5),
    y2 = c(2, 1, 2.5, 1.2, 2.2),
    id = c("a", "b", "c", "a", "b"),
    stringsAsFactors = FALSE
  )

  err <- tryCatch(
    hsquared:::hs_build_model_spec(
      cbind(y1, y2) ~ animal(1 | id, pedigree = ped) + permanent(1 | id),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    error = function(e) e
  )
  expect_s3_class(err, "hsquared_unsupported_syntax")
  msg <- conditionMessage(err)
  expect_match(msg, "hsquared#237", fixed = TRUE)
  expect_match(msg, "target = \"repeatability\"", fixed = TRUE)
  expect_match(msg, "drop `permanent()`", fixed = TRUE)
  expect_match(msg, "public_covered_count stays 7", fixed = TRUE)
  expect_match(msg, "not implemented", fixed = TRUE)
})

test_that("SG2: matrix_free target is engine-only named abort (D2)", {
  for (tgt in c("matrix_free", "matrix_free_reml", "matfree")) {
    err <- tryCatch(
      hsquared:::hs_validate_julia_target(tgt),
      error = function(e) e
    )
    expect_s3_class(err, "hsquared_unsupported_syntax")
    msg <- conditionMessage(err)
    expect_match(msg, "engine-only", fixed = TRUE)
    expect_match(msg, "fit_matrix_free_reml", fixed = TRUE)
    expect_match(msg, "scale_method = \"auto\"", fixed = TRUE)
    expect_match(msg, "public_covered_count stays 7", fixed = TRUE)
  }
})

test_that("D1: factor_analytic and lowrank stay named planned errors", {
  for (gs in c("factor_analytic", "lowrank")) {
    err <- tryCatch(
      hsquared:::hs_validate_genetic_structure_control(
        hs_control(
          engine = "julia",
          engine_control = list(
            target = "multivariate",
            genetic_structure = gs
          )
        ),
        "multivariate"
      ),
      error = function(e) e
    )
    expect_s3_class(err, "hsquared_unsupported_syntax")
    msg <- conditionMessage(err)
    if (identical(gs, "factor_analytic")) {
      expect_match(msg, "factor_analytic", fixed = TRUE)
      expect_match(msg, "not an R-public factor-analytic fit", fixed = TRUE)
    } else {
      expect_match(msg, "lowrank", fixed = TRUE)
      expect_match(msg, "planned", fixed = TRUE)
    }
  }
})

test_that("SG3: RR + permanent remains a named planned reject", {
  ped <- hs_parity09_tiny_ped()
  dat <- data.frame(
    y = c(1, 2, 3, 1.5, 2.5),
    id = c("a", "b", "c", "a", "b"),
    age = c(1, 1, 1, 2, 2),
    stringsAsFactors = FALSE
  )
  err <- tryCatch(
    hsquared:::hs_build_model_spec(
      y ~ animal(rr(age, order = 2) | id, pedigree = ped) + permanent(1 | id),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    error = function(e) e
  )
  expect_s3_class(err, "hsquared_unsupported_syntax")
  expect_match(
    conditionMessage(err),
    "planned, not implemented",
    fixed = TRUE
  )
})
