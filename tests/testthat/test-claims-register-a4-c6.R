# A29 A4 + A21 C6 leftover: claims-register wording only.
# docs/ is .Rbuildignore'd, so this is source-tree-only (skipped under R CMD check).

test_that("claims register states covered multivariate scope after G10", {
  path <- file.path(
    testthat::test_path("..", ".."),
    "docs",
    "design",
    "06-public-claims-register.md"
  )
  skip_if_not(
    file.exists(path),
    "06-public-claims-register.md is not in the build tarball"
  )
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")

  expect_true(
    grepl("Covered numeric claim is scoped to k=2 unstructured", text, fixed = TRUE),
    info = "A4: post-G10 covered-scope sentence is unconditional"
  )
  expect_false(
    grepl("Covered numeric claim (when flipped) is scoped", text, fixed = TRUE),
    info = "A4: pre-flip (when flipped) conditional is gone"
  )
  expect_true(
    grepl("| multivariate Gaussian animal model | covered (validation-scale) |", text, fixed = TRUE) ||
      grepl("multivariate Gaussian animal model | covered (validation-scale)", text, fixed = TRUE),
    info = "A4: claims register status cell is covered"
  )
})

test_that("claims register distinguishes the two maternal m2 denominators", {
  path <- file.path(
    testthat::test_path("..", ".."),
    "docs",
    "design",
    "06-public-claims-register.md"
  )
  skip_if_not(
    file.exists(path),
    "06-public-claims-register.md is not in the build tarball"
  )
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")

  expect_true(
    grepl("maternal_proportion()", text, fixed = TRUE),
    info = "C6: two-effect fence names maternal_proportion()"
  )
  expect_true(
    grepl("no covariance", text, fixed = TRUE),
    info = "C6: two-effect fence names the no-covariance m2 denominator"
  )
  expect_true(
    grepl("omits the covariance", text, fixed = TRUE),
    info = "C6: direct-maternal fence names the other m2 denominator"
  )
  expect_true(
    grepl("m2_maternal", text, fixed = TRUE),
    info = "C6: direct-maternal fence names m2_maternal"
  )
})
