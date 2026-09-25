# Honesty lock: R design docs that describe the R↔Julia bridge must keep
# engine-covered FA/SS rows cross-linked without promoting R-public status.
# No Julia required. Does not flip capability rows or public_covered_count.

test_that("bridge-gap design notes engine-covered V4-FA / V2-SSHINV without R flips", {
  path <- testthat::test_path("..", "..", "docs", "design", "19-on-main-bridge-gap.md")
  skip_if_not(file.exists(path), "design-19 not present in the check copy")
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")

  expect_match(text, "Update 2026-09-05", fixed = TRUE)
  expect_match(text, "V4-FA", fixed = TRUE)
  expect_match(text, "engine-covered", fixed = TRUE)
  expect_match(text, "V2-SSHINV", fixed = TRUE)
  expect_match(text, "public_covered_count` stays **7**", fixed = TRUE)
  expect_match(text, "12-bridge-compatibility.md", fixed = TRUE)
  expect_match(text, "itchyshin/HSquared.jl", fixed = TRUE)
  expect_match(text, "R FA grammar stays **planned**", fixed = TRUE)
  expect_match(text, "R surface stays opt-in partial", fixed = TRUE)
  expect_false(grepl("R FA grammar stays \\*\\*covered\\*\\*", text))
})

test_that("single-step bridge design keeps R partial after engine V2-SSHINV cover", {
  path <- testthat::test_path(
    "..",
    "..",
    "docs",
    "design",
    "25-single-step-construction-bridge.md"
  )
  skip_if_not(file.exists(path), "design-25 not present in the check copy")
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")

  expect_match(text, "V2-SSHINV", fixed = TRUE)
  expect_match(text, "engine-covered", fixed = TRUE)
  expect_match(text, "cf2a9bbf", fixed = TRUE)
  expect_match(text, "opt-in partial", fixed = TRUE)
  expect_match(text, "public_covered_count` stays **7**", fixed = TRUE)
  expect_match(text, "12-bridge-compatibility.md", fixed = TRUE)
  expect_false(grepl("mirrors the twin `V2-SSHINV` \\(partial\\)", text))
})

test_that("FA eigenbasis contract notes engine V4-FA cover without R promotion", {
  path <- testthat::test_path(
    "..",
    "..",
    "docs",
    "design",
    "29-structured-covariance-eigenbasis-bridge-contract.md"
  )
  skip_if_not(file.exists(path), "design-29 not present in the check copy")
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")

  expect_match(text, "Cross-link (2026-09-05)", fixed = TRUE)
  expect_match(text, "V4-FA", fixed = TRUE)
  expect_match(text, "engine-covered", fixed = TRUE)
  expect_match(text, "60895208", fixed = TRUE)
  expect_match(text, "stays **planned**", fixed = TRUE)
  expect_match(text, "public_covered_count` stays **7**", fixed = TRUE)
  expect_match(text, "12-bridge-compatibility.md", fixed = TRUE)
})
