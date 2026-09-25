# ?hsquared is a user front door: live animal-model path + limits link.
# It must not reopen as a genomic / Laplace / Willham audit page.

test_that("?hsquared description is a short live-path front door", {
  rd <- testthat::test_path("..", "..", "man", "hsquared.Rd")
  skip_if_not(file.exists(rd), "man/hsquared.Rd not present in the check copy")
  text <- paste(readLines(rd, warn = FALSE), collapse = "\n")

  desc <- sub("^.*\\\\description\\{", "", text)
  desc <- sub("\\}\\s*\\\\details\\{.*$", "", desc)

  expect_match(desc, "animal(1 | id", fixed = TRUE)
  expect_match(desc, "validate", fixed = TRUE)
  expect_match(text, "current-limits", fixed = TRUE)
  expect_match(text, "Can I fit and report this", fixed = TRUE)
  expect_match(text, "validation_status", fixed = TRUE)

  expect_no_match(text, "genomic_variance_ratio", fixed = TRUE)
  expect_no_match(text, "VanRaden", fixed = TRUE)
  expect_no_match(text, "Laplace", ignore.case = TRUE)
  expect_no_match(text, "Willham", fixed = TRUE)
  expect_no_match(text, "first planned v0.1", fixed = TRUE)
  expect_no_match(text, "The v0.1 parser accepts only", fixed = TRUE)

  collapsed <- gsub("\\s+", " ", desc)
  expect_lt(nchar(collapsed), 800L)
})
