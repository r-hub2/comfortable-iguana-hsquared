test_that("hs_require_suggests is a no-op when the package is installed", {
  expect_true(hs_require_suggests("testthat", gate = "unit"))
})

test_that("hs_require_suggests fails loudly under NOT_CRAN when Suggests missing", {
  withr::local_envvar(c(NOT_CRAN = "true"))
  expect_error(
    hs_require_suggests("hsquaredFakeSuggestsMissingPkg", gate = "MV-1"),
    "MV-1 requires Suggests package 'hsquaredFakeSuggestsMissingPkg'",
    fixed = TRUE
  )
})

test_that("hs_require_suggests skips quietly when NOT_CRAN is unset and Suggests missing", {
  withr::local_envvar(c(NOT_CRAN = "false"))
  cnd <- tryCatch(
    hs_require_suggests("hsquaredFakeSuggestsMissingPkg", gate = "MV-1"),
    skip = identity
  )
  expect_s3_class(cnd, "skip")
})
