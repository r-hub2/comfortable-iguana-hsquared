test_that("hs_skip_live_julia skips on CRAN lane unless HSQUARED_JULIA_TESTS=true", {
  withr::local_envvar(c(NOT_CRAN = "false", HSQUARED_JULIA_TESTS = NA))
  cnd <- tryCatch(hs_skip_live_julia(), skip = identity)
  expect_s3_class(cnd, "skip")

  withr::local_envvar(c(NOT_CRAN = "false", HSQUARED_JULIA_TESTS = "true"))
  expect_silent(hs_skip_live_julia())

  withr::local_envvar(c(NOT_CRAN = "true", HSQUARED_JULIA_TESTS = NA))
  expect_silent(hs_skip_live_julia())
})
