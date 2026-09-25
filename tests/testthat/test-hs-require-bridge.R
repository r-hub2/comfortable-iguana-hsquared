test_that("hs_require_bridge is a no-op when the bridge is available", {
  skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "Live Julia bridge required for the positive-path unit."
  )
  withr::local_envvar(c(
    NOT_CRAN = "true",
    HSQUARED_JULIA_TESTS = "true",
    HSQUARED_REQUIRE_BRIDGE = "true"
  ))
  expect_true(hs_require_bridge("unit"))
})

test_that("hs_require_bridge fails loudly when REQUIRE_BRIDGE and bridge missing", {
  missing_project <- tempfile("hsquared-no-julia-project-")
  dir.create(missing_project)
  withr::defer(unlink(missing_project, recursive = TRUE))
  withr::local_envvar(c(
    NOT_CRAN = "true",
    HSQUARED_JULIA_TESTS = "true",
    HSQUARED_REQUIRE_BRIDGE = "true",
    HSQUARED_JULIA_PROJECT = missing_project
  ))
  expect_error(
    hs_require_bridge("A26 parity"),
    "HSQUARED_REQUIRE_BRIDGE=true was set",
    fixed = TRUE
  )
})

test_that("hs_require_bridge skips quietly when REQUIRE_BRIDGE is unset and bridge missing", {
  missing_project <- tempfile("hsquared-no-julia-project-")
  dir.create(missing_project)
  withr::defer(unlink(missing_project, recursive = TRUE))
  withr::local_envvar(c(
    NOT_CRAN = "true",
    HSQUARED_JULIA_TESTS = "true",
    HSQUARED_REQUIRE_BRIDGE = "false",
    HSQUARED_JULIA_PROJECT = missing_project
  ))
  cnd <- tryCatch(
    hs_require_bridge("A26 parity"),
    skip = identity
  )
  expect_s3_class(cnd, "skip")
})
