test_that("package metadata is available", {
  expect_equal(utils::packageDescription("hsquared")$Package, "hsquared")
})
