# Experimental repeatability-coefficient CI (engine row V3-REPEAT-REML, partial).
# R-side extractor + normalizer tested with fixtures; the live engine path is
# exercised opportunistically by the repeatability bridge when Julia is present.

test_that("repeatability_interval() returns the interval when present", {
  result <- list(
    repeatability_interval = data.frame(
      estimate = 0.55,
      lower = 0.34,
      upper = 0.74,
      level = 0.95,
      se = 0.10,
      stringsAsFactors = FALSE
    )
  )
  fit <- hsquared:::hs_new_fit(
    call = quote(hsquared(
      y ~ animal(1 | id, pedigree = ped) + permanent(1 | id),
      data = dat
    )),
    spec = list(method = "REML", family = list(family = "gaussian")),
    payload = list(y = seq_len(10)),
    result = result
  )

  expect_s3_class(fit, "hsquared_fit")
  expect_equal(repeatability_interval(fit), result$repeatability_interval)
})

test_that("repeatability_interval() errors clearly without the field or object", {
  expect_error(
    repeatability_interval(seq_len(10)),
    "requires an `hsquared_fit` object"
  )

  fit_no_ri <- hsquared:::hs_new_fit(
    call = quote(hsquared(
      y ~ animal(1 | id, pedigree = ped) + permanent(1 | id),
      data = dat
    )),
    spec = list(method = "REML", family = list(family = "gaussian")),
    payload = list(y = seq_len(10)),
    result = list(repeatability = 0.55)
  )
  expect_error(
    repeatability_interval(fit_no_ri),
    "experimental repeatability confidence interval"
  )
})

test_that("hs_normalize_repeatability_interval() builds a one-row data frame", {
  out <- hsquared:::hs_normalize_repeatability_interval(
    list(
      repeatability = 0.55,
      lower = 0.34,
      upper = 0.74,
      level = 0.95,
      se = 0.10
    )
  )
  expect_equal(nrow(out), 1L)
  expect_equal(out$estimate, 0.55)
  expect_true(out$lower <= out$estimate && out$estimate <= out$upper)
  expect_equal(out$se, 0.10)
  # The repeatability CI is for t, not h2; there is no `method` column.
  expect_false("method" %in% names(out))
})
