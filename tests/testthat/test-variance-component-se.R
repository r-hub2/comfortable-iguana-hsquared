# Experimental variance-component / heritability standard errors (engine row
# V1-HERIT-CI, partial). R-side extractor + normalizer tested with fixtures; the
# live engine path is exercised opportunistically by the bridge when Julia is
# present.

test_that("SE extractors return the fields when present", {
  result <- list(
    variance_component_se = data.frame(
      component = c("animal", "residual"),
      se = c(0.12, 0.18),
      stringsAsFactors = FALSE
    ),
    heritability = data.frame(
      term = "animal",
      estimate = 0.4,
      stringsAsFactors = FALSE
    ),
    heritability_se = 0.07
  )
  fit <- hsquared:::hs_new_fit(
    call = quote(hsquared(y ~ animal(1 | id, pedigree = ped), data = dat)),
    spec = list(method = "REML", family = list(family = "gaussian")),
    payload = list(y = seq_len(10)),
    result = result
  )

  expect_s3_class(fit, "hsquared_fit")
  expect_equal(
    variance_component_standard_errors(fit),
    result$variance_component_se
  )
  # Public shape is data.frame(term, se) so it merges with heritability()
  # (hsquared#236). Internal storage may still be a bare numeric.
  h2se <- heritability_standard_error(fit)
  expect_s3_class(h2se, "data.frame")
  expect_equal(names(h2se), c("term", "se"))
  expect_equal(h2se$term, "animal")
  expect_equal(h2se$se, 0.07)
  expect_equal(
    merge(heritability(fit), h2se, by = "term")$se,
    0.07
  )
})

test_that("SE extractors error clearly without the field or object", {
  expect_error(
    variance_component_standard_errors(seq_len(10)),
    "requires an `hsquared_fit` object"
  )
  expect_error(
    heritability_standard_error(seq_len(10)),
    "requires an `hsquared_fit` object"
  )

  fit_no_se <- hsquared:::hs_new_fit(
    call = quote(hsquared(y ~ animal(1 | id, pedigree = ped), data = dat)),
    spec = list(method = "REML", family = list(family = "gaussian")),
    payload = list(y = seq_len(10)),
    result = list(heritability = data.frame(term = "animal", estimate = 0.4))
  )
  expect_error(
    variance_component_standard_errors(fit_no_se),
    "variance-component standard errors"
  )
  expect_error(
    heritability_standard_error(fit_no_se),
    "heritability standard error"
  )
})

test_that("hs_normalize_variance_component_se() builds a component/se data frame", {
  out <- hsquared:::hs_normalize_variance_component_se(
    list(sigma_a2 = 0.12, sigma_e2 = 0.18)
  )
  expect_equal(out$component, c("animal", "residual"))
  expect_equal(out$se, c(0.12, 0.18))
})

test_that("heritability_standard_error accepts a term/se frame already on the fit", {
  result <- list(
    heritability = data.frame(
      term = "animal",
      estimate = 0.35,
      stringsAsFactors = FALSE
    ),
    heritability_se = data.frame(
      term = "animal",
      se = 0.05,
      stringsAsFactors = FALSE
    )
  )
  fit <- hsquared:::hs_new_fit(
    call = quote(hsquared(y ~ animal(1 | id, pedigree = ped), data = dat)),
    spec = list(method = "REML", family = list(family = "gaussian")),
    payload = list(y = seq_len(10)),
    result = result
  )
  expect_equal(heritability_standard_error(fit), result$heritability_se)
})
