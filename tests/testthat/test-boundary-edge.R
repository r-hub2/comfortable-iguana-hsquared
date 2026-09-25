# Edge cases for the variance-component boundary surfacing:
# (#7) a NEGATIVE minimum share must be reported as an inadmissible negative
#      variance, distinct from the benign at/near-zero boundary; and
# (#8) a malformed (atomic-vector) `variance_components` must degrade without
#      crashing summary() / fit_diagnostics().

test_that("a negative animal variance is reported as inadmissible, not at/near zero", {
  negative <- hsquared:::hs_new_fit(
    spec = list(
      method = "REML",
      family = list(family = "gaussian"),
      target = "univariate"
    ),
    payload = list(y = 1:3),
    result = list(
      variance_components = data.frame(
        component = c("animal", "residual"),
        estimate = c(-0.2, 1.0)
      ),
      converged = TRUE
    )
  )

  # summary() print wording: negative / inadmissible, NOT "at/near zero".
  summ <- summary(negative)
  expect_identical(summ$at_boundary_class, "negative")
  expect_true(isTRUE(summ$at_boundary))
  printed <- paste(capture.output(print(summ)), collapse = "\n")
  expect_match(printed, "NEGATIVE")
  expect_match(printed, "inadmissible")
  expect_false(grepl("at/near zero", printed))

  # fit_diagnostics() condition row carries the same distinction.
  diag <- fit_diagnostics(negative)
  expect_equal(diag$value[diag$metric == "at_boundary"], "TRUE")
  expect_equal(
    diag$value[diag$metric == "at_boundary_condition"],
    "negative (inadmissible variance)"
  )
})

test_that("an atomic-vector variance_components degrades without crashing", {
  malformed <- hsquared:::hs_new_fit(
    spec = list(
      method = "REML",
      family = list(family = "gaussian"),
      target = "univariate"
    ),
    payload = list(y = 1:3),
    result = list(
      # Wrong shape: an atomic numeric vector, not a list/data.frame. The old
      # code crashed here with "$ operator is invalid for atomic vectors".
      variance_components = c(animal = 0.5, residual = 0.5),
      converged = TRUE
    )
  )

  expect_silent(summ <- summary(malformed))
  expect_null(summ$at_boundary)
  expect_null(summ$at_boundary_class)
  expect_silent(capture.output(print(summ)))

  expect_silent(diag <- fit_diagnostics(malformed))
  # Unavailable -> NULL contract: the boundary rows are dropped, not present.
  expect_false("at_boundary" %in% diag$metric)
  expect_false("at_boundary_condition" %in% diag$metric)
})

test_that("a matrix variance_components also degrades without crashing", {
  malformed <- hsquared:::hs_new_fit(
    spec = list(
      method = "REML",
      family = list(family = "gaussian"),
      target = "univariate"
    ),
    payload = list(y = 1:3),
    result = list(
      variance_components = matrix(c(0.5, 0.5), ncol = 1),
      converged = TRUE
    )
  )

  expect_silent(diag <- fit_diagnostics(malformed))
  expect_false("at_boundary" %in% diag$metric)
})

test_that("a positive near-zero animal variance still flags the zero boundary", {
  boundary <- hsquared:::hs_new_fit(
    spec = list(
      method = "REML",
      family = list(family = "gaussian"),
      target = "univariate"
    ),
    payload = list(y = 1:3),
    result = list(
      variance_components = data.frame(
        component = c("animal", "residual"),
        estimate = c(1e-8, 1.0)
      ),
      converged = TRUE
    )
  )

  summ <- summary(boundary)
  expect_identical(summ$at_boundary_class, "zero")
  expect_true(isTRUE(summ$at_boundary))
  printed <- paste(capture.output(print(summ)), collapse = "\n")
  expect_match(printed, "at/near zero")
  expect_false(grepl("NEGATIVE", printed))

  diag <- fit_diagnostics(boundary)
  expect_equal(diag$value[diag$metric == "at_boundary"], "TRUE")
  expect_equal(
    diag$value[diag$metric == "at_boundary_condition"],
    "at/near zero (clean boundary)"
  )
})
