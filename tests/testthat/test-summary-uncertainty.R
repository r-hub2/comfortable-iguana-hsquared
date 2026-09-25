# summary()/print() surfaces the experimental uncertainty fields (CIs / SEs)
# when an hsquared_fit carries them, with honest experimental labelling.

make_fit <- function(result) {
  hsquared:::hs_new_fit(
    call = quote(hsquared(y ~ animal(1 | id, pedigree = ped), data = dat)),
    spec = list(method = "REML", family = list(family = "gaussian")),
    payload = list(y = seq_len(10)),
    result = result
  )
}

test_that("summary() carries the experimental uncertainty fields when present", {
  fit <- make_fit(list(
    heritability = data.frame(term = "animal", estimate = 0.42),
    heritability_se = 0.07,
    heritability_interval = data.frame(
      estimate = 0.42,
      lower = 0.21,
      upper = 0.66,
      level = 0.95,
      se = 0.07,
      method = "delta",
      stringsAsFactors = FALSE
    ),
    variance_component_se = data.frame(
      component = c("animal", "residual"),
      se = c(0.12, 0.18),
      stringsAsFactors = FALSE
    )
  ))
  s <- summary(fit)
  expect_equal(s$heritability_se, 0.07)
  expect_equal(s$heritability_interval$lower, 0.21)
  expect_equal(s$variance_component_se$se, c(0.12, 0.18))

  out <- paste(capture.output(print(s)), collapse = "\n")
  expect_match(out, "experimental")
  expect_match(out, "95% CI")
  expect_match(out, "SE")
})

test_that("print(summary()) omits the uncertainty block when fields are absent", {
  fit <- make_fit(list(
    heritability = data.frame(term = "animal", estimate = 0.4),
    converged = TRUE
  ))
  out <- paste(capture.output(print(summary(fit))), collapse = "\n")
  expect_false(grepl("uncertainty", out))
})

test_that("print(summary()) shows the experimental repeatability interval", {
  fit <- make_fit(list(
    repeatability_interval = data.frame(
      estimate = 0.55,
      lower = 0.34,
      upper = 0.74,
      level = 0.95,
      se = 0.10,
      stringsAsFactors = FALSE
    )
  ))
  out <- paste(capture.output(print(summary(fit))), collapse = "\n")
  expect_match(out, "repeatability uncertainty \\(experimental")
  expect_match(out, "\\[0.34, 0.74\\]")
})
