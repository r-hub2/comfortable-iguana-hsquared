# plot.hsquared_fit: base-graphics diagnostics. Tests draw to a null device and
# assert the call runs without error and returns the fit invisibly.

make_fit <- function(result, y = NULL) {
  hsquared:::hs_new_fit(
    call = quote(hsquared(y ~ animal(1 | id, pedigree = ped), data = dat)),
    spec = list(method = "REML", family = list(family = "gaussian")),
    payload = list(y = y %||% seq_len(10)),
    result = result
  )
}
`%||%` <- function(a, b) if (is.null(a)) b else a

with_null_device <- function(code) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(code)
}

test_that("plot(type = 'variance') draws with and without SE whiskers", {
  fit_se <- make_fit(list(
    variance_components = data.frame(
      component = c("animal", "residual"),
      estimate = c(0.4, 0.6)
    ),
    variance_component_se = data.frame(
      component = c("animal", "residual"),
      se = c(0.1, 0.12),
      stringsAsFactors = FALSE
    )
  ))
  with_null_device({
    out <- plot(fit_se, type = "variance")
    expect_identical(out, fit_se)
  })

  fit_plain <- make_fit(list(
    variance_components = data.frame(
      component = c("animal", "residual"),
      estimate = c(0.4, 0.6)
    )
  ))
  with_null_device(expect_invisible(plot(fit_plain, type = "variance")))
})

test_that("plot(type = 'residuals') draws when fitted + response exist", {
  fit <- make_fit(
    list(predictions = data.frame(.fitted = seq_len(10) + 0.5)),
    y = seq_len(10)
  )
  with_null_device(expect_invisible(plot(fit, type = "residuals")))
})

test_that("plot() errors clearly when the needed fields are absent", {
  with_null_device({
    expect_error(
      plot(make_fit(list(heritability = 0.4)), type = "variance"),
      "no variance-component estimates"
    )
    expect_error(
      plot(
        make_fit(list(
          variance_components = data.frame(
            component = "animal",
            estimate = 0.4
          )
        )),
        type = "residuals"
      ),
      "no fitted values"
    )
  })
})
