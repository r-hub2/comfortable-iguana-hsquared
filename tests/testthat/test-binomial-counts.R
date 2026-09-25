# Binomial count responses: cbind(successes, failures) ~ ... with family =
# binomial() must be a non-Gaussian binomial-trials model, NOT silently coerced
# into a 2-trait multivariate Gaussian (the family-blind cbind bug) or down to a
# binary Bernoulli (dropping the trial counts). Per-record trial counts
# (successes + failures) are parsed and carried as a per-record integer vector
# that may vary across records; a constant vector is the common-trial special
# case and an all-ones vector reduces to Bernoulli. Live varying-trial fitting is
# verified through the engine round-trip ([live] tests below); the pure-R tests
# cover the R-side parsing, payload carry, and family classification.

ped4 <- function() {
  data.frame(
    id = c("s", "d", "a", "b"),
    sire = c(NA, NA, "s", "s"),
    dam = c(NA, NA, "d", "d"),
    stringsAsFactors = FALSE
  )
}

test_that("cbind(successes, failures) + binomial() builds a binomial-counts spec, not multivariate", {
  ped <- ped4()
  # equal row totals (3 trials each) -> a constant per-record vector c(3,3,3,3)
  dat <- data.frame(
    succ = c(1, 2, 3, 0),
    fail = c(2, 1, 0, 3),
    id = c("s", "d", "a", "b")
  )
  spec <- hsquared:::hs_build_model_spec(
    cbind(succ, fail) ~ animal(1 | id, pedigree = ped),
    data = dat,
    family = stats::binomial(),
    REML = TRUE,
    allow_families = c("gaussian", "poisson", "binomial")
  )
  # NOT multivariate; it is a single-response binomial-counts model
  expect_false(isTRUE(spec$response$multivariate))
  expect_true(isTRUE(spec$response$binomial_counts))
  # per-record vector (constant here); equal totals are the common-trial case
  expect_identical(spec$response$n_trials, rep(3L, 4))
  # the response values are the success counts
  expect_equal(as.numeric(spec$response$values), c(1, 2, 3, 0))
})

test_that("cbind(successes, failures) + binomial() with varying row totals builds a per-record vector", {
  ped <- ped4()
  dat <- data.frame(
    succ = c(1, 2, 3, 0),
    fail = c(2, 1, 1, 3), # totals 3,3,4,3 -> per-record varying
    id = c("s", "d", "a", "b")
  )
  spec <- hsquared:::hs_build_model_spec(
    cbind(succ, fail) ~ animal(1 | id, pedigree = ped),
    data = dat,
    family = stats::binomial(),
    REML = TRUE,
    allow_families = c("gaussian", "poisson", "binomial")
  )
  expect_true(isTRUE(spec$response$binomial_counts))
  expect_false(isTRUE(spec$response$multivariate))
  # per-record trial counts: an order-preserving integer vector, not collapsed
  expect_identical(spec$response$n_trials, c(3L, 3L, 4L, 3L))
  expect_equal(as.numeric(spec$response$values), c(1, 2, 3, 0))
})

test_that("cbind(t1, t2) + gaussian() is still multivariate (no regression)", {
  ped <- ped4()
  dat <- data.frame(
    t1 = c(1.1, 2.2, 3.0, 2.5),
    t2 = c(0.5, 1.5, 2.0, 1.0),
    id = c("s", "d", "a", "b")
  )
  spec <- hsquared:::hs_build_model_spec(
    cbind(t1, t2) ~ animal(1 | id, pedigree = ped),
    data = dat,
    family = stats::gaussian(),
    REML = TRUE
  )
  expect_true(isTRUE(spec$response$multivariate))
})

test_that("a binary 0/1 binomial response stays Bernoulli (vector, not counts)", {
  ped <- ped4()
  dat <- data.frame(y = c(0, 1, 1, 0), id = c("s", "d", "a", "b"))
  spec <- hsquared:::hs_build_model_spec(
    y ~ animal(1 | id, pedigree = ped),
    data = dat,
    family = stats::binomial(),
    REML = TRUE,
    allow_families = c("gaussian", "poisson", "binomial")
  )
  expect_false(isTRUE(spec$response$multivariate))
  expect_false(isTRUE(spec$response$binomial_counts))
})

test_that("the family-symbol mapper distinguishes Bernoulli from Binomial(n_trials)", {
  expect_equal(
    hsquared:::hs_nongaussian_family_symbol(stats::binomial()),
    "bernoulli"
  )
  expect_equal(
    hsquared:::hs_nongaussian_family_symbol(stats::binomial(), n_trials = 1L),
    "bernoulli"
  )
  expect_equal(
    hsquared:::hs_nongaussian_family_symbol(stats::binomial(), n_trials = 5L),
    "binomial"
  )
  expect_equal(
    hsquared:::hs_nongaussian_family_symbol(stats::poisson()),
    "poisson"
  )
  # vector-safe: any per-record trial > 1 is Binomial, even if the first is 1
  expect_equal(
    hsquared:::hs_nongaussian_family_symbol(
      stats::binomial(),
      n_trials = c(1L, 4L, 5L)
    ),
    "binomial"
  )
  # an all-ones per-record vector reduces to Bernoulli
  expect_equal(
    hsquared:::hs_nongaussian_family_symbol(
      stats::binomial(),
      n_trials = c(1L, 1L, 1L)
    ),
    "bernoulli"
  )
})

test_that("the live bridge fits a balanced binomial-counts model [live]", {
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for the live bridge."
  )

  set.seed(5)
  ped <- data.frame(
    id = c("s1", "s2", "d1", "d2", paste0("a", 1:16)),
    sire = c(NA, NA, NA, NA, rep(c("s1", "s2"), 8)),
    dam = c(NA, NA, NA, NA, rep(c("d1", "d2"), 8))
  )
  n <- nrow(ped)
  trials <- 10L
  # A flat-probability rbinom() draw carries zero additive genetic signal by
  # construction: the fitted sigma_a2 then rides the search bracket's lower
  # rail instead of estimating anything (HSquared.jl#342 refuses that as a
  # boundary fit). Draw the success probability from a pedigree-based
  # breeding value with a real sigma_a2 = 1 instead (hsquared#225 -- the R
  # bridge cannot forward a non-default `initial`, so every live nongaussian
  # fit starts its search at sigma_a2 = 1).
  a <- hs_sim_genedrop_bv(ped, sigma_a2 = 1, seed = 5)
  p <- stats::plogis(stats::qlogis(0.4) + a)
  set.seed(5)
  succ <- rbinom(n, trials, p)
  dat <- data.frame(succ = succ, fail = trials - succ, id = ped$id)

  fit <- hsquared(
    cbind(succ, fail) ~ animal(1 | id, pedigree = ped),
    data = dat,
    family = stats::binomial(),
    REML = TRUE,
    control = hs_control(
      engine = "julia",
      engine_control = list(target = "nongaussian")
    )
  )
  expect_s3_class(fit, "hsquared_fit")
  expect_equal(fit$result$family, "binomial")
  vc <- variance_components(fit)
  expect_identical(vc$component, c("V_A", "V_RE", "V_O"))
  expect_true(all(is.finite(vc$estimate)))
  # V_A landing inside the search bracket's interior confirms the fit is not
  # riding the boundary refusal it would hit at sa0*exp(-+6) (HSquared.jl#342;
  # hsquared#225 -- the bridge always starts the search at sa0 = 1).
  sa2 <- vc$estimate[match("V_A", vc$component)]
  expect_gt(sa2, exp(-5.9))
  expect_lt(sa2, exp(5.9))
  expect_identical(vc$estimate[match("V_RE", vc$component)], 0)
  expect_identical(vc$estimate[match("V_O", vc$component)], 0)
  expect_equal(nrow(breeding_values(fit)), n)
  h2 <- heritability(fit)
  expect_identical(
    h2$field,
    c("h2_latent", "h2_liability", "h2_observation")
  )
  expect_true(all(is.finite(h2$estimate)))

  # parity: the R binomial-counts fit matches a direct engine fit_laplace_reml
  # with family = :binomial and the common n_trials (the bridge left hsq_*).
  direct_sa2 <- JuliaCall::julia_eval(
    "HSquared.fit_laplace_reml(hsq_y, hsq_X, hsq_Z, hsq_Ainv; family = :binomial, n_trials = Int(hsq_n_trials), ids = hsq_ped.ids).variance_components.sigma_a2"
  )
  expect_equal(
    vc$estimate[match("V_A", vc$component)],
    direct_sa2,
    tolerance = 1e-6
  )
})

test_that("a cbind binomial with one trial reduces to the Bernoulli fit [live]", {
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for the live bridge."
  )

  set.seed(106)
  ped <- data.frame(
    id = c("s1", "s2", "d1", "d2", paste0("a", 1:16)),
    sire = c(NA, NA, NA, NA, rep(c("s1", "s2"), 8)),
    dam = c(NA, NA, NA, NA, rep(c("d1", "d2"), 8))
  )
  n <- nrow(ped)
  # A flat-probability draw carries zero additive genetic signal by
  # construction, so the fitted sigma_a2 rides the search bracket's lower
  # rail instead of estimating anything (HSquared.jl#342 refuses that as a
  # boundary fit). Draw from a pedigree-based breeding value with a real
  # sigma_a2 = 2 instead (hsquared#225 -- the R bridge cannot forward a
  # non-default `initial`, so every live nongaussian fit starts its search
  # at sigma_a2 = 1). A single Bernoulli trial per animal is only weakly
  # informative about sigma_a2, so seed 106 was checked directly against
  # HSquared.jl::fit_laplace_reml to land comfortably interior for both the
  # cbind and binary routes below (sigma_a2 ~= 0.26, > 100x the lower rail).
  a <- hs_sim_genedrop_bv(ped, sigma_a2 = 2, seed = 106)
  p <- stats::plogis(a) # qlogis(0.5) == 0
  set.seed(106)
  y01 <- rbinom(n, 1L, p)
  ng_control <- hs_control(
    engine = "julia",
    engine_control = list(target = "nongaussian")
  )

  # cbind(successes, failures) with one trial each: n_trials = 1 -> Bernoulli
  fit_cbind <- hsquared(
    cbind(y, no) ~ animal(1 | id, pedigree = ped),
    data = data.frame(y = y01, no = 1L - y01, id = ped$id),
    family = stats::binomial(),
    REML = TRUE,
    control = ng_control
  )
  # the same data as a binary 0/1 Bernoulli response
  fit_binary <- hsquared(
    y ~ animal(1 | id, pedigree = ped),
    data = data.frame(y = y01, id = ped$id),
    family = stats::binomial(),
    REML = TRUE,
    control = ng_control
  )
  expect_equal(fit_cbind$result$family, "bernoulli")
  expect_equal(
    variance_components(fit_cbind)$estimate,
    variance_components(fit_binary)$estimate,
    tolerance = 1e-8
  )
  # V_A landing inside the search bracket's interior confirms the fit is not
  # riding the boundary refusal it would hit at sa0*exp(-+6) (HSquared.jl#342;
  # hsquared#225 -- the bridge always starts the search at sa0 = 1).
  vc_cbind <- variance_components(fit_cbind)
  sa2 <- vc_cbind$estimate[match("V_A", vc_cbind$component)]
  expect_gt(sa2, exp(-5.9))
  expect_lt(sa2, exp(5.9))
})

test_that("the live bridge fits a binomial-counts model with per-record varying trials [live]", {
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for the live bridge."
  )

  set.seed(7)
  ped <- data.frame(
    id = c("s1", "s2", "d1", "d2", paste0("a", 1:16)),
    sire = c(NA, NA, NA, NA, rep(c("s1", "s2"), 8)),
    dam = c(NA, NA, NA, NA, rep(c("d1", "d2"), 8))
  )
  n <- nrow(ped)
  trials <- sample(2:12, n, replace = TRUE) # per-record VARYING trials
  # A flat-probability draw carries zero additive genetic signal by
  # construction, so the fitted sigma_a2 rides the search bracket's lower
  # rail instead of estimating anything (HSquared.jl#342 refuses that as a
  # boundary fit). Draw from a pedigree-based breeding value with a real
  # sigma_a2 = 1 instead (hsquared#225 -- the R bridge cannot forward a
  # non-default `initial`, so every live nongaussian fit starts its search
  # at sigma_a2 = 1).
  a <- hs_sim_genedrop_bv(ped, sigma_a2 = 1, seed = 7)
  p <- stats::plogis(stats::qlogis(0.4) + a)
  set.seed(7)
  succ <- rbinom(n, trials, p)
  dat <- data.frame(succ = succ, fail = trials - succ, id = ped$id)

  fit <- hsquared(
    cbind(succ, fail) ~ animal(1 | id, pedigree = ped),
    data = dat,
    family = stats::binomial(),
    REML = TRUE,
    control = hs_control(
      engine = "julia",
      engine_control = list(target = "nongaussian")
    )
  )
  expect_s3_class(fit, "hsquared_fit")
  expect_equal(fit$result$family, "binomial")
  vc <- variance_components(fit)
  expect_identical(vc$component, c("V_A", "V_RE", "V_O"))
  expect_true(all(is.finite(vc$estimate)))
  # V_A landing inside the search bracket's interior confirms the fit is not
  # riding the boundary refusal it would hit at sa0*exp(-+6) (HSquared.jl#342;
  # hsquared#225 -- the bridge always starts the search at sa0 = 1).
  sa2 <- vc$estimate[match("V_A", vc$component)]
  expect_gt(sa2, exp(-5.9))
  expect_lt(sa2, exp(5.9))
  expect_identical(vc$estimate[match("V_RE", vc$component)], 0)
  expect_identical(vc$estimate[match("V_O", vc$component)], 0)
  expect_equal(nrow(breeding_values(fit)), n)
  h2 <- heritability(fit)
  expect_identical(
    h2$field,
    c("h2_latent", "h2_liability", "h2_observation")
  )
  expect_true(all(is.finite(h2$estimate[h2$field != "h2_observation"])))
  expect_true(is.nan(h2$estimate[h2$field == "h2_observation"]))
  expect_identical(
    h2$undefined_reason[h2$field == "h2_observation"],
    "varying_trials_no_scalar_estimand"
  )
  expect_identical(
    fit$result$h2_observation_undefined_reason,
    "varying_trials_no_scalar_estimand"
  )

  # parity: matches a direct engine fit with the per-record n_trials VECTOR
  # (the bridge left hsq_* in the Julia session; varying totals -> Vector{Int}).
  direct_sa2 <- JuliaCall::julia_eval(
    "HSquared.fit_laplace_reml(hsq_y, hsq_X, hsq_Z, hsq_Ainv; family = :binomial, n_trials = Vector{Int}(hsq_n_trials), ids = hsq_ped.ids).variance_components.sigma_a2"
  )
  expect_equal(
    vc$estimate[match("V_A", vc$component)],
    direct_sa2,
    tolerance = 1e-6
  )
})
