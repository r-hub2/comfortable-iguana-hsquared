# The 0.9 three-field non-Gaussian contract is deliberately tested without a
# Julia session.  It is a private wire-normalization and admission boundary;
# the legacy `nongaussian_result_payload` path remains a separate compatibility
# surface.

ng09_payload <- function(n = 4L, trials = NULL, predictors = FALSE) {
  X <- matrix(1, nrow = n, ncol = if (predictors) 2L else 1L)
  colnames(X) <- if (predictors) c("(Intercept)", "x") else "(Intercept)"
  structure(
    list(
      y = c(0, 1, 1, 0)[seq_len(n)],
      X = X,
      n_trials = trials,
      metadata = list(fixed_colnames = colnames(X))
    ),
    class = c("hs_bridge_payload", "list")
  )
}

ng09_raw <- function(family = "poisson", method = "laplace", n_trials = NULL) {
  components <- c(V_A = 0.4, V_RE = 0, V_O = 0)
  v_eta <- sum(components)
  mu <- 0.3
  varying_trials <- identical(family, "binomial") &&
    length(n_trials) > 1L && length(unique(n_trials)) > 1L
  list(
    schema = "nongaussian_three_field_v09",
    family = family,
    method = method,
    loglik = -12.5,
    converged = TRUE,
    boundary = FALSE,
    breeding_ids = c("a", "b"),
    breeding_values = c(0.1, -0.1),
    components = components,
    fixed_effects = c("(Intercept)" = mu),
    h2_latent = components[["V_A"]] / v_eta,
    h2_liability = if (identical(family, "poisson")) {
      NULL
    } else {
      components[["V_A"]] / (v_eta + pi^2 / 3)
    },
    h2_observation = if (identical(family, "poisson")) {
      components[["V_A"]] / (expm1(v_eta) + exp(-(mu + v_eta / 2)))
    } else if (isTRUE(varying_trials)) {
      NaN
    } else {
      0.125
    },
    h2_observation_undefined_reason = if (identical(family, "poisson")) {
      NULL
    } else if (isTRUE(varying_trials)) {
      "varying_trials_no_scalar_estimand"
    } else {
      NULL
    },
    n_trials = n_trials
  )
}

ng09_julia_raw <- function(family = "poisson", method = "laplace", n_trials = NULL) {
  raw <- ng09_raw(family = family, method = method, n_trials = n_trials)
  raw$components <- list(
    names = names(raw$components), values = unname(raw$components)
  )
  raw$fixed_effects <- list(
    names = names(raw$fixed_effects), values = unname(raw$fixed_effects)
  )
  raw
}

test_that("the v0.9 normalizer retains the three ratified Poisson fields", {
  raw <- ng09_raw()
  result <- hsquared:::hs_normalize_nongaussian_three_field_v09(
    raw,
    ng09_payload()
  )

  expect_equal(result$family, "poisson")
  expect_equal(result$marginal_method, "laplace")
  expect_equal(result$loglik_kind, "laplace marginal loglik")
  expect_identical(result$h2_latent, raw$h2_latent)
  expect_identical(result$h2_observation, raw$h2_observation)
  expect_identical(
    result$h2_observation_label,
    "count-scale observation h2 (conditional)"
  )
  expect_false("h2_liability" %in% names(result))
  expect_false("h2_observation_undefined_reason" %in% names(result))
  expect_false("n_trials" %in% names(result))
  expect_identical(result$breeding_values$id, c("a", "b"))
  expect_identical(result$breeding_values$value, c(0.1, -0.1))
  expect_identical(result$random_effects$animal, result$breeding_values)
})

test_that("the v0.9 normalizer consumes Julia names-values containers", {
  raw <- ng09_julia_raw()
  result <- hsquared:::hs_normalize_nongaussian_three_field_v09(
    raw,
    ng09_payload()
  )

  expect_identical(result$variance_components$component, c("V_A", "V_RE", "V_O"))
  expect_identical(result$fixed_effects, c("(Intercept)" = 0.3))

  reordered <- ng09_julia_raw()
  reordered$components <- reordered$components[c("values", "names")]
  reordered$fixed_effects <- reordered$fixed_effects[c("values", "names")]
  expect_silent(hsquared:::hs_normalize_nongaussian_three_field_v09(
    reordered,
    ng09_payload()
  ))
})

test_that("the v0.9 normalizer requires a truthful converged wire field", {
  raw <- ng09_julia_raw()
  result <- hsquared:::hs_normalize_nongaussian_three_field_v09(
    raw,
    ng09_payload()
  )
  expect_true(result$converged)

  bad <- raw
  bad$converged <- FALSE
  expect_error(
    hsquared:::hs_normalize_nongaussian_three_field_v09(bad, ng09_payload()),
    "converged"
  )
  bad <- raw
  bad$converged <- NULL
  expect_error(
    hsquared:::hs_normalize_nongaussian_three_field_v09(bad, ng09_payload()),
    "missing required `converged`"
  )
  bad <- raw
  bad$converged <- "yes"
  expect_error(
    hsquared:::hs_normalize_nongaussian_three_field_v09(bad, ng09_payload()),
    "converged"
  )
})

test_that("the v0.9 normalizer preserves scalar and varying Binomial trials", {
  scalar <- hsquared:::hs_normalize_nongaussian_three_field_v09(
    ng09_julia_raw("binomial", n_trials = 3L),
    ng09_payload()
  )
  expect_identical(scalar$n_trials, 3L)
  expect_identical(scalar$h2_observation, 0.125)
  expect_false("h2_observation_undefined_reason" %in% names(scalar))

  common_vector <- hsquared:::hs_normalize_nongaussian_three_field_v09(
    ng09_julia_raw("binomial", n_trials = rep(3L, 4L)),
    ng09_payload()
  )
  expect_identical(common_vector$n_trials, rep(3L, 4L))
  expect_identical(common_vector$h2_observation, 0.125)
  expect_false("h2_observation_undefined_reason" %in% names(common_vector))

  varying_trials <- c(2L, 3L, 4L, 5L)
  varying <- hsquared:::hs_normalize_nongaussian_three_field_v09(
    ng09_julia_raw("binomial", n_trials = varying_trials),
    ng09_payload()
  )
  expect_identical(varying$n_trials, varying_trials)
  expect_true(is.nan(varying$h2_observation))
  expect_identical(
    varying$h2_observation_undefined_reason,
    "varying_trials_no_scalar_estimand"
  )
})

test_that("the v0.9 normalizer exposes only ratified labelled scale rows", {
  poisson <- hsquared:::hs_normalize_nongaussian_three_field_v09(
    ng09_julia_raw(), ng09_payload()
  )
  expect_identical(
    poisson$heritability$label,
    c("latent-scale h2 (conditional)", "count-scale observation h2 (conditional)")
  )

  logit <- hsquared:::hs_normalize_nongaussian_three_field_v09(
    ng09_julia_raw("bernoulli"), ng09_payload()
  )
  expect_identical(
    logit$heritability$label,
    c(
      "latent-scale h2 (conditional)",
      "liability-scale h2 (conditional)",
      "observation-scale h2 (conditional)"
    )
  )
  expect_identical(logit$heritability$estimate[[3L]], 0.125)
  expect_true(is.na(logit$heritability$undefined_reason[[3L]]))
})

test_that("the existing non-Gaussian target validates before Julia setup", {
  payload <- ng09_payload()
  payload$pedigree <- list(id = character(), sire = character(), dam = character())
  expect_error(
    hsquared:::hs_fit_julia_nongaussian_payload(
      payload,
      project = tempfile(),
      family = stats::Gamma()
    ),
    "admits only poisson\\(log\\) and binomial\\(logit\\)"
  )
})

test_that("the public non-Gaussian target rejects weights and dots before Julia", {
  ped <- data.frame(
    id = c("s", "d", "a", "b"),
    sire = c(NA, NA, "s", "s"),
    dam = c(NA, NA, "d", "d")
  )
  dat <- data.frame(y = c(0, 1, 2, 3), id = ped$id)
  control <- hs_control(
    engine = "julia",
    engine_control = list(target = "nongaussian", julia_project = tempfile())
  )
  expect_error(
    hsquared(
      y ~ animal(1 | id, pedigree = ped),
      data = dat,
      family = stats::poisson(),
      control = control,
      weights = rep(1, nrow(dat))
    ),
    "weights"
  )
  expect_error(
    hsquared(
      y ~ animal(1 | id, pedigree = ped),
      data = dat,
      family = stats::poisson(),
      control = control,
      offset = 0
    ),
    "does not accept"
  )
})

test_that("the v0.9 Julia command builds a complete explicit envelope", {
  command <- hsquared:::hs_nongaussian_three_field_julia_command(
    family_symbol = "binomial",
    marginal = "variational",
    n_trials = c(2L, 3L, 4L, 5L)
  )

  expect_match(command, "HSquared\\.nongaussian_three_field_payload")
  expect_match(command, "predictor_variance = 0\\.0")
  expect_match(command, "response_length = length\\(hsq_y\\)")
  expect_match(command, "n_trials = Vector\\{Int\\}\\(hsq_n_trials\\)")
  for (key in c(
    "schema", "family", "method", "loglik", "components", "fixed_effects",
    "h2_latent", "h2_liability", "h2_observation",
    "h2_observation_undefined_reason", "n_trials", "breeding_ids",
    "breeding_values", "converged"
  )) {
    expect_match(command, paste0("\\\"", key, "\\\" =>"), fixed = FALSE)
  }
  expect_match(command, "\\\"converged\\\" => hsq_fit\\.converged")
})

test_that("the v0.9 normalizer distinguishes defined and non-scalar logit observations", {
  raw <- ng09_raw("bernoulli", method = "VA")
  result <- hsquared:::hs_normalize_nongaussian_three_field_v09(
    raw,
    ng09_payload()
  )

  expect_equal(result$marginal_method, "variational")
  expect_equal(result$loglik_kind, "elbo (variational lower bound)")
  expect_identical(result$h2_liability, raw$h2_liability)
  expect_identical(result$h2_observation, 0.125)
  expect_identical(result$h2_observation_label, "observation-scale h2 (conditional)")
  expect_false("h2_observation_undefined_reason" %in% names(result))

  bad <- raw
  bad$h2_observation <- NA_real_
  expect_error(
    hsquared:::hs_normalize_nongaussian_three_field_v09(bad, ng09_payload()),
    "finite numeric scalar"
  )
  bad <- raw
  bad$h2_observation_undefined_reason <- ""
  expect_error(
    hsquared:::hs_normalize_nongaussian_three_field_v09(bad, ng09_payload()),
    "present `nothing`"
  )
  bad <- raw
  bad$h2_observation <- 1.01
  expect_error(
    hsquared:::hs_normalize_nongaussian_three_field_v09(bad, ng09_payload()),
    "in \\[0, 1\\]"
  )

  varying_raw <- ng09_julia_raw("binomial", n_trials = c(2L, 3L, 4L, 5L))
  varying_raw$h2_observation <- 0.125
  expect_error(
    hsquared:::hs_normalize_nongaussian_three_field_v09(
      varying_raw,
      ng09_payload()
    ),
    "literal NaN"
  )
})

test_that("the v0.9 normalizer rejects schema and estimand mutations exactly", {
  raw <- ng09_raw()
  bad <- raw
  bad$schema <- "nongaussian_result_payload"
  expect_error(
    hsquared:::hs_normalize_nongaussian_three_field_v09(bad, ng09_payload()),
    "nongaussian_three_field_v09"
  )

  bad <- raw
  bad$h2_latent <- bad$h2_latent + 1e-15
  expect_error(
    hsquared:::hs_normalize_nongaussian_three_field_v09(bad, ng09_payload()),
    "h2_latent"
  )
  bad <- raw
  bad$h2_observation <- raw$components[["V_A"]] /
    (expm1(raw$components[["V_A"]]) +
      exp(-raw$fixed_effects[["(Intercept)"]]))
  expect_error(
    hsquared:::hs_normalize_nongaussian_three_field_v09(bad, ng09_payload()),
    "h2_observation"
  )
  bad <- raw
  bad$h2_liability <- 0
  expect_error(
    hsquared:::hs_normalize_nongaussian_three_field_v09(bad, ng09_payload()),
    "Poisson"
  )
})

test_that("the v0.9 normalizer rejects non-structural V_RE and V_O", {
  raw <- ng09_raw()
  raw$components[["V_RE"]] <- 0.01
  raw$h2_latent <- raw$components[["V_A"]] / sum(raw$components)
  raw$h2_observation <- raw$components[["V_A"]] / (
    expm1(sum(raw$components)) +
      exp(-(raw$fixed_effects[["(Intercept)"]] + sum(raw$components) / 2))
  )
  expect_error(
    hsquared:::hs_normalize_nongaussian_three_field_v09(raw, ng09_payload()),
    "structural zero"
  )

  raw <- ng09_raw()
  raw$components[["V_O"]] <- 0.01
  raw$h2_latent <- raw$components[["V_A"]] / sum(raw$components)
  raw$h2_observation <- raw$components[["V_A"]] / (
    expm1(sum(raw$components)) +
      exp(-(raw$fixed_effects[["(Intercept)"]] + sum(raw$components) / 2))
  )
  expect_error(
    hsquared:::hs_normalize_nongaussian_three_field_v09(raw, ng09_payload()),
    "structural zero"
  )
})

test_that("the v0.9 admission helper classifies all four ratified input cells", {
  poisson <- ng09_payload()
  poisson$y <- c(0, 1, 2, 3)
  expect_identical(
    hsquared:::hs_validate_nongaussian_three_field_v09_admission(
      poisson,
      stats::poisson(),
      marginal = "LA"
    ),
    list(family = "poisson", method = "laplace", n_trials = NULL)
  )

  bernoulli <- ng09_payload()
  expect_identical(
    hsquared:::hs_validate_nongaussian_three_field_v09_admission(
      bernoulli,
      stats::binomial(),
      marginal = "Variational"
    ),
    list(family = "bernoulli", method = "variational", n_trials = NULL)
  )

  common <- ng09_payload(trials = rep(3L, 4L))
  common$y <- c(0, 1, 2, 3)
  expect_identical(
    hsquared:::hs_validate_nongaussian_three_field_v09_admission(
      common,
      stats::binomial()
    ),
    list(family = "binomial", method = "laplace", n_trials = 3L)
  )

  varying <- ng09_payload(trials = c(2L, 3L, 4L, 5L))
  varying$y <- c(0, 1, 2, 3)
  expect_identical(
    hsquared:::hs_validate_nongaussian_three_field_v09_admission(
      varying,
      stats::binomial()
    ),
    list(
      family = "binomial", method = "laplace",
      n_trials = c(2L, 3L, 4L, 5L)
    )
  )
})

test_that("the v0.9 admission helper rejects predictor, trial, weight, and dot shortcuts", {
  expect_error(
    hsquared:::hs_validate_nongaussian_three_field_v09_admission(
      ng09_payload(predictors = TRUE), stats::poisson()
    ),
    "one intercept"
  )
  expect_error(
    hsquared:::hs_validate_nongaussian_three_field_v09_admission(
      ng09_payload(), stats::poisson(), predictor_variance = 0.01
    ),
    "predictor_variance"
  )
  bad_trials <- ng09_payload(trials = c(3L, 3L))
  bad_trials$y <- c(0, 1, 2, 3)
  expect_error(
    hsquared:::hs_validate_nongaussian_three_field_v09_admission(
      bad_trials, stats::binomial()
    ),
    "response length"
  )
  expect_error(
    hsquared:::hs_validate_nongaussian_three_field_v09_admission(
      ng09_payload(), stats::binomial(), weights = rep(1, 4)
    ),
    "weights"
  )
  expect_error(
    hsquared:::hs_validate_nongaussian_three_field_v09_admission(
      ng09_payload(), stats::binomial(), dots = list(offset = 1)
    ),
    "does not accept"
  )
})

test_that("the legacy non-Gaussian normalizer remains a separate compatibility path", {
  expect_true(exists("hs_normalize_nongaussian_result", envir = asNamespace("hsquared")))
  expect_false(identical(
    hsquared:::hs_normalize_nongaussian_result,
    hsquared:::hs_normalize_nongaussian_three_field_v09
  ))
})

test_that("hsquared#222: `hs_ng09_boundary` passes FALSE and refuses TRUE with a classed error", {
  raw <- ng09_raw()
  expect_false(hsquared:::hs_ng09_boundary(raw))

  boundary_raw <- raw
  boundary_raw$boundary <- TRUE
  err <- tryCatch(
    hsquared:::hs_ng09_boundary(boundary_raw),
    error = function(e) e
  )
  expect_s3_class(err, "hsquared_error")
  expect_match(conditionMessage(err), "boundary", ignore.case = TRUE)
  expect_match(conditionMessage(err), "initial")
  expect_match(conditionMessage(err), "restart_check")

  bad <- raw
  bad$boundary <- NA
  expect_error(hsquared:::hs_ng09_boundary(bad), "boundary")
  bad <- raw
  bad$boundary <- "true"
  expect_error(hsquared:::hs_ng09_boundary(bad), "boundary")
  bad <- raw
  bad$boundary <- NULL
  expect_error(hsquared:::hs_ng09_boundary(bad), "missing required `boundary`")
})

test_that("hsquared#222: the v0.9 normalizer surfaces `boundary` next to `converged` and refuses boundary = TRUE", {
  raw <- ng09_raw()
  result <- hsquared:::hs_normalize_nongaussian_three_field_v09(raw, ng09_payload())
  expect_identical(result$boundary, FALSE)

  bad <- raw
  bad$boundary <- TRUE
  err <- tryCatch(
    hsquared:::hs_normalize_nongaussian_three_field_v09(bad, ng09_payload()),
    error = function(e) e
  )
  expect_s3_class(err, "hsquared_error")
  expect_match(conditionMessage(err), "boundary", ignore.case = TRUE)
})

test_that("the v0.9 Poisson route carries the three-field result through the live bridge", {
  # This is deliberately a tiny deterministic integration check, not a
  # calibration run.  Resolve the explicitly configured Julia project so Tier-1
  # CI and local callers exercise the intended checkout rather than a retired
  # temporary worktree.
  project <- hsquared:::hs_default_julia_project()
  hs_require_bridge("A3 v0.9 non-Gaussian bridge", project = project)

  # Two three-progeny sire families with clearly separated Poisson means
  # (~6.5 vs ~0.5) so sigma_a2 is genuinely identifiable. The original n = 4,
  # one-full-sib-pair fixture carried no information for sigma_a2 separate
  # from Poisson noise, so the fitted value rode the search bracket's lower
  # rail (HSquared.jl#342 refuses that as a boundary fit); the R bridge
  # cannot forward a non-default `initial` (hsquared#225), so every live
  # nongaussian fit starts its search at sigma_a2 = 1. This stays a tiny
  # deterministic fixture, not a calibration run.
  pedigree <- data.frame(
    id = c("s1", "d1", "s2", "d2", "a1", "a2", "a3", "b1", "b2", "b3"),
    sire = c(NA, NA, NA, NA, "s1", "s1", "s1", "s2", "s2", "s2"),
    dam = c(NA, NA, NA, NA, "d1", "d1", "d1", "d2", "d2", "d2")
  )
  data <- data.frame(
    y = c(6, 7, 0, 1, 6, 7, 5, 0, 1, 0),
    id = pedigree$id
  )
  fit <- hsquared(
    y ~ animal(1 | id, pedigree = pedigree),
    data = data,
    family = stats::poisson(),
    control = hs_control(
      engine = "julia",
      engine_control = list(
        target = "nongaussian",
        julia_project = project,
        iterations = 20L
      )
    )
  )

  expect_s3_class(fit, "hsquared_fit")
  expect_identical(fit$engine, "HSquared.jl")
  expect_identical(fit$spec$target, "nongaussian")
  expect_identical(fit$result$family, "poisson")
  expect_true(isTRUE(fit$result$converged))
  expect_identical(
    fit$result$heritability$field,
    c("h2_latent", "h2_observation")
  )
  expect_identical(
    fit$result$heritability$label,
    c(
      "latent-scale h2 (conditional)",
      "count-scale observation h2 (conditional)"
    )
  )
  expect_true(all(is.finite(fit$result$heritability$estimate)))
  # V_A landing inside the search bracket's interior confirms the fit is not
  # riding the boundary refusal it would hit at sa0*exp(-+6) (HSquared.jl#342;
  # hsquared#225 -- the bridge always starts the search at sa0 = 1).
  vc <- variance_components(fit)
  sa2 <- vc$estimate[match("V_A", vc$component)]
  expect_gt(sa2, exp(-5.9))
  expect_lt(sa2, exp(5.9))
  expect_false("h2_liability" %in% names(fit$result))
  expect_false("n_trials" %in% names(fit$result))
  expect_identical(
    heritability(fit)$field,
    fit$result$heritability$field
  )

  # The live route selects only the versioned normalizer; the legacy public
  # compatibility normalizer remains a distinct function.
  expect_false(identical(
    hsquared:::hs_normalize_nongaussian_result,
    hsquared:::hs_normalize_nongaussian_three_field_v09
  ))
})

test_that("the v0.9 Binomial route carries A4-1 observation semantics through the live bridge", {
  # This tiny configured bridge check proves the Julia producer and R consumer
  # agree on the A4-1 field distinction; it is not a calibration campaign.
  project <- hsquared:::hs_default_julia_project()
  hs_require_bridge("A4-1 v0.9 Binomial bridge", project = project)

  pedigree <- data.frame(
    id = c("s1", "d1", "s2", "d2", "a1", "a2", "b1", "b2"),
    sire = c(NA, NA, NA, NA, "s1", "s1", "s2", "s2"),
    dam = c(NA, NA, NA, NA, "d1", "d1", "d2", "d2")
  )
  control <- hs_control(
    engine = "julia",
    engine_control = list(
      target = "nongaussian",
      julia_project = project,
      iterations = 40L
    )
  )

  common <- data.frame(
    successes = c(0, 1, 2, 3, 0, 1, 2, 3),
    failures = c(3, 2, 1, 0, 3, 2, 1, 0),
    id = pedigree$id
  )
  common_fit <- hsquared(
    cbind(successes, failures) ~ animal(1 | id, pedigree = pedigree),
    data = common,
    family = stats::binomial(),
    control = control
  )
  expect_identical(common_fit$result$family, "binomial")
  expect_true(is.finite(common_fit$result$h2_observation))
  expect_true(common_fit$result$h2_observation >= 0)
  expect_true(common_fit$result$h2_observation <= 1)
  expect_false("h2_observation_undefined_reason" %in% names(common_fit$result))

  # The original perfect 0/1 alternation had exactly one Bernoulli trial per
  # animal-random-effect level (8 records, 8 animals): that is genuine
  # quasi-complete separation -- every 0/1 pattern at that n_animals ==
  # n_records ratio drives sigma_a2 -> the search bracket's upper rail
  # (HSquared.jl#342 refuses that as a boundary fit), independent of the
  # pattern chosen. Give each animal 3 replicate Bernoulli records instead
  # (still all-ones trials, so the family classification below stays
  # "bernoulli"), with mixed (not uniformly identical) outcomes within both
  # full-sib families, so the fit is no longer separable.
  all_one_succ <- c(
    1, 1, 0, # s1
    1, 0, 1, # d1
    0, 0, 1, # s2
    0, 1, 0, # d2
    1, 1, 1, # a1
    0, 1, 1, # a2
    0, 0, 1, # b1
    1, 0, 0 # b2
  )
  all_one <- data.frame(
    successes = all_one_succ,
    failures = 1L - all_one_succ,
    id = rep(pedigree$id, each = 3)
  )
  all_one_fit <- hsquared(
    cbind(successes, failures) ~ animal(1 | id, pedigree = pedigree),
    data = all_one,
    family = stats::binomial(),
    control = control
  )
  expect_identical(all_one_fit$result$family, "bernoulli")
  expect_true(is.finite(all_one_fit$result$h2_observation))
  expect_false("h2_observation_undefined_reason" %in% names(all_one_fit$result))
  # V_A landing inside the search bracket's interior confirms the fit is not
  # riding the boundary refusal it would hit at sa0*exp(-+6) (HSquared.jl#342;
  # hsquared#225 -- the bridge always starts the search at sa0 = 1).
  all_one_vc <- variance_components(all_one_fit)
  all_one_sa2 <- all_one_vc$estimate[match("V_A", all_one_vc$component)]
  expect_gt(all_one_sa2, exp(-5.9))
  expect_lt(all_one_sa2, exp(5.9))

  varying <- data.frame(
    successes = c(0, 1, 3, 4, 1, 2, 1, 5),
    failures = c(2, 2, 1, 1, 1, 1, 3, 0),
    id = pedigree$id
  )
  varying_fit <- hsquared(
    cbind(successes, failures) ~ animal(1 | id, pedigree = pedigree),
    data = varying,
    family = stats::binomial(),
    control = control
  )
  expect_identical(varying_fit$result$family, "binomial")
  expect_true(is.nan(varying_fit$result$h2_observation))
  expect_identical(
    varying_fit$result$h2_observation_undefined_reason,
    "varying_trials_no_scalar_estimand"
  )
})

# hsquared#225: engine_control$initial and $restart_check on target =
# "nongaussian" -- live evidence that the bridge actually forwards them.

ng225_pedigree <- function() {
  data.frame(
    id = c("s1", "s2", "d1", "d2", paste0("a", 1:16)),
    sire = c(NA, NA, NA, NA, rep(c("s1", "s2"), 8)),
    dam = c(NA, NA, NA, NA, rep(c("d1", "d2"), 8))
  )
}

test_that("hsquared#225: supplying `initial` on a healthy fixture fits and agrees with the default", {
  project <- hsquared:::hs_default_julia_project()
  hs_require_bridge("hsquared#225 initial forwarding", project = project)

  ped <- ng225_pedigree()
  n <- nrow(ped)
  # A real pedigree-based breeding value (not a flat-probability draw) so
  # sigma_a2 is genuinely identifiable and lands well inside the search
  # bracket's interior at the default start -- the same fixture shape as the
  # binomial-counts live tests above (hsquared#227).
  a <- hs_sim_genedrop_bv(ped, sigma_a2 = 2, seed = 106)
  p <- stats::plogis(a)
  set.seed(106)
  y01 <- rbinom(n, 1L, p)
  dat <- data.frame(y = y01, id = ped$id)

  fit_default <- hsquared(
    y ~ animal(1 | id, pedigree = ped),
    data = dat,
    family = stats::binomial(),
    control = hs_control(
      engine = "julia",
      engine_control = list(target = "nongaussian", julia_project = project)
    )
  )
  fit_initial <- hsquared(
    y ~ animal(1 | id, pedigree = ped),
    data = dat,
    family = stats::binomial(),
    control = hs_control(
      engine = "julia",
      engine_control = list(
        target = "nongaussian",
        julia_project = project,
        initial = list(sigma_a2 = 0.05)
      )
    )
  )

  sa2_default <- variance_components(fit_default)$estimate[
    match("V_A", variance_components(fit_default)$component)
  ]
  sa2_initial <- variance_components(fit_initial)$estimate[
    match("V_A", variance_components(fit_initial)$component)
  ]
  expect_true(is.finite(sa2_default))
  expect_true(is.finite(sa2_initial))
  expect_false(isTRUE(fit_default$result$boundary))
  expect_false(isTRUE(fit_initial$result$boundary))
  # Agreement to 1e-4 shows a supplied `initial` does not perturb a
  # well-identified fit. It is NOT evidence that `initial` reached the Julia
  # call -- this assertion also passes if the forwarding is reverted. The
  # forwarding itself is pinned by the command-string test in
  # test-nongaussian-boundary-and-initial.R.
  expect_equal(sa2_default, sa2_initial, tolerance = 1e-4)
})

# hsquared#230: `fit_diagnostics()` must surface the bridge's `boundary` flag
# as `search_boundary` for a real non-Gaussian fit, not only via
# `fit$result$boundary` (the pre-#230 access path documented in ?hs_control).
test_that("hsquared#230: fit_diagnostics() reports search_boundary = FALSE for a healthy live non-Gaussian fit", {
  project <- hsquared:::hs_default_julia_project()
  hs_require_bridge("hsquared#230 search_boundary row", project = project)

  ped <- ng225_pedigree()
  n <- nrow(ped)
  a <- hs_sim_genedrop_bv(ped, sigma_a2 = 2, seed = 106)
  p <- stats::plogis(a)
  set.seed(106)
  y01 <- rbinom(n, 1L, p)
  dat <- data.frame(y = y01, id = ped$id)

  fit <- hsquared(
    y ~ animal(1 | id, pedigree = ped),
    data = dat,
    family = stats::binomial(),
    control = hs_control(
      engine = "julia",
      engine_control = list(target = "nongaussian", julia_project = project)
    )
  )

  diag <- fit_diagnostics(fit)
  expect_false(isTRUE(fit$result$boundary))
  expect_equal(diag$value[diag$metric == "search_boundary"], "FALSE")
  condition_value <- diag$value[diag$metric == "search_boundary_condition"]
  expect_equal(length(condition_value), 1L)
  expect_equal(condition_value, "interior")
})

test_that("hsquared#222/#225: a boundary-riding fixture fails the DEFAULT fit with a classed, actionable error", {
  project <- hsquared:::hs_default_julia_project()
  hs_require_bridge("hsquared#222/#225 boundary refusal", project = project)

  ped <- ng225_pedigree()
  n <- nrow(ped)
  # A flat-probability Bernoulli draw carries zero additive genetic signal by
  # construction -- exactly the pattern hsquared#227 diagnosed and moved the
  # OTHER live fixtures away from. Kept here deliberately: this is the
  # boundary-riding case the DEFAULT (initial = NULL, restart_check = FALSE)
  # fit must still fail on, so #222's translated-error contract and #225's
  # advice ("retry with restart_check = TRUE or a different initial") are
  # exercised end to end.
  set.seed(42)
  y_flat <- rbinom(n, 1L, 0.4)
  dat_flat <- data.frame(y = y_flat, id = ped$id)
  control <- hs_control(
    engine = "julia",
    engine_control = list(target = "nongaussian", julia_project = project)
  )

  err <- tryCatch(
    hsquared(
      y ~ animal(1 | id, pedigree = ped),
      data = dat_flat,
      family = stats::binomial(),
      control = control
    ),
    error = function(e) e
  )
  expect_s3_class(err, "hsquared_julia_error")
  expect_s3_class(err, "hsquared_error")
  expect_match(conditionMessage(err), "boundary", ignore.case = TRUE)
  expect_match(conditionMessage(err), "initial")
  expect_match(conditionMessage(err), "restart_check")

  # restart_check = TRUE on the SAME healthy fixture from the test above
  # returns normally -- the lever works without requiring boundary-riding data.
  a <- hs_sim_genedrop_bv(ped, sigma_a2 = 2, seed = 106)
  p <- stats::plogis(a)
  set.seed(106)
  y01 <- rbinom(n, 1L, p)
  dat_healthy <- data.frame(y = y01, id = ped$id)
  fit_restart <- hsquared(
    y ~ animal(1 | id, pedigree = ped),
    data = dat_healthy,
    family = stats::binomial(),
    control = hs_control(
      engine = "julia",
      engine_control = list(
        target = "nongaussian",
        julia_project = project,
        restart_check = TRUE
      )
    )
  )
  expect_s3_class(fit_restart, "hsquared_fit")
  expect_false(isTRUE(fit_restart$result$boundary))
})
