# v0.7 genomic GREML activation: a single genomic random effect from either a
# supplied precision or the frozen sample-p VanRaden-1 marker construction.
# Gaussian REML remains an explicit Julia target after the activation pilot
# stopped; the default fit path must reject it.

hs_test_ginv <- function(ids) {
  n <- length(ids)
  g <- diag(n)
  for (i in seq_len(n - 1L)) {
    g[i, i + 1L] <- g[i + 1L, i] <- 0.2
  }
  ginv <- solve(g)
  dimnames(ginv) <- list(ids, ids)
  ginv
}

test_that("the parser accepts genomic(1 | id, Ginv = Ginv) as a primary effect", {
  ids <- paste0("g", 1:4)
  Ginv <- hs_test_ginv(ids)
  dat <- data.frame(y = c(1, 2, 3, 4), id = ids)

  spec <- hsquared:::hs_build_model_spec(
    y ~ genomic(1 | id, Ginv = Ginv),
    data = dat,
    family = stats::gaussian(),
    REML = TRUE
  )

  expect_named(spec$random, "genomic")
  expect_equal(spec$random$genomic$type, "genomic")
  expect_equal(spec$random$genomic$group, "id")
  expect_equal(spec$random$genomic$relationship, "genomic")
  expect_equal(spec$random$genomic$ids, ids)
  expect_match(spec$bridge$target, "Ginv", fixed = TRUE)
})

test_that("a formula must contain exactly one primary effect", {
  ids <- paste0("g", 1:2)
  Ginv <- hs_test_ginv(ids)
  ped <- data.frame(id = ids, sire = c(NA, NA), dam = c(NA, NA))
  dat <- data.frame(y = c(1, 2), id = ids)
  # animal() AND genomic() together is rejected (one primary effect only)
  expect_error(
    hsquared:::hs_build_model_spec(
      y ~ animal(1 | id, pedigree = ped) + genomic(1 | id, Ginv = Ginv),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "one primary",
    fixed = TRUE
  )
})

test_that("genomic() requires a Ginv argument", {
  ids <- paste0("g", 1:2)
  dat <- data.frame(y = c(1, 2), id = ids)
  expect_error(
    hsquared:::hs_build_model_spec(
      y ~ genomic(1 | id),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "requires a `Ginv`",
    fixed = TRUE
  )
})

test_that("genomic() ids must be in the Ginv dimnames", {
  ids <- paste0("g", 1:3)
  Ginv <- hs_test_ginv(ids)
  dat <- data.frame(y = c(1, 2), id = c("g1", "ghost"))
  expect_error(
    hsquared:::hs_build_model_spec(
      y ~ genomic(1 | id, Ginv = Ginv),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "not in the `Ginv`",
    fixed = TRUE
  )
})

test_that("genomic is a valid opt-in julia target", {
  expect_equal(hsquared:::hs_validate_julia_target("genomic"), "genomic")
})

test_that("the default fit path keeps genomic() opt-in", {
  ids <- paste0("g", 1:3)
  Ginv <- hs_test_ginv(ids)
  dat <- data.frame(y = c(1, 2, 3), id = ids)

  spec <- suppressMessages(hsquared(
    y ~ genomic(1 | id, Ginv = Ginv),
    data = dat,
    family = stats::gaussian(),
    control = hs_control(engine = "validate")
  ))
  expect_match(spec$bridge$target, "Ginv", fixed = TRUE)

  err <- tryCatch(
    hsquared(
      y ~ genomic(1 | id, Ginv = Ginv),
      data = dat,
      family = stats::gaussian()
    ),
    error = function(e) conditionMessage(e)
  )
  expect_true(grepl(
    "`genomic()` is not on the default",
    err,
    fixed = TRUE
  ))
  expect_true(grepl("Ordinary-call / default-route", err, fixed = TRUE))
  expect_true(grepl("BOUNDARY_HOLDOUT_FAIL", err, fixed = TRUE))
  expect_true(grepl("Closest working call", err, fixed = TRUE))
  expect_true(grepl("target = \"genomic\"", err, fixed = TRUE))
  expect_true(grepl("genomic_variance_ratio", err, fixed = TRUE))
  expect_true(grepl("public_covered_count stays 7", err, fixed = TRUE))
  expect_true(grepl("y ~ genomic(1 | id, Ginv = Ginv)", err, fixed = TRUE))
  expect_true(grepl("data = dat", err, fixed = TRUE))
  expect_true(grepl(
    "not a silent ordinary-route promotion",
    err,
    fixed = TRUE
  ))

  # Marker form gets the same ordinary-call hold (no silent default route).
  markers <- matrix(
    c(0, 1, 2, 1, 0, 2, 2, 1, 0),
    nrow = 3L,
    dimnames = list(ids, paste0("m", 1:3))
  )
  err_m <- tryCatch(
    hsquared(
      y ~ genomic(1 | id, markers = markers),
      data = dat,
      family = stats::gaussian()
    ),
    error = function(e) conditionMessage(e)
  )
  expect_true(grepl(
    "`genomic()` is not on the default",
    err_m,
    fixed = TRUE
  ))
  expect_true(grepl("target = \"genomic\"", err_m, fixed = TRUE))
  expect_true(grepl("snp_blup", err_m, fixed = TRUE))

  # Targetless engine = "julia" also refuses silent genomic routing.
  err_j <- tryCatch(
    hsquared(
      y ~ genomic(1 | id, Ginv = Ginv),
      data = dat,
      family = stats::gaussian(),
      control = hs_control(engine = "julia")
    ),
    error = function(e) conditionMessage(e)
  )
  expect_true(grepl("target = \"genomic\"", err_j, fixed = TRUE))
  expect_true(
    grepl("fit_animal_model", err_j, fixed = TRUE) ||
      grepl("single additive-genetic", err_j, fixed = TRUE)
  )

  # REML = FALSE stays rejected on the explicit genomic route (design-44).
  expect_error(
    hsquared(
      y ~ genomic(1 | id, Ginv = Ginv),
      data = dat,
      family = stats::gaussian(),
      REML = FALSE
    ),
    "REML = FALSE",
    fixed = TRUE
  )
})

test_that("formula_status names the ordinary-call genomic hold", {
  status <- formula_status()
  ginv_row <- status[
    status$term == "genomic(1 | id, Ginv = Ginv)",
    ,
    drop = FALSE
  ]
  expect_length(ginv_row$current_behavior, 1L)
  expect_match(
    ginv_row$current_behavior,
    "NOT on the default",
    fixed = TRUE
  )
  expect_match(
    ginv_row$current_behavior,
    "BOUNDARY_HOLDOUT_FAIL",
    fixed = TRUE
  )
  expect_match(
    ginv_row$current_behavior,
    "public_covered_count remains 7",
    fixed = TRUE
  )

  markers_row <- status[
    status$term == "genomic(1 | id, markers = M)",
    ,
    drop = FALSE
  ]
  expect_length(markers_row$current_behavior, 1L)
  expect_match(
    markers_row$current_behavior,
    "ordinary-call genomic activation remains held",
    fixed = TRUE
  )
})

test_that("the genomic bridge requires an internal payload", {
  expect_error(
    hsquared:::hs_fit_julia_genomic_payload(list()),
    "`payload` must be an internal `hs_bridge_payload`.",
    fixed = TRUE
  )
})

test_that("the explicit supplied-Ginv route fits [live]", {
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for live GREML."
  )

  set.seed(5)
  na <- 8
  ids <- paste0("g", seq_len(na))
  m <- matrix(stats::rbinom(na * 60, 2, 0.3), na, 60)
  mc <- scale(m, scale = FALSE)
  g <- tcrossprod(mc)
  g <- g / mean(diag(g)) + diag(na) * 0.01
  Ginv <- solve(g)
  dimnames(Ginv) <- list(ids, ids)

  n <- 24
  rec <- rep(ids, length.out = n)
  dat <- data.frame(
    y = 3 + stats::rnorm(n, 0, 1),
    id = rec
  )

  fit <- hsquared(
    y ~ genomic(1 | id, Ginv = Ginv),
    data = dat,
    family = stats::gaussian(),
    control = hs_control(
      engine = "julia",
      engine_control = list(target = "genomic")
    )
  )

  expect_s3_class(fit, "hsquared_fit")
  expect_null(fit$result$genomic_boundary)
  expect_equal(fit$spec$target, "genomic")
  vc <- variance_components(fit)
  expect_equal(vc$component, c("genomic", "residual"))
  expect_true(all(is.finite(vc$estimate)) && all(vc$estimate > 0))
  ratio <- heritability(fit)
  expect_equal(ratio$term, "genomic")
  expect_equal(ratio$component, "genomic_variance_ratio")
  expect_equal(ratio$relationship_source, "supplied_Ginv")
  expect_equal(ratio$relationship_scale, "inverse_of_supplied_precision")
  expect_true(is.na(ratio$relationship_method))
  expect_true(is.na(ratio$allele_frequency_source))
  expect_true(is.na(ratio$ridge))
  h2 <- ratio$estimate
  expect_true(is.finite(h2) && h2 > 0 && h2 < 1)
  expect_equal(nrow(breeding_values(fit)), na)
  expect_equal(breeding_values(fit)$id, ids)
  expect_match(
    fit$result$relationship_provenance$precision_fingerprint,
    "^[0-9a-f]{64}$"
  )
  expect_equal(
    fit_diagnostics(fit)$value[
      fit_diagnostics(fit)$metric == "variance_components_source"
    ],
    "estimated_genomic_ai_reml"
  )
})

test_that("one-record genomic bridge surfaces a scientific lower endpoint [live]", {
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for live boundary GREML."
  )
  set.seed(99)
  n <- 8L
  y <- stats::rnorm(n)
  A <- matrix(stats::rnorm(n * n), n, n)
  K <- crossprod(A) / n + diag(n) * 0.05
  ids <- paste0("b", seq_len(n))
  Q <- solve(K)
  dimnames(Q) <- list(ids, ids)
  dat <- data.frame(y = y, id = ids)
  fit <- hsquared(
    y ~ genomic(1 | id, Ginv = Q),
    data = dat,
    family = stats::gaussian(),
    control = hs_control(
      engine = "julia",
      engine_control = list(target = "genomic")
    )
  )
  expect_true(fit$result$converged)
  expect_identical(fit$result$genomic_boundary$status, "boundary_lower")
  expect_identical(heritability(fit)$estimate, 0)
  expect_identical(heritability(fit)$numerical_estimate, 1e-7)
  expect_null(fit$result[["breeding_values"]])
  expect_null(fit$result$prediction_error_variance)
  expect_null(fit$result$reliability)
  expect_error(breeding_values(fit), "does not contain")
  printed <- capture.output(print(fit))
  expect_true(any(grepl("scientific endpoint.*0", printed)))
  expect_true(any(grepl("numerical MME.*1e-07", printed)))
})

test_that("genomic() accepts a marker matrix to build the relationship", {
  ids <- paste0("g", 1:5)
  set.seed(4)
  M <- matrix(stats::rbinom(5 * 20, 2, 0.3), 5, 20)
  rownames(M) <- ids
  dat <- data.frame(y = c(1, 2, 3, 4, 5), id = ids)

  spec <- hsquared:::hs_build_model_spec(
    y ~ genomic(1 | id, markers = M),
    data = dat,
    family = stats::gaussian(),
    REML = TRUE
  )

  expect_equal(spec$random$genomic$type, "genomic")
  expect_equal(spec$random$genomic$source, "markers")
  expect_equal(spec$random$genomic$ids, ids)
  expect_true(is.matrix(spec$random$genomic$markers))
  expect_null(spec$random$genomic$ginv)
  expect_match(spec$bridge$target, "genomic_relationship_inverse", fixed = TRUE)
})

test_that("marker-based genomic builds a markers bridge payload", {
  ids <- paste0("g", 1:5)
  set.seed(4)
  M <- matrix(stats::rbinom(5 * 20, 2, 0.3), 5, 20)
  rownames(M) <- ids
  dat <- data.frame(y = c(1, 2, 3, 4, 5), id = ids)

  spec <- hsquared:::hs_build_model_spec(
    y ~ genomic(1 | id, markers = M),
    data = dat,
    family = stats::gaussian(),
    REML = TRUE
  )
  payload <- hsquared:::hs_build_bridge_payload(spec)

  expect_equal(payload$relationship_source, "markers")
  expect_null(payload$Ginv)
  expect_true(is.matrix(payload$markers))
  expect_equal(dim(payload$markers), c(5L, 20L))
  expect_equal(payload$ridge, 0.01)
  expect_equal(
    payload$relationship_provenance[c(
      "relationship_source",
      "relationship_method",
      "allele_frequency_source",
      "ridge",
      "relationship_scale"
    )],
    list(
      relationship_source = "markers",
      relationship_method = "vanraden1",
      allele_frequency_source = "sample",
      ridge = 0.01,
      relationship_scale = "K_lambda"
    )
  )
})

test_that("supplied-Ginv genomic still builds a supplied bridge payload", {
  ids <- paste0("g", 1:3)
  Ginv <- hs_test_ginv(ids)
  dat <- data.frame(y = c(1, 2, 3), id = ids)

  spec <- hsquared:::hs_build_model_spec(
    y ~ genomic(1 | id, Ginv = Ginv),
    data = dat,
    family = stats::gaussian(),
    REML = TRUE
  )
  payload <- hsquared:::hs_build_bridge_payload(spec)

  expect_equal(payload$relationship_source, "supplied_Ginv")
  expect_null(payload$markers)
  expect_true(is.matrix(payload$Ginv))
  expect_true(is.na(payload$ridge))
  expect_equal(
    payload$relationship_provenance$relationship_scale,
    "inverse_of_supplied_precision"
  )
  expect_true(is.na(payload$relationship_provenance$relationship_method))
  expect_true(is.na(payload$relationship_provenance$allele_frequency_source))
  expect_true(is.na(payload$relationship_provenance$ridge))
  expect_true(is.na(payload$relationship_provenance$scale_denominator))
})

test_that("genomic() takes exactly one of Ginv or markers", {
  ids <- paste0("g", 1:3)
  Ginv <- hs_test_ginv(ids)
  M <- matrix(0, 3, 5)
  rownames(M) <- ids
  dat <- data.frame(y = c(1, 2, 3), id = ids)
  expect_error(
    hsquared:::hs_build_model_spec(
      y ~ genomic(1 | id, Ginv = Ginv, markers = M),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "exactly one of `Ginv` or `markers`",
    fixed = TRUE
  )
})

test_that("genomic() marker ids must cover the data ids", {
  ids <- paste0("g", 1:3)
  M <- matrix(c(0, 1, 2, 1, 0, 2, 2, 1, 0, 0, 2, 1, 1, 2, 0), 3, 5)
  rownames(M) <- ids
  dat <- data.frame(y = c(1, 2), id = c("g1", "ghost"))
  expect_error(
    hsquared:::hs_build_model_spec(
      y ~ genomic(1 | id, markers = M),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "not in the `markers`",
    fixed = TRUE
  )
})

test_that("marker validation freezes dosage, names, and polymorphism", {
  ids <- paste0("g", 1:3)
  M <- matrix(c(0, 1, 2, 0, 2, 1), 3, 2)
  rownames(M) <- ids

  # Column names remain optional for backward compatibility.
  expect_equal(hsquared:::hs_validate_genomic_markers(M), M)

  bad <- M
  bad[1, 1] <- NA_real_
  expect_error(
    hsquared:::hs_validate_genomic_markers(bad),
    "finite, non-missing",
    fixed = TRUE
  )
  bad <- M
  bad[1, 1] <- 2.1
  expect_error(
    hsquared:::hs_validate_genomic_markers(bad),
    "[0, 2]",
    fixed = TRUE
  )
  # The dosage-range guard is a classed hsquared_error (#216), naming the
  # 0/1/2 convention and the -1/0/1 recode, not a bare stop().
  cnd <- tryCatch(
    hsquared:::hs_validate_genomic_markers(bad),
    error = function(e) e
  )
  expect_s3_class(cnd, "hsquared_error")
  expect_match(cnd$message, "recoded", fixed = TRUE)
  bad <- M
  rownames(bad)[1] <- ""
  expect_error(
    hsquared:::hs_validate_genomic_markers(bad),
    "nonempty row names",
    fixed = TRUE
  )
  bad <- M
  colnames(bad) <- c("m", "m")
  expect_error(
    hsquared:::hs_validate_genomic_markers(bad),
    "column names",
    fixed = TRUE
  )
  monomorphic <- matrix(0, 3, 4, dimnames = list(ids, NULL))
  expect_error(
    hsquared:::hs_validate_genomic_markers(monomorphic),
    "denominator `k` must be positive",
    fixed = TRUE
  )
})

test_that("supplied genomic precision is finite, symmetric, and positive definite", {
  ids <- paste0("g", 1:3)
  Q <- hs_test_ginv(ids)
  expect_equal(hsquared:::hs_validate_genomic_ginv(Q), Q)

  bad <- Q
  bad[1, 2] <- bad[1, 2] + 0.1
  expect_error(
    hsquared:::hs_validate_genomic_ginv(bad),
    "must be symmetric",
    fixed = TRUE
  )
  bad <- diag(c(1, 1, -1))
  dimnames(bad) <- list(ids, ids)
  expect_error(
    hsquared:::hs_validate_genomic_ginv(bad),
    "positive definite",
    fixed = TRUE
  )
  bad <- Q
  bad[1, 1] <- Inf
  expect_error(
    hsquared:::hs_validate_genomic_ginv(bad),
    "finite, non-missing",
    fixed = TRUE
  )
})

test_that("genomic payload preserves repeated records and extra genotyped ids", {
  ids <- paste0("g", 1:4)
  M <- matrix(
    c(
      0,
      1,
      2,
      0,
      2,
      1,
      0,
      1,
      0,
      2,
      1,
      2
    ),
    4,
    3
  )
  rownames(M) <- ids
  dat <- data.frame(
    y = c(1, 2, 3, 4),
    id = c("g2", "g1", "g2", "g1")
  )
  spec <- hsquared:::hs_build_model_spec(
    y ~ genomic(1 | id, markers = M),
    data = dat,
    family = stats::gaussian(),
    REML = TRUE
  )
  payload <- hsquared:::hs_build_bridge_payload(spec)

  expect_equal(payload$ids, ids)
  expect_equal(payload$metadata$observed_ids, dat$id)
  expect_equal(payload$metadata$observed_id_index, c(2L, 1L, 2L, 1L))
  expect_equal(dim(payload$Z), c(4L, 4L))
  expect_equal(
    unname(as.matrix(payload$Z)),
    rbind(c(0, 1, 0, 0), c(1, 0, 0, 0), c(0, 1, 0, 0), c(1, 0, 0, 0))
  )
})

test_that("genomic grammar rejects knobs and additional random effects", {
  ids <- paste0("g", 1:3)
  M <- matrix(c(0, 1, 2, 2, 1, 0), 3, 2, dimnames = list(ids, NULL))
  dat <- data.frame(y = 1:3, id = ids, pen = letters[1:3])

  expect_error(
    hsquared:::hs_build_model_spec(
      y ~ genomic(1 | id, markers = M, ridge = 0.02),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "`ridge`",
    fixed = TRUE
  )
  expect_error(
    hsquared:::hs_build_model_spec(
      y ~ genomic(1 | id, markers = M) + permanent(1 | pen),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "requires an `animal()` primary",
    fixed = TRUE
  )
  expect_error(
    hsquared(
      y ~ genomic(1 | id, markers = M),
      data = dat,
      family = stats::poisson()
    ),
    "not fitted on this path",
    fixed = TRUE
  )
})

test_that("genomic provenance normalizer enforces the frozen contract", {
  sha <- function(letter) paste(rep(letter, 64), collapse = "")
  marker <- list(
    relationship_source = "markers",
    relationship_method = "vanraden1",
    allele_frequency_source = "sample",
    ridge = 0.01,
    scale_denominator = 12.5,
    relationship_scale = "K_lambda",
    id_order_fingerprint = sha("a"),
    marker_content_fingerprint = sha("b"),
    kernel_fingerprint = sha("c"),
    precision_fingerprint = sha("d")
  )
  expect_equal(
    hsquared:::hs_normalize_genomic_provenance(marker),
    marker
  )

  mutated <- marker
  mutated$ridge <- 0.02
  expect_error(
    hsquared:::hs_normalize_genomic_provenance(mutated),
    "frozen VanRaden-1",
    fixed = TRUE
  )
  mutated <- marker
  mutated$kernel_fingerprint <- "not-a-sha"
  expect_error(
    hsquared:::hs_normalize_genomic_provenance(mutated),
    "lowercase SHA-256",
    fixed = TRUE
  )
  for (field in c("relationship_source", "relationship_method", "ridge")) {
    mutated <- marker
    mutated[[field]] <- if (identical(field, "ridge")) {
      NA_real_
    } else {
      NA_character_
    }
    expect_error(
      hsquared:::hs_normalize_genomic_provenance(mutated),
      "Internal bridge error",
      fixed = TRUE
    )
  }
  mutated <- marker
  mutated$relationship_source <- "bogus"
  expect_error(
    hsquared:::hs_normalize_genomic_provenance(mutated),
    "source must be exactly",
    fixed = TRUE
  )

  supplied <- list(
    relationship_source = "supplied_Ginv",
    relationship_method = NA_character_,
    allele_frequency_source = NA_character_,
    ridge = NA_real_,
    scale_denominator = NA_real_,
    relationship_scale = "inverse_of_supplied_precision",
    id_order_fingerprint = sha("e"),
    marker_content_fingerprint = NA_character_,
    kernel_fingerprint = NA_character_,
    precision_fingerprint = sha("f")
  )
  expect_equal(
    hsquared:::hs_normalize_genomic_provenance(supplied),
    supplied
  )
  mutated <- supplied
  mutated$ridge <- 0.01
  expect_error(
    hsquared:::hs_normalize_genomic_provenance(mutated),
    "must leave its construction",
    fixed = TRUE
  )
})

test_that("genomic boundary metadata distinguishes endpoint from MME representation", {
  contract <- hsquared:::hs_v07_genomic_boundary_contract()
  expect_identical(contract$doc46_commit, "fe96a147")
  expect_identical(
    contract$doc46_sha256,
    "283ab00bab3da925f0ac2916959efacaa7fb711c5da4dce09dd49ea568eef030"
  )
  expect_identical(
    contract$julia_implementation_commit,
    "ecc058f380be71058c9cfde373c345ab7a2f6aba"
  )
  expect_identical(contract$candidate_id, "v07_genomic_closed_boundary_v1")
  lower <- list(
    status = "boundary_lower",
    reason = "boundary_lower",
    profile_ratio = 0,
    numerical_ratio = 1e-7,
    boundary_epsilon = 1e-7,
    profile_loglik = -10,
    lower_derivative_per_observation = -0.1,
    upper_derivative_per_observation = -0.2
  )
  normalized <- hsquared:::hs_normalize_genomic_boundary(lower)
  expect_identical(normalized$status, "boundary_lower")
  expect_identical(normalized$profile_ratio, 0)
  expect_identical(normalized$numerical_ratio, 1e-7)

  upper <- lower
  upper$status <- upper$reason <- "boundary_upper"
  upper$profile_ratio <- 1
  upper$numerical_ratio <- 1 - 1e-7
  expect_identical(
    hsquared:::hs_normalize_genomic_boundary(upper)$profile_ratio,
    1
  )

  interior <- lower
  interior$status <- "interior"
  interior$reason <- "ai_interior"
  interior$profile_ratio <- 0.4
  interior$numerical_ratio <- 0.4
  expect_identical(
    hsquared:::hs_normalize_genomic_boundary(interior)$status,
    "interior"
  )

  unresolved <- lapply(lower, function(x) NA)
  unresolved$status <- "boundary_unresolved"
  unresolved$reason <- "endpoint_pair_tie"
  unresolved$boundary_epsilon <- 1e-7
  expect_true(is.na(
    hsquared:::hs_normalize_genomic_boundary(unresolved)$profile_ratio
  ))

  mutated <- lower
  mutated$boundary_epsilon <- 1e-6
  expect_error(
    hsquared:::hs_normalize_genomic_boundary(mutated),
    "epsilon drift"
  )
  mutated <- lower
  mutated$profile_ratio <- 1e-7
  expect_error(
    hsquared:::hs_normalize_genomic_boundary(mutated),
    "lower-boundary ratio contract drift"
  )
  mutated <- lower
  mutated$status <- "ordinary_convergence"
  expect_error(
    hsquared:::hs_normalize_genomic_boundary(mutated),
    "unknown genomic boundary status"
  )
})

test_that("genomic ratio keeps term compatibility and fences uncertainty", {
  fit <- structure(
    list(
      spec = list(target = "genomic"),
      result = list(
        heritability = data.frame(
          term = "genomic",
          component = "genomic_variance_ratio",
          estimate = 0.4,
          relationship_scale = "K_lambda",
          relationship_source = "markers",
          relationship_method = "vanraden1",
          allele_frequency_source = "sample",
          ridge = 0.01
        ),
        heritability_interval = data.frame(lower = 0.2, upper = 0.6),
        heritability_se = 0.1
      )
    ),
    class = "hsquared_fit"
  )

  out <- heritability(fit)
  expect_equal(out$term, "genomic")
  expect_equal(out$component, "genomic_variance_ratio")
  expect_equal(out$relationship_scale, "K_lambda")
  expect_error(
    heritability_interval(fit),
    "not available for genomic fits",
    fixed = TRUE
  )
  expect_error(
    heritability_standard_error(fit),
    "not available for genomic fits",
    fixed = TRUE
  )
})

test_that("the genomic target fixture pins VanRaden GBLUP and SNP-BLUP routes", {
  fixture <- testthat::test_path(
    "fixtures",
    "genomic_gblup_snpblup_target"
  )
  phenotypes <- utils::read.csv(
    file.path(fixture, "phenotypes.csv"),
    stringsAsFactors = FALSE
  )
  markers_df <- utils::read.csv(
    file.path(fixture, "markers.csv"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  frequencies <- utils::read.csv(
    file.path(fixture, "allele_frequencies.csv"),
    stringsAsFactors = FALSE
  )
  expected_g <- utils::read.csv(
    file.path(fixture, "expected_genomic_relationship.csv"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  expected_ginv <- utils::read.csv(
    file.path(fixture, "expected_genomic_precision.csv"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  expected_beta <- utils::read.csv(
    file.path(fixture, "expected_beta.csv"),
    stringsAsFactors = FALSE
  )
  expected_gebv <- utils::read.csv(
    file.path(fixture, "expected_gebv.csv"),
    stringsAsFactors = FALSE
  )
  expected_marker_effects <- utils::read.csv(
    file.path(fixture, "expected_marker_effects.csv"),
    stringsAsFactors = FALSE
  )
  metadata <- utils::read.csv(
    file.path(fixture, "expected_metadata.csv"),
    stringsAsFactors = FALSE
  )
  meta <- stats::setNames(metadata$value, metadata$key)

  ids <- markers_df$id
  markers <- as.matrix(markers_df[, -1, drop = FALSE])
  rownames(markers) <- ids
  storage.mode(markers) <- "numeric"

  expect_equal(phenotypes$id, ids)
  expect_equal(frequencies$marker, colnames(markers))
  expect_equal(expected_marker_effects$marker, colnames(markers))
  expect_equal(expected_gebv$id, ids)
  expect_equal(as.numeric(meta[["sigma_g2"]]), 2)
  expect_equal(as.numeric(meta[["sigma_e2"]]), 1)

  centered <- sweep(markers, 2, 2 * frequencies$frequency, "-")
  vanraden_scale <- 2 * sum(frequencies$frequency * (1 - frequencies$frequency))
  G <- tcrossprod(centered) / vanraden_scale
  rownames(G) <- colnames(G) <- ids

  G_expected <- as.matrix(expected_g[, -1, drop = FALSE])
  rownames(G_expected) <- expected_g$id
  storage.mode(G_expected) <- "numeric"
  expect_equal(rownames(G_expected), ids)
  expect_equal(colnames(G_expected), ids)
  expect_equal(G, G_expected, tolerance = 1e-12)
  expect_true(all(eigen(G, symmetric = TRUE, only.values = TRUE)$values > 0))
  expect_equal(as.numeric(meta[["k"]]), vanraden_scale, tolerance = 1e-12)

  Ginv <- solve(G)
  Ginv_expected <- as.matrix(expected_ginv[, -1, drop = FALSE])
  rownames(Ginv_expected) <- expected_ginv$id
  storage.mode(Ginv_expected) <- "numeric"
  expect_equal(rownames(Ginv_expected), ids)
  expect_equal(colnames(Ginv_expected), ids)
  expect_equal(Ginv, Ginv_expected, tolerance = 1e-12)
  expect_equal(unname(G %*% Ginv), diag(length(ids)), tolerance = 1e-12)

  X <- matrix(1, nrow(phenotypes), 1)
  sigma_ratio <- as.numeric(meta[["sigma_e2"]]) / as.numeric(meta[["sigma_g2"]])
  mme <- rbind(
    cbind(crossprod(X), t(X)),
    cbind(X, diag(length(ids)) + sigma_ratio * Ginv)
  )
  solution <- solve(mme, c(crossprod(X, phenotypes$y), phenotypes$y))
  beta <- unname(solution[1])
  gblup <- stats::setNames(unname(solution[-1]), ids)

  expect_equal(beta, expected_beta$value[1], tolerance = 1e-12)
  expect_equal(unname(gblup), expected_gebv$gblup, tolerance = 1e-12)
  expect_equal(
    max(abs(expected_gebv$gblup - expected_gebv$snp_blup)),
    as.numeric(meta[["gblup_snp_blup_max_abs_gebv_diff"]]),
    tolerance = 1e-12
  )

  snp_gebv <- centered %*% expected_marker_effects$effect
  expect_equal(as.numeric(snp_gebv), expected_gebv$snp_blup, tolerance = 1e-12)

  perturbed_g <- G_expected
  perturbed_g[1, 1] <- perturbed_g[1, 1] + 0.01
  expect_gt(max(abs(G - perturbed_g)), 0.001)
  perturbed_effects <- expected_marker_effects$effect
  perturbed_effects[1] <- perturbed_effects[1] + 0.01
  expect_gt(
    max(abs(
      as.numeric(centered %*% perturbed_effects) - expected_gebv$snp_blup
    )),
    0.001
  )
})

test_that("explicit marker and exact supplied-Q routes agree [live]", {
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for live marker GREML."
  )

  set.seed(6)
  na <- 10
  ids <- paste0("g", seq_len(na))
  M <- matrix(stats::rbinom(na * 80, 2, 0.3), na, 80)
  rownames(M) <- ids

  # Eight phenotyped individuals, four repeated, and two additional genotyped
  # individuals without records. Include a nonconstant fixed covariate so this
  # live identity gate exercises X and Z marshalling together.
  record_rows <- c(seq_len(8), seq_len(4))
  rec <- ids[record_rows]
  x <- seq(-1, 1, length.out = length(rec))
  dat <- data.frame(
    y = 3 + 0.4 * x + stats::rnorm(length(rec), 0, 1),
    x = x,
    id = rec
  )

  fit <- hsquared(
    y ~ x + genomic(1 | id, markers = M),
    data = dat,
    family = stats::gaussian(),
    control = hs_control(
      engine = "julia",
      engine_control = list(target = "genomic")
    )
  )

  expect_s3_class(fit, "hsquared_fit")
  expect_equal(fit$spec$target, "genomic")
  vc <- variance_components(fit)
  expect_equal(vc$component, c("genomic", "residual"))
  expect_true(all(is.finite(vc$estimate)) && all(vc$estimate > 0))
  expect_equal(nrow(breeding_values(fit)), na)
  expect_equal(breeding_values(fit)$id, ids)
  ratio <- heritability(fit)
  expect_equal(ratio$term, "genomic")
  expect_equal(ratio$component, "genomic_variance_ratio")
  expect_equal(ratio$relationship_source, "markers")
  expect_equal(ratio$relationship_method, "vanraden1")
  expect_equal(ratio$allele_frequency_source, "sample")
  expect_equal(ratio$ridge, 0.01)

  # Reuse the engine's exact frozen construction as a supplied precision. This
  # proves the two R routes reach the same supplied-Q estimator and the same
  # canonical precision fingerprint.
  JuliaCall::julia_assign("hsq_test_markers", unname(M))
  JuliaCall::julia_assign("hsq_test_ids", ids)
  JuliaCall::julia_command(paste(
    "hsq_test_activation = HSquared._genomic_activation_construction(",
    "hsq_test_markers, hsq_test_ids; marker_names=nothing, ridge=0.01);"
  ))
  Q <- JuliaCall::julia_eval("hsq_test_activation.Q")
  dimnames(Q) <- list(ids, ids)
  fit_q <- hsquared(
    y ~ x + genomic(1 | id, Ginv = Q),
    data = dat,
    family = stats::gaussian(),
    control = hs_control(
      engine = "julia",
      engine_control = list(target = "genomic")
    )
  )

  expect_equal(breeding_values(fit_q)$id, ids)

  expect_equal(
    variance_components(fit),
    variance_components(fit_q),
    tolerance = 1e-8
  )
  expect_equal(
    heritability(fit)$estimate,
    heritability(fit_q)$estimate,
    tolerance = 1e-8
  )
  expect_equal(fixef(fit), fixef(fit_q), tolerance = 1e-8)
  expect_equal(
    breeding_values(fit)$value,
    breeding_values(fit_q)$value,
    tolerance = 1e-8
  )
  expect_equal(
    fit$result$relationship_provenance$precision_fingerprint,
    fit_q$result$relationship_provenance$precision_fingerprint
  )
})

test_that("frozen activation fixture matches base R and both public routes [live]", {
  hs_skip_live_julia()
  project <- hsquared:::hs_default_julia_project()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(project),
    "JuliaCall, Julia, and local HSquared.jl are required for live fixture parity."
  )
  fixture <- file.path(
    project,
    "test",
    "fixtures",
    "genomic_public_activation_target"
  )
  testthat::skip_if_not(
    dir.exists(fixture),
    "The HSquared.jl genomic activation fixture is required."
  )

  marker_frame <- utils::read.csv(
    file.path(fixture, "markers.csv"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  ids <- marker_frame$id
  M <- as.matrix(marker_frame[, -1L, drop = FALSE])
  storage.mode(M) <- "double"
  rownames(M) <- ids
  dat <- utils::read.csv(
    file.path(fixture, "phenotypes.csv"),
    stringsAsFactors = FALSE
  )

  # Independent base-R construction: no package or Julia construction helper.
  p <- colMeans(M) / 2
  W <- sweep(M, 2L, 2 * p, "-")
  k <- 2 * sum(p * (1 - p))
  G <- tcrossprod(W) / k
  K <- G + diag(0.01, nrow(G))
  Q <- solve(K)
  dimnames(Q) <- list(ids, ids)

  hsquared:::hs_julia_setup(project)
  JuliaCall::julia_assign("hsq_fixture_markers", unname(M))
  JuliaCall::julia_assign("hsq_fixture_ids", ids)
  JuliaCall::julia_assign("hsq_fixture_marker_names", colnames(M))
  JuliaCall::julia_command(paste(
    "hsq_fixture_construction = HSquared._genomic_activation_construction(",
    "hsq_fixture_markers, hsq_fixture_ids;",
    "marker_names=hsq_fixture_marker_names, ridge=0.01);"
  ))
  Q_julia <- JuliaCall::julia_eval("hsq_fixture_construction.Q")
  expect_lte(max(abs(Q - Q_julia)), 1e-10)
  dimnames(Q_julia) <- list(ids, ids)
  genomic_control <- hs_control(
    engine = "julia",
    engine_control = list(target = "genomic")
  )

  fit_markers <- hsquared(
    y ~ x + genomic(1 | id, markers = M),
    data = dat,
    family = stats::gaussian(),
    control = genomic_control
  )
  fit_q_base_r <- hsquared(
    y ~ x + genomic(1 | id, Ginv = Q),
    data = dat,
    family = stats::gaussian(),
    control = genomic_control
  )
  fit_q <- hsquared(
    y ~ x + genomic(1 | id, Ginv = Q_julia),
    data = dat,
    family = stats::gaussian(),
    control = genomic_control
  )

  expect_equal(
    variance_components(fit_markers),
    variance_components(fit_q),
    tolerance = 1e-8
  )
  expect_equal(
    variance_components(fit_markers),
    variance_components(fit_q_base_r),
    tolerance = 1e-8
  )
  expect_equal(
    heritability(fit_markers)$estimate,
    heritability(fit_q)$estimate,
    tolerance = 1e-8
  )
  expect_equal(fixef(fit_markers), fixef(fit_q), tolerance = 1e-8)
  expect_gte(
    stats::cor(
      breeding_values(fit_markers)$value,
      breeding_values(fit_q)$value
    ),
    0.99999999
  )
  expect_lte(
    max(abs(
      breeding_values(fit_markers)$value - breeding_values(fit_q)$value
    )),
    1e-8
  )
  expect_equal(breeding_values(fit_markers)$id, ids)
  expect_equal(breeding_values(fit_q)$id, ids)
  expect_equal(
    fit_markers$result$relationship_provenance$precision_fingerprint,
    fit_q$result$relationship_provenance$precision_fingerprint
  )
})
