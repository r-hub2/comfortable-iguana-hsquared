# Diagonal-G multivariate structure: fixture-verified unpack parity + an
# engine-guarded live bridge leg. The R unpack and the diagonal-vs-unstructured
# LRT are checked WITHOUT a live engine by feeding the twin's serialized
# `structured_covariance_parity` target through `hs_normalize_multivariate_result()`
# (the twin's #61 / #42 diagonal payload contract); the live leg actually fits
# `genetic_structure = "diagonal"` through the bridge and is skipped unless a
# local Julia + HSquared.jl is available. Engine row V4-MV-REML (partial).

hs_scp_path <- function(dir, file) {
  testthat::test_path("fixtures", dir, file)
}

hs_scp_csv <- function(dir, file) {
  utils::read.csv(
    hs_scp_path(dir, file),
    stringsAsFactors = FALSE,
    na.strings = c("", "NA")
  )
}

hs_scp_matrix <- function(dir, file) {
  d <- hs_scp_csv(dir, file)
  out <- as.matrix(d[-1])
  storage.mode(out) <- "double"
  rownames(out) <- d[[1]]
  out
}

hs_scp_meta <- function(dir) {
  d <- hs_scp_csv(dir, "expected_metadata.csv")
  stats::setNames(d$value, d$key)
}

test_that("R consumes the shared diagonal-G multivariate parity fixture", {
  dir <- "structured_covariance_parity"
  ped <- hs_scp_csv(dir, "pedigree.csv")
  names(ped)[names(ped) == "animal"] <- "id"
  ped$sire[ped$sire == "0"] <- NA
  ped$dam[ped$dam == "0"] <- NA

  pheno <- hs_scp_csv(dir, "phenotypes.csv")
  G0 <- hs_scp_matrix(dir, "expected_genetic_covariance.csv")
  R0 <- hs_scp_matrix(dir, "expected_residual_covariance.csv")
  beta <- hs_scp_matrix(dir, "expected_beta.csv")
  h2 <- hs_scp_csv(dir, "expected_heritability.csv")
  ebv <- hs_scp_csv(dir, "expected_ebv.csv")
  metadata <- hs_scp_meta(dir)

  # The fixture's defining property: a DIAGONAL genetic covariance (the
  # off-diagonal genetic covariances are exactly zero, and diag(G0) are the
  # per-trait genetic variances the metadata records).
  expect_equal(G0[1, 2], 0)
  expect_equal(G0[2, 1], 0)
  expect_equal(
    diag(G0),
    c(
      as.numeric(metadata[["genetic_variance_trait1"]]),
      as.numeric(metadata[["genetic_variance_trait2"]])
    ),
    tolerance = 1e-10,
    ignore_attr = TRUE
  )
  expect_identical(unname(metadata[["genetic_structure"]]), "diagonal")
  expect_identical(as.integer(metadata[["n_genetic_params"]]), 2L)

  spec <- hsquared:::hs_build_model_spec(
    cbind(trait1, trait2) ~ x + animal(1 | animal, pedigree = ped),
    data = pheno,
    family = stats::gaussian(),
    REML = TRUE
  )
  payload <- hsquared:::hs_build_bridge_payload(spec)

  # The diagonal payload the bridge would marshal at the fixture's values.
  raw <- list(
    genetic_covariance = G0,
    residual_covariance = R0,
    genetic_correlation = stats::cov2cor(G0), # diagonal G0 -> identity
    residual_correlation = stats::cov2cor(R0),
    heritability = h2$h2,
    beta = beta,
    breeding_ids = ebv$animal,
    breeding_traits = c("trait1", "trait2"),
    breeding_values = as.matrix(ebv[c("trait1", "trait2")]),
    loglik = as.numeric(metadata[["loglik"]]),
    converged = identical(tolower(metadata[["converged"]]), "true"),
    iterations = as.integer(metadata[["iterations"]]),
    traits = c("trait1", "trait2"),
    genetic_structure = "diagonal",
    n_genetic_params = as.integer(metadata[["n_genetic_params"]])
  )
  result <- hsquared:::hs_normalize_multivariate_result(raw, payload)
  fit <- hsquared:::hs_new_fit(
    spec = list(
      method = "REML",
      family = list(family = "gaussian"),
      target = "multivariate"
    ),
    payload = payload,
    result = result
  )

  expect_equal(genetic_covariance(fit), G0, tolerance = 1e-10)
  expect_equal(residual_covariance(fit), R0, tolerance = 1e-10)
  expect_equal(heritability(fit)$estimate, h2$h2, tolerance = 1e-10)
  # Diagonal G0 implies an identity genetic correlation.
  expect_equal(
    unname(genetic_correlation(fit)),
    diag(2),
    tolerance = 1e-10
  )

  expected_ebv <- data.frame(
    id = rep(ebv$animal, times = 2L),
    trait = rep(c("trait1", "trait2"), each = nrow(ebv)),
    value = c(ebv$trait1, ebv$trait2),
    stringsAsFactors = FALSE
  )
  expect_equal(breeding_values(fit), expected_ebv, tolerance = 1e-10)
  expect_equal(as.numeric(stats::logLik(fit)), as.numeric(metadata[["loglik"]]))

  # The contract fields the structure LRT reads.
  expect_identical(fit$result$genetic_structure, "diagonal")
  expect_identical(fit$result$n_genetic_params, 2L)
  expect_identical(
    fit_diagnostics(fit)$value[
      fit_diagnostics(fit)$metric == "genetic_structure"
    ],
    "diagonal"
  )

  # When the engine omits n_genetic_params, the normalizer derives it from the
  # structure (diagonal = t); this fallback feeds the structure-LRT df.
  raw_no_np <- raw
  raw_no_np$n_genetic_params <- NULL
  expect_identical(
    hsquared:::hs_normalize_multivariate_result(
      raw_no_np,
      payload
    )$n_genetic_params,
    2L
  )
})

test_that("the live Julia bridge fits diagonal G and the structure LRT end-to-end", {
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for the live diagonal multivariate bridge."
  )

  dir <- "structured_covariance_parity"
  ped <- hs_scp_csv(dir, "pedigree.csv")
  names(ped)[names(ped) == "animal"] <- "id"
  ped$sire[ped$sire == "0"] <- NA
  ped$dam[ped$dam == "0"] <- NA
  pheno <- hs_scp_csv(dir, "phenotypes.csv")
  G0 <- hs_scp_matrix(dir, "expected_genetic_covariance.csv")
  R0 <- hs_scp_matrix(dir, "expected_residual_covariance.csv")
  h2 <- hs_scp_csv(dir, "expected_heritability.csv")

  form <- cbind(trait1, trait2) ~ x + animal(1 | animal, pedigree = ped)

  diag_fit <- hsquared(
    form,
    data = pheno,
    family = stats::gaussian(),
    REML = TRUE,
    control = hs_control(
      engine = "julia",
      engine_control = list(
        target = "multivariate",
        genetic_structure = "diagonal",
        iterations = 2000L
      )
    )
  )
  full_fit <- hsquared(
    form,
    data = pheno,
    family = stats::gaussian(),
    REML = TRUE,
    control = hs_control(
      engine = "julia",
      engine_control = list(target = "multivariate", iterations = 2000L)
    )
  )

  # The diagonal fit estimates the off-diagonal genetic covariance as exactly
  # zero and recovers the serialized diagonal target.
  Ghat <- genetic_covariance(diag_fit)
  expect_equal(Ghat[1, 2], 0, tolerance = 1e-8)
  expect_equal(diag(Ghat), diag(G0), tolerance = 1e-4, ignore_attr = TRUE)
  expect_equal(
    residual_covariance(diag_fit),
    R0,
    tolerance = 1e-4,
    ignore_attr = TRUE
  )
  expect_equal(heritability(diag_fit)$estimate, h2$h2, tolerance = 1e-4)
  expect_identical(diag_fit$result$genetic_structure, "diagonal")
  expect_identical(diag_fit$result$n_genetic_params, 2L)

  # Diagonal-in-unstructured is the interior null: df = t(t-1)/2 = 1, and the
  # unstructured fit cannot do worse than the diagonal one.
  lrt <- covariance_structure_lrt(diag_fit, full_fit)
  expect_equal(lrt$df, 1L)
  expect_false(lrt$boundary)
  expect_gte(lrt$statistic, -1e-6)
  expect_identical(lrt$constrained, "diagonal")
  expect_identical(lrt$full, "unstructured")
})
