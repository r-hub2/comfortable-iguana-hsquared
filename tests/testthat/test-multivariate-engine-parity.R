# A26: R<->engine element-wise multivariate parity at the k = 2 promotion fixture.
#
# Tolerances are PREDECLARED in
# docs/dev-log/check-log.d/2026-09-02-h2-a26-mv-bridge-parity-predeclaration.md
# and committed before the first run (A25 Rose blocker 1). They are absolute and
# asserted as max|delta| so the declared number means what the predeclaration
# says; a miss is a finding, not a tolerance to widen.
#
# What this file adds that test-multivariate.R does not: the existing fixture
# test injects the expected_*.csv values into a hand-built `raw` list and checks
# the R normalizer surfaces them -- it executes no Julia. The Julia lane's
# testset evaluates the supplied-covariance MME at the stored covariances and
# never re-runs the optimizer. Neither crosses the bridge. These tests fit the
# fixture THROUGH the bridge and confront the serialized target element-wise.
#
# Passing this file does NOT authorize the R multivariate covered flip: it
# discharges one of nine acceptance criteria. R multivariate is covered at validation scale (G10); this file still pins local parity only.

hs_a26_fixture <- function(file) {
  utils::read.csv(
    testthat::test_path("fixtures", "phase4_multitrait_parity", file),
    stringsAsFactors = FALSE,
    na.strings = c("", "NA")
  )
}

hs_a26_matrix <- function(file) {
  dat <- hs_a26_fixture(file)
  out <- as.matrix(dat[-1])
  storage.mode(out) <- "double"
  rownames(out) <- dat[[1]]
  out
}

hs_a26_pedigree <- function() {
  ped <- hs_a26_fixture("pedigree.csv")
  names(ped)[names(ped) == "animal"] <- "id"
  ped$sire[ped$sire == "0"] <- NA
  ped$dam[ped$dam == "0"] <- NA
  ped
}

# The engine's own default start is data-scaled: 0.5 * per-trait phenotypic
# variance (src/multivariate.jl). The bridge never reaches it -- it always sends
# identity -- so matching it is what makes the start-matched leg meaningful.
hs_a26_engine_start <- function(pheno) {
  phen <- c(stats::var(pheno$trait1), stats::var(pheno$trait2))
  list(G0 = diag(0.5 * phen), R0 = diag(0.5 * phen))
}

hs_a26_fit <- function(ped, pheno, initial = NULL) {
  engine_control <- list(target = "multivariate", iterations = 2000L)
  if (!is.null(initial)) {
    engine_control$initial <- initial
  }
  hsquared(
    cbind(trait1, trait2) ~ x + animal(1 | animal, pedigree = ped),
    data = pheno,
    family = stats::gaussian(),
    REML = TRUE,
    control = hs_control(engine = "julia", engine_control = engine_control)
  )
}

# EBVs in the fixture's animal order, as a q x 2 matrix.
hs_a26_ebv_matrix <- function(fit, animals) {
  bv <- breeding_values(fit)
  vapply(
    c("trait1", "trait2"),
    function(tr) {
      rows <- bv[bv$trait == tr, , drop = FALSE]
      rows$value[match(animals, rows$id)]
    },
    numeric(length(animals))
  )
}

hs_a26_skip_bridge <- function() {
  hs_require_bridge("A26 parity")
}

test_that("Tier A: the payload R sends matches the engine's construction rule", {
  # Julia-free. Re-derives the Julia testset's X / Y / Z construction
  # (test/runtests.jl "Phase 4 shared multi-trait parity fixture") independently
  # in R and requires the emitted payload to equal it element-wise.
  ped <- hs_a26_pedigree()
  pheno <- hs_a26_fixture("phenotypes.csv")
  ebv <- hs_a26_fixture("expected_ebv.csv")

  spec <- hsquared:::hs_build_model_spec(
    cbind(trait1, trait2) ~ x + animal(1 | animal, pedigree = ped),
    data = pheno,
    family = stats::gaussian(),
    REML = TRUE
  )
  payload <- hsquared:::hs_build_bridge_payload(spec)

  # The fixture serializes EBVs in Julia's normalize_pedigree(...).ids order
  # (the Julia testset asserts `ebv_ids == ped.ids`). Requiring the R emitter's
  # id order to equal it pins the cross-lane animal ordering with no Julia at
  # runtime -- so a divergence between R's hs_topological_pedigree and Julia's
  # _topological_order is caught on the CRAN lane too.
  expect_identical(as.character(payload$ids), as.character(ebv$animal))

  # X == [1 x]; compare values only (model.matrix carries an `assign` attribute
  # the engine never sees).
  expect_equal(
    as.vector(payload$X),
    as.vector(cbind(1, pheno$x)),
    tolerance = 0
  )
  expect_equal(dim(payload$X), c(nrow(pheno), 2L))
  expect_equal(
    as.vector(payload$Y),
    as.vector(cbind(pheno$trait1, pheno$trait2)),
    tolerance = 0
  )

  z_expected <- matrix(0, nrow(pheno), length(payload$ids))
  z_expected[cbind(
    seq_len(nrow(pheno)),
    match(pheno$animal, payload$ids)
  )] <- 1
  expect_equal(
    as.vector(as.matrix(payload$Z)),
    as.vector(z_expected),
    tolerance = 0
  )
  # One animal per record: an incidence matrix, not a covariate block.
  expect_true(all(Matrix::rowSums(payload$Z) == 1))
})

test_that("Tier B1/B2: the bridge reproduces the serialized k=2 target element-wise", {
  hs_a26_skip_bridge()

  ped <- hs_a26_pedigree()
  pheno <- hs_a26_fixture("phenotypes.csv")
  G0 <- hs_a26_matrix("expected_genetic_covariance.csv")
  R0 <- hs_a26_matrix("expected_residual_covariance.csv")
  beta <- hs_a26_matrix("expected_beta.csv")
  h2 <- hs_a26_fixture("expected_heritability.csv")$h2
  ebv <- hs_a26_fixture("expected_ebv.csv")
  ebv_target <- as.matrix(ebv[c("trait1", "trait2")])
  meta_raw <- hs_a26_fixture("expected_metadata.csv")
  meta <- stats::setNames(meta_raw$value, meta_raw$key)
  loglik <- as.numeric(meta[["loglik"]])

  # ---- Tier B1: start-matched. atol 1e-8. -------------------------------
  # Same estimator, same data, same starting simplex, deterministic unseeded
  # Nelder-Mead -- so this is a marshalling/determinism pin. A transposed
  # matrix, a mis-ordered Z column, a broken sparse-CSC hand-off, or an
  # NA->NaN slip cannot survive it.
  fit <- hs_a26_fit(ped, pheno, initial = hs_a26_engine_start(pheno))
  expect_true(fit$result$converged)

  expect_lt(max(abs(genetic_covariance(fit) - G0)), 1e-8)
  expect_lt(max(abs(residual_covariance(fit) - R0)), 1e-8)
  expect_lt(max(abs(genetic_correlation(fit) - stats::cov2cor(G0))), 1e-8)
  expect_lt(max(abs(residual_correlation(fit) - stats::cov2cor(R0))), 1e-8)
  expect_lt(max(abs(heritability(fit)$estimate - h2)), 1e-8)
  expect_lt(
    max(abs(matrix(fixef(fit)$estimate, nrow = nrow(beta)) - beta)),
    1e-8
  )
  expect_lt(
    max(abs(hs_a26_ebv_matrix(fit, ebv$animal) - ebv_target)),
    1e-8
  )
  expect_lt(abs(fit$result$loglik - loglik), 1e-8)
  expect_lt(
    abs(
      stats::cov2cor(genetic_covariance(fit))[1, 2] -
        as.numeric(meta[["genetic_correlation_trait1_trait2"]])
    ),
    1e-8
  )

  # Structural contract (Tier A, live half).
  expect_equal(
    dimnames(genetic_covariance(fit)),
    list(c("trait1", "trait2"), c("trait1", "trait2"))
  )
  expect_equal(heritability(fit)$trait, c("trait1", "trait2"))
  expect_equal(nrow(breeding_values(fit)), nrow(ebv) * 2L)
  expect_equal(stats::nobs(fit), nrow(pheno) * 2L)
  expect_equal(attr(stats::logLik(fit), "df"), 2L * 2L + 2L * 3L)
  # The engine's iteration count survives the bridge (it is carried in
  # diagnostics, not at the top level of `result`) -- cap exhaustion must stay
  # visible to an R user, since the optimizer is derivative-free.
  diagnostics <- fit_diagnostics(fit)
  iterations <- diagnostics$value[diagnostics$metric == "iterations"]
  expect_equal(length(iterations), 1L)
  expect_gt(as.integer(iterations), 0L)

  # The recovered off-diagonal is genuine, so r_g is a real quantity here and
  # not a degenerate 0/1 that any implementation would reproduce.
  expect_gt(abs(stats::cov2cor(G0)[1, 2]), 0.05)
  expect_lt(abs(stats::cov2cor(G0)[1, 2]), 0.99)

  # ---- Tier B2: the bridge's OWN default start. --------------------------
  # hs_validate_multivariate_initial(NULL, t) sends identity, never the
  # engine's data-scaled default, so this leg measures what the shipped
  # default path actually recovers. Predeclared: 5e-4 covariance / h2 /
  # correlation / beta, 5e-3 EBV, 1e-3 loglik.
  fit_default <- hs_a26_fit(ped, pheno)
  expect_true(fit_default$result$converged)

  expect_lt(max(abs(genetic_covariance(fit_default) - G0)), 5e-4)
  expect_lt(max(abs(residual_covariance(fit_default) - R0)), 5e-4)
  expect_lt(
    max(abs(genetic_correlation(fit_default) - stats::cov2cor(G0))),
    5e-4
  )
  expect_lt(
    max(abs(residual_correlation(fit_default) - stats::cov2cor(R0))),
    5e-4
  )
  expect_lt(max(abs(heritability(fit_default)$estimate - h2)), 5e-4)
  expect_lt(
    max(abs(
      matrix(fixef(fit_default)$estimate, nrow = nrow(beta)) - beta
    )),
    5e-4
  )
  expect_lt(
    max(abs(hs_a26_ebv_matrix(fit_default, ebv$animal) - ebv_target)),
    5e-3
  )
  expect_lt(abs(fit_default$result$loglik - loglik), 1e-3)
})

test_that("Tier B3: the fit is invariant to pedigree row order", {
  hs_a26_skip_bridge()

  # The R emitter builds Z columns over R's hs_topological_pedigree order; Julia
  # rebuilds Ainv over its own normalize_pedigree order. Both are DFS post-order
  # (sire subtree, dam subtree, emit), so they agree by parallel implementation
  # and are verified by no test. The fixture pedigree is ALREADY topologically
  # sorted, so it cannot discriminate -- reversing and shuffling the rows forces
  # a genuine re-sort on both sides. A Z<->Ainv misalignment surfaces here and
  # nowhere else in the suite. Predeclared atol 1e-6.
  ped <- hs_a26_pedigree()
  pheno <- hs_a26_fixture("phenotypes.csv")
  G0 <- hs_a26_matrix("expected_genetic_covariance.csv")
  R0 <- hs_a26_matrix("expected_residual_covariance.csv")
  ebv <- hs_a26_fixture("expected_ebv.csv")
  ebv_target <- as.matrix(ebv[c("trait1", "trait2")])
  initial <- hs_a26_engine_start(pheno)

  ped_reversed <- ped[rev(seq_len(nrow(ped))), ]
  row.names(ped_reversed) <- NULL

  withr::with_seed(26, {
    ped_shuffled <- ped[sample.int(nrow(ped)), ]
  })
  row.names(ped_shuffled) <- NULL

  # Guard the guard: the permutation must actually change the emitted order,
  # otherwise this test would silently assert nothing.
  emitted_ids <- function(p) {
    hsquared:::hs_build_bridge_payload(hsquared:::hs_build_model_spec(
      cbind(trait1, trait2) ~ x + animal(1 | animal, pedigree = p),
      data = pheno,
      family = stats::gaussian(),
      REML = TRUE
    ))$ids
  }
  expect_false(identical(emitted_ids(ped_reversed), emitted_ids(ped)))
  expect_false(identical(emitted_ids(ped_shuffled), emitted_ids(ped)))

  for (permuted in list(ped_reversed, ped_shuffled)) {
    fit <- hs_a26_fit(permuted, pheno, initial = initial)
    expect_true(fit$result$converged)
    expect_lt(max(abs(genetic_covariance(fit) - G0)), 1e-6)
    expect_lt(max(abs(residual_covariance(fit) - R0)), 1e-6)
    # ID-matched: the EBV for animal `a` must be a's EBV whatever row a sat in.
    expect_lt(
      max(abs(hs_a26_ebv_matrix(fit, ebv$animal) - ebv_target)),
      1e-6
    )
  }
})
