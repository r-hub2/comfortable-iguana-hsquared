# Opt-in, EXPERIMENTAL supplied-relationship primary effects. `relmat(1 | id,
# K = K)` supplies a dense symmetric positive-definite relationship/covariance
# matrix K (ID-keyed via row/column names); the parser marshals the supplied
# inverse Kinv = solve(K). `precision(1 | id, Q = Q)` supplies the precision
# (inverse) directly. Both dispatch through the SAME covered supplied-
# relationship-inverse engine path as genomic()'s Ginv (fit_ai_reml on a
# relationship-inverse animal_model_spec). NOT covered, NOT the default, does not
# move public_covered_count.

# A small symmetric positive-definite relationship matrix, ID-keyed.
hs_test_relmat_K <- function(ids) {
  n <- length(ids)
  k <- diag(n)
  for (i in seq_len(n - 1L)) {
    k[i, i + 1L] <- k[i + 1L, i] <- 0.25
  }
  dimnames(k) <- list(ids, ids)
  k
}

# Additive relationship matrix A for a topologically ordered pedigree (founders
# first, NA parents unknown). Mrode's tabular recursion; used as the reduction
# anchor K = A_pedigree and Q = solve(A) = Ainv.
hs_test_pedigree_A <- function(ped) {
  ids <- as.character(ped$id)
  n <- length(ids)
  sire <- match(as.character(ped$sire), ids)
  dam <- match(as.character(ped$dam), ids)
  A <- matrix(0, n, n)
  for (i in seq_len(n)) {
    s <- sire[[i]]
    d <- dam[[i]]
    A[i, i] <- 1 + if (!is.na(s) && !is.na(d)) 0.5 * A[s, d] else 0
    for (j in seq_len(i - 1L)) {
      a_js <- if (!is.na(s)) A[j, s] else 0
      a_jd <- if (!is.na(d)) A[j, d] else 0
      A[i, j] <- A[j, i] <- 0.5 * (a_js + a_jd)
    }
  }
  dimnames(A) <- list(ids, ids)
  A
}

# --- Parser (pure-R, no Julia) -------------------------------------------------

test_that("relmat(1 | id, K = K) parses as a supplied-inverse primary effect", {
  ids <- paste0("a", 1:4)
  K <- hs_test_relmat_K(ids)
  dat <- data.frame(y = c(1, 2, 3, 4), id = ids)

  spec <- hsquared:::hs_build_model_spec(
    y ~ relmat(1 | id, K = K),
    data = dat,
    family = stats::gaussian(),
    REML = TRUE
  )

  expect_named(spec$random, "relmat")
  expect_equal(spec$random$relmat$type, "relmat")
  expect_equal(spec$random$relmat$group, "id")
  expect_equal(spec$random$relmat$relationship, "relmat")
  expect_equal(spec$random$relmat$source, "supplied")
  expect_equal(spec$random$relmat$ids, ids)
  # The parser marshals the supplied INVERSE (like genomic's Ginv).
  expect_equal(spec$random$relmat$ginv, solve(K))
  expect_match(spec$bridge$target, "Kinv", fixed = TRUE)
})

test_that("precision(1 | id, Q = Q) parses and supplies the inverse directly", {
  ids <- paste0("a", 1:3)
  K <- hs_test_relmat_K(ids)
  Q <- solve(K)
  dat <- data.frame(y = c(1, 2, 3), id = ids)

  spec <- hsquared:::hs_build_model_spec(
    y ~ precision(1 | id, Q = Q),
    data = dat,
    family = stats::gaussian(),
    REML = TRUE
  )

  expect_equal(spec$random$precision$type, "precision")
  expect_equal(spec$random$precision$relationship, "precision")
  expect_equal(spec$random$precision$source, "supplied")
  # The user supplies the precision (inverse) DIRECTLY; it is carried verbatim.
  expect_equal(spec$random$precision$ginv, Q)
  expect_match(spec$bridge$target, "Kinv", fixed = TRUE)
})

test_that("precision() accepts a Kinv = alias for Q", {
  ids <- paste0("a", 1:3)
  Q <- solve(hs_test_relmat_K(ids))
  dat <- data.frame(y = c(1, 2, 3), id = ids)

  spec <- hsquared:::hs_build_model_spec(
    y ~ precision(1 | id, Kinv = Q),
    data = dat,
    family = stats::gaussian(),
    REML = TRUE
  )
  expect_equal(spec$random$precision$ginv, Q)
})

test_that("relmat() requires a K argument", {
  ids <- paste0("a", 1:3)
  dat <- data.frame(y = c(1, 2, 3), id = ids)
  expect_error(
    hsquared:::hs_build_model_spec(
      y ~ relmat(1 | id),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "requires a `K`",
    fixed = TRUE
  )
})

test_that("precision() requires a Q (or Kinv) argument", {
  ids <- paste0("a", 1:3)
  dat <- data.frame(y = c(1, 2, 3), id = ids)
  expect_error(
    hsquared:::hs_build_model_spec(
      y ~ precision(1 | id),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "requires a `Q`",
    fixed = TRUE
  )
})

test_that("a non-symmetric K is rejected", {
  ids <- paste0("a", 1:3)
  K <- hs_test_relmat_K(ids)
  K[1, 2] <- 0.9 # break symmetry
  dat <- data.frame(y = c(1, 2, 3), id = ids)
  expect_error(
    hsquared:::hs_build_model_spec(
      y ~ relmat(1 | id, K = K),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "symmetric",
    fixed = TRUE
  )
})

test_that("a non-positive-definite K is rejected", {
  ids <- paste0("a", 1:3)
  # A symmetric matrix with a negative eigenvalue (not PD).
  K <- matrix(c(1, 2, 0, 2, 1, 0, 0, 0, 1), 3, 3)
  dimnames(K) <- list(ids, ids)
  dat <- data.frame(y = c(1, 2, 3), id = ids)
  expect_error(
    hsquared:::hs_build_model_spec(
      y ~ relmat(1 | id, K = K),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "positive definite",
    fixed = TRUE
  )
})

test_that("a non-finite K is rejected", {
  ids <- paste0("a", 1:3)
  K <- hs_test_relmat_K(ids)
  K[1, 1] <- Inf
  dat <- data.frame(y = c(1, 2, 3), id = ids)
  expect_error(
    hsquared:::hs_build_model_spec(
      y ~ relmat(1 | id, K = K),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "finite",
    fixed = TRUE
  )
})

test_that("a K without row/column names is rejected", {
  ids <- paste0("a", 1:3)
  K <- hs_test_relmat_K(ids)
  dimnames(K) <- NULL
  dat <- data.frame(y = c(1, 2, 3), id = ids)
  expect_error(
    hsquared:::hs_build_model_spec(
      y ~ relmat(1 | id, K = K),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "row/column names",
    fixed = TRUE
  )
})

test_that("relmat() ids must be present in the K dimnames", {
  ids <- paste0("a", 1:3)
  K <- hs_test_relmat_K(ids)
  dat <- data.frame(y = c(1, 2), id = c("a1", "ghost"))
  expect_error(
    hsquared:::hs_build_model_spec(
      y ~ relmat(1 | id, K = K),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "not in the `K`",
    fixed = TRUE
  )
})

test_that("a non-square K is rejected", {
  ids <- paste0("a", 1:3)
  # A 3x2 (non-square) matrix carrying id row names.
  K <- matrix(0, nrow = 3L, ncol = 2L)
  rownames(K) <- ids
  colnames(K) <- ids[1:2]
  dat <- data.frame(y = c(1, 2, 3), id = ids)
  expect_error(
    hsquared:::hs_build_model_spec(
      y ~ relmat(1 | id, K = K),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "square",
    fixed = TRUE
  )
})

test_that("a K whose row and column names differ is rejected", {
  ids <- paste0("a", 1:3)
  K <- hs_test_relmat_K(ids)
  # Valid, unique row names but mismatched (different) column names.
  colnames(K) <- paste0("b", seq_along(ids))
  dat <- data.frame(y = c(1, 2, 3), id = ids)
  expect_error(
    hsquared:::hs_build_model_spec(
      y ~ relmat(1 | id, K = K),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "row and column names must match",
    fixed = TRUE
  )
})

test_that("a non-numeric K is rejected", {
  ids <- paste0("a", 1:3)
  K <- matrix("x", nrow = 3L, ncol = 3L)
  dimnames(K) <- list(ids, ids)
  dat <- data.frame(y = c(1, 2, 3), id = ids)
  expect_error(
    hsquared:::hs_build_model_spec(
      y ~ relmat(1 | id, K = K),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "numeric matrix",
    fixed = TRUE
  )
})

test_that("a non-matrix K is rejected", {
  ids <- paste0("a", 1:3)
  # A bare numeric vector is not a matrix.
  K <- c(a1 = 1, a2 = 1, a3 = 1)
  dat <- data.frame(y = c(1, 2, 3), id = ids)
  expect_error(
    hsquared:::hs_build_model_spec(
      y ~ relmat(1 | id, K = K),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "numeric matrix",
    fixed = TRUE
  )
})

test_that("a malformed precision() Q (non-symmetric) is rejected", {
  ids <- paste0("a", 1:3)
  Q <- solve(hs_test_relmat_K(ids))
  Q[1, 2] <- Q[1, 2] + 0.5 # break symmetry
  dat <- data.frame(y = c(1, 2, 3), id = ids)
  expect_error(
    hsquared:::hs_build_model_spec(
      y ~ precision(1 | id, Q = Q),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "symmetric",
    fixed = TRUE
  )
})

test_that("relmat(1 | id, Kinv = X) fits through the direct-inverse path", {
  ids <- paste0("a", 1:3)
  Q <- solve(hs_test_relmat_K(ids))
  dat <- data.frame(y = c(1, 2, 3), id = ids)

  # relmat(Kinv = X) and relmat(Q = X) supply the inverse directly (no solve),
  # exactly like precision(Q = X). All three must produce the SAME ginv.
  spec_relmat_kinv <- hsquared:::hs_build_model_spec(
    y ~ relmat(1 | id, Kinv = Q),
    data = dat,
    family = stats::gaussian(),
    REML = TRUE
  )
  spec_relmat_q <- hsquared:::hs_build_model_spec(
    y ~ relmat(1 | id, Q = Q),
    data = dat,
    family = stats::gaussian(),
    REML = TRUE
  )
  spec_precision <- hsquared:::hs_build_model_spec(
    y ~ precision(1 | id, Q = Q),
    data = dat,
    family = stats::gaussian(),
    REML = TRUE
  )

  # The supplied inverse is carried verbatim (no solve()), matching precision().
  expect_equal(spec_relmat_kinv$random$relmat$ginv, Q)
  expect_equal(spec_relmat_q$random$relmat$ginv, Q)
  expect_equal(
    spec_relmat_kinv$random$relmat$ginv,
    spec_precision$random$precision$ginv
  )
})

test_that("a formula must contain exactly one primary effect (animal + relmat)", {
  ids <- paste0("a", 1:2)
  K <- hs_test_relmat_K(ids)
  ped <- data.frame(id = ids, sire = c(NA, NA), dam = c(NA, NA))
  dat <- data.frame(y = c(1, 2), id = ids)
  expect_error(
    hsquared:::hs_build_model_spec(
      y ~ animal(1 | id, pedigree = ped) + relmat(1 | id, K = K),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "one primary",
    fixed = TRUE
  )
})

test_that("relmat and precision are valid opt-in julia targets", {
  expect_equal(hsquared:::hs_validate_julia_target("relmat"), "relmat")
  expect_equal(hsquared:::hs_validate_julia_target("precision"), "precision")
})

test_that("the default engine = \"fit\" path rejects a relmat() formula", {
  ids <- paste0("a", 1:3)
  K <- hs_test_relmat_K(ids)
  dat <- data.frame(y = c(1, 2, 3), id = ids)
  expect_error(
    hsquared(
      y ~ relmat(1 | id, K = K),
      data = dat,
      family = stats::gaussian()
    ),
    "experimental and opt-in",
    fixed = TRUE
  )
})

test_that("the default engine = \"fit\" path rejects a precision() formula", {
  ids <- paste0("a", 1:3)
  Q <- solve(hs_test_relmat_K(ids))
  dat <- data.frame(y = c(1, 2, 3), id = ids)
  expect_error(
    hsquared(
      y ~ precision(1 | id, Q = Q),
      data = dat,
      family = stats::gaussian()
    ),
    "experimental and opt-in",
    fixed = TRUE
  )
})

# --- Bridge payload (pure-R) ---------------------------------------------------

test_that("relmat builds a supplied-relationship-inverse bridge payload", {
  ids <- paste0("a", 1:4)
  K <- hs_test_relmat_K(ids)
  # Records reference the ids out of order to check id alignment.
  rec <- c("a3", "a1", "a4", "a2", "a1")
  dat <- data.frame(y = c(5, 1, 8, 3, 2), id = rec)

  spec <- hsquared:::hs_build_model_spec(
    y ~ relmat(1 | id, K = K),
    data = dat,
    family = stats::gaussian(),
    REML = TRUE
  )
  payload <- hsquared:::hs_build_bridge_payload(spec)

  expect_equal(payload$relationship_source, "supplied")
  expect_equal(payload$relationship, "relmat")
  expect_null(payload$markers)
  # The supplied INVERSE is marshalled (provenance), not an estimate.
  expect_equal(payload$Ginv, unname(solve(K)))
  expect_equal(payload$ids, ids)
  # Z is a record incidence on the id columns, in id order.
  expect_equal(dim(payload$Z), c(length(rec), length(ids)))
  expect_equal(as.integer(payload$Z %*% seq_along(ids)), match(rec, ids))
})

test_that("precision carries the supplied precision verbatim in the payload", {
  ids <- paste0("a", 1:3)
  Q <- solve(hs_test_relmat_K(ids))
  dat <- data.frame(y = c(1, 2, 3), id = ids)

  spec <- hsquared:::hs_build_model_spec(
    y ~ precision(1 | id, Q = Q),
    data = dat,
    family = stats::gaussian(),
    REML = TRUE
  )
  payload <- hsquared:::hs_build_bridge_payload(spec)

  expect_equal(payload$relationship_source, "supplied")
  expect_equal(payload$relationship, "precision")
  expect_equal(payload$Ginv, unname(Q))
})

# --- Live reduction / NULL cells (skip-guarded; require Julia) ------------------

test_that("relmat with K = A_pedigree fits identically to the animal model", {
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for live fits."
  )

  set.seed(101)
  # Topologically ordered pedigree: two founders, then offspring.
  ped <- data.frame(
    id = paste0("a", 1:8),
    sire = c(NA, NA, "a1", "a1", "a3", "a3", "a5", "a5"),
    dam = c(NA, NA, "a2", "a2", "a4", "a4", "a6", "a6"),
    stringsAsFactors = FALSE
  )
  A <- hs_test_pedigree_A(ped)

  # Simulate with genuine genetic signal (correlated breeding values from A) so
  # the REML optimum is interior and well identified; a pure-noise response sits
  # at the h2 = 0 boundary, where the reduction is degenerate.
  bv <- as.numeric(t(chol(A)) %*% stats::rnorm(nrow(A)))
  names(bv) <- ped$id
  n <- 64
  rec <- rep(ped$id, length.out = n)
  dat <- data.frame(y = 3 + bv[rec] + stats::rnorm(n, 0, 0.7), id = rec)

  fit_animal <- hsquared(
    y ~ animal(1 | id, pedigree = ped),
    data = dat,
    family = stats::gaussian()
  )
  fit_relmat <- hsquared(
    y ~ relmat(1 | id, K = A),
    data = dat,
    family = stats::gaussian(),
    control = hs_control(
      engine = "julia",
      engine_control = list(target = "relmat")
    )
  )

  vc_animal <- variance_components(fit_animal)
  vc_relmat <- variance_components(fit_relmat)
  # Same estimand up to component labelling; compare the sorted estimates.
  expect_equal(
    sort(vc_relmat$estimate),
    sort(vc_animal$estimate),
    tolerance = 1e-4
  )
  # The animal ratio is renamed to the relmat relationship, not "animal".
  expect_true("relmat" %in% vc_relmat$component)
})

test_that("precision with Q = Ainv matches the animal model", {
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for live fits."
  )

  set.seed(102)
  ped <- data.frame(
    id = paste0("a", 1:8),
    sire = c(NA, NA, "a1", "a1", "a3", "a3", "a5", "a5"),
    dam = c(NA, NA, "a2", "a2", "a4", "a4", "a6", "a6"),
    stringsAsFactors = FALSE
  )
  A <- hs_test_pedigree_A(ped)
  Ainv <- solve(A)

  bv <- as.numeric(t(chol(A)) %*% stats::rnorm(nrow(A)))
  names(bv) <- ped$id
  n <- 64
  rec <- rep(ped$id, length.out = n)
  dat <- data.frame(y = 3 + bv[rec] + stats::rnorm(n, 0, 0.7), id = rec)

  fit_animal <- hsquared(
    y ~ animal(1 | id, pedigree = ped),
    data = dat,
    family = stats::gaussian()
  )
  fit_prec <- hsquared(
    y ~ precision(1 | id, Q = Ainv),
    data = dat,
    family = stats::gaussian(),
    control = hs_control(
      engine = "julia",
      engine_control = list(target = "precision")
    )
  )

  expect_equal(
    sort(variance_components(fit_prec)$estimate),
    sort(variance_components(fit_animal)$estimate),
    tolerance = 1e-4
  )
})

test_that("relmat with an arbitrary PD K fits and returns sane components", {
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for live fits."
  )

  set.seed(13)
  na <- 8
  ids <- paste0("a", seq_len(na))
  # An arbitrary symmetric positive-definite relationship matrix.
  B <- matrix(stats::rnorm(na * na), na, na)
  K <- crossprod(B) / na + diag(na)
  dimnames(K) <- list(ids, ids)

  n <- 32
  rec <- rep(ids, length.out = n)
  dat <- data.frame(y = 3 + stats::rnorm(n, 0, 1), id = rec)

  fit <- hsquared(
    y ~ relmat(1 | id, K = K),
    data = dat,
    family = stats::gaussian(),
    control = hs_control(
      engine = "julia",
      engine_control = list(target = "relmat")
    )
  )

  expect_s3_class(fit, "hsquared_fit")
  expect_equal(fit$spec$target, "relmat")
  vc <- variance_components(fit)
  expect_equal(vc$component, c("relmat", "residual"))
  expect_true(all(is.finite(vc$estimate)) && all(vc$estimate > 0))
  h2 <- heritability(fit)$estimate
  expect_true(is.finite(h2) && h2 > 0 && h2 < 1)
  expect_equal(nrow(breeding_values(fit)), na)
})

test_that("relmat with K = I reduces to an independent iid animal fit", {
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for live fits."
  )

  set.seed(14)
  na <- 8
  ids <- paste0("a", seq_len(na))
  # Identity relationship: relmat(K = I) is an iid random intercept on id.
  K <- diag(na)
  dimnames(K) <- list(ids, ids)

  # Independent anchor: an ALL-FOUNDERS pedigree (all parents unknown) has A = I,
  # so the COVERED default animal-model path fits the SAME iid model through a
  # DIFFERENT engine route (fit_animal_model on a pedigree Ainv, not the supplied-
  # inverse fit_ai_reml). This verifies the null reduction against a genuinely
  # independent covered path, not just against precision(Q = I).
  ped <- data.frame(
    id = ids,
    sire = rep(NA_character_, na),
    dam = rep(NA_character_, na),
    stringsAsFactors = FALSE
  )

  # Genuine id-level (iid) signal so the REML optimum is interior and identified.
  id_effect <- stats::rnorm(na)
  names(id_effect) <- ids
  n <- 64
  rec <- rep(ids, length.out = n)
  dat <- data.frame(
    y = 3 + id_effect[rec] + stats::rnorm(n, 0, 0.7),
    id = rec
  )

  fit_relmat <- hsquared(
    y ~ relmat(1 | id, K = K),
    data = dat,
    family = stats::gaussian(),
    control = hs_control(
      engine = "julia",
      engine_control = list(target = "relmat")
    )
  )
  # Covered default engine = "fit" path; A = I from the all-founders pedigree.
  fit_animal <- hsquared(
    y ~ animal(1 | id, pedigree = ped),
    data = dat,
    family = stats::gaussian()
  )

  vc_relmat <- variance_components(fit_relmat)
  vc_animal <- variance_components(fit_animal)
  expect_true(
    all(is.finite(vc_relmat$estimate)) && all(vc_relmat$estimate > 0)
  )
  # Same estimand up to component labelling; compare sorted estimates and h2.
  expect_equal(
    sort(vc_relmat$estimate),
    sort(vc_animal$estimate),
    tolerance = 1e-4
  )
  expect_equal(
    heritability(fit_relmat)$estimate,
    heritability(fit_animal)$estimate,
    tolerance = 1e-4
  )
})

test_that("the supplied matrix is provenance, never reported as an estimate", {
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for live fits."
  )

  set.seed(15)
  na <- 6
  ids <- paste0("a", seq_len(na))
  K <- hs_test_relmat_K(ids)

  n <- 24
  rec <- rep(ids, length.out = n)
  dat <- data.frame(y = 3 + stats::rnorm(n, 0, 1), id = rec)

  fit <- hsquared(
    y ~ relmat(1 | id, K = K),
    data = dat,
    family = stats::gaussian(),
    control = hs_control(
      engine = "julia",
      engine_control = list(target = "relmat")
    )
  )

  # The supplied inverse travels in the payload as provenance (verbatim), and the
  # variance components are estimated scalars, NOT the supplied matrix.
  expect_equal(fit$payload$Ginv, unname(solve(K)))
  vc <- variance_components(fit)
  expect_equal(nrow(vc), 2L)
  expect_true(is.numeric(vc$estimate) && length(vc$estimate) == 2L)
})
