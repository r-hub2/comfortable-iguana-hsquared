hs_tiny_animal_validation_fixture <- function() {
  pedigree <- data.frame(
    id = c("calf", "sire", "dam"),
    sire = c("sire", NA, NA),
    dam = c("dam", NA, NA),
    stringsAsFactors = FALSE
  )
  data <- data.frame(
    y = c(1.2, 1.8, 2.6),
    age = c(0, 1, 2),
    id = c("sire", "dam", "calf"),
    stringsAsFactors = FALSE
  )

  formula <- y ~ age + animal(1 | id, pedigree = pedigree)
  ids <- c("sire", "dam", "calf")
  Ainv <- matrix(
    c(
      1.5,
      0.5,
      -1,
      0.5,
      1.5,
      -1,
      -1,
      -1,
      2
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(ids, ids)
  )

  list(
    name = "tiny_henderson_calf",
    description = paste(
      "Three-animal Henderson-style pedigree fixture for parser, bridge,",
      "and sparse Ainv validation."
    ),
    formula = formula,
    data = data,
    pedigree = pedigree,
    expected = list(
      ids = ids,
      sire_index = c(0L, 0L, 1L),
      dam_index = c(0L, 0L, 2L),
      Z = diag(3),
      Ainv = Ainv
    )
  )
}

hs_mrode9_pedigree_validation_fixture <- function() {
  if (!requireNamespace("nadiv", quietly = TRUE)) {
    stop(
      "The optional `nadiv` package is required for the Mrode9 validation ",
      "fixture.",
      call. = FALSE
    )
  }

  env <- new.env(parent = emptyenv())
  utils::data("Mrode9", package = "nadiv", envir = env)
  pedigree <- env$Mrode9[, c("pig", "sire", "dam")]
  names(pedigree) <- c("id", "sire", "dam")
  pedigree$id <- as.character(pedigree$id)
  pedigree$sire <- as.character(pedigree$sire)
  pedigree$dam <- as.character(pedigree$dam)
  pedigree$sire[is.na(pedigree$sire)] <- NA_character_
  pedigree$dam[is.na(pedigree$dam)] <- NA_character_

  ainv <- nadiv::makeAinv(pedigree)$Ainv
  colnames(ainv) <- rownames(ainv)

  list(
    name = "mrode9_nadiv_pedigree",
    description = paste(
      "Pedigree adapted from Mrode example 9.1 as shipped by nadiv; used",
      "for optional sparse Ainv comparator validation."
    ),
    source = paste(
      "nadiv::Mrode9, documented as adapted from example 9.1 of Mrode",
      "(2005), Linear Models for the Prediction of Animal Breeding Values."
    ),
    pedigree = pedigree,
    expected = list(Ainv = ainv)
  )
}

hs_henderson_mme_validation_fixture <- function() {
  ids <- c("founder_a", "founder_b", "animal_1", "animal_2", "animal_3")
  pedigree <- data.frame(
    id = ids,
    sire = c(NA, NA, "founder_a", "founder_a", "animal_1"),
    dam = c(NA, NA, "founder_b", "founder_b", "animal_2"),
    stringsAsFactors = FALSE
  )
  data <- data.frame(
    y = c(3.2, 4.1, 5.4, 5.9),
    x = c(0, 1, 0, 1),
    id = c("animal_1", "animal_2", "animal_3", "animal_3"),
    stringsAsFactors = FALSE
  )
  Ainv <- matrix(
    c(
      2,
      1,
      -1,
      -1,
      0,
      1,
      2,
      -1,
      -1,
      0,
      -1,
      -1,
      2.5,
      0.5,
      -1,
      -1,
      -1,
      0.5,
      2.5,
      -1,
      0,
      0,
      -1,
      -1,
      2
    ),
    nrow = 5,
    byrow = TRUE,
    dimnames = list(ids, ids)
  )

  list(
    name = "henderson_supplied_variance_mme",
    description = paste(
      "Five-animal supplied-variance Henderson mixed-model-equation",
      "fixture for fixed effects, EBVs, fitted values, and heritability."
    ),
    formula = y ~ x + animal(1 | id, pedigree = pedigree),
    data = data,
    pedigree = pedigree,
    sigma_a2 = 1.2,
    sigma_e2 = 0.8,
    expected = list(
      ids = ids,
      Ainv = Ainv,
      fixed_effects = c(
        "(Intercept)" = 3.898701298701298,
        x = 0.6454545454545471
      ),
      breeding_values = data.frame(
        id = ids,
        value = c(
          0,
          0,
          -0.054545454545454695,
          0.05454545454545385,
          0.8571428571428561
        ),
        stringsAsFactors = FALSE
      ),
      fitted = data.frame(
        .fitted = c(
          3.844155844155843,
          4.5987012987012985,
          4.755844155844154,
          5.401298701298701
        )
      ),
      heritability = 0.6
    )
  )
}

hs_reml_likelihood_validation_fixture <- function() {
  ids <- c("a", "b", "c")
  pedigree <- data.frame(
    id = ids,
    sire = c(NA, NA, NA),
    dam = c(NA, NA, NA),
    stringsAsFactors = FALSE
  )
  data <- data.frame(
    y = c(1, 2, 3),
    id = ids,
    stringsAsFactors = FALSE
  )

  list(
    name = "tiny_sparse_reml_likelihood_identity",
    description = paste(
      "Three-founder supplied-variance Gaussian likelihood fixture for",
      "dense REML, sparse REML, and ML identity checks."
    ),
    formula = y ~ animal(1 | id, pedigree = pedigree),
    data = data,
    pedigree = pedigree,
    sigma_a2 = 1,
    sigma_e2 = 1,
    expected = list(
      ids = ids,
      sire_index = c(0L, 0L, 0L),
      dam_index = c(0L, 0L, 0L),
      Z = diag(3),
      Ainv = diag(3),
      fixed_effects = c("(Intercept)" = 2),
      ml_loglik = -0.5 * (3 * log(2 * pi) + 3 * log(2) + 1),
      reml_loglik = -0.5 * (2 * log(2 * pi) + 3 * log(2) + log(1.5) + 1)
    )
  )
}

hs_mrode_supplied_variance_validation_fixture <- function() {
  ids <- as.character(c(2, 4, 1, 3, 5, 6, 7, 8, 9, 10, 11, 12))
  pedigree <- data.frame(
    id = ids,
    sire = c(NA, NA, NA, NA, "1", "3", "6", NA, "3", "3", "6", "6"),
    dam = c(NA, NA, NA, NA, "2", "4", "5", "5", "8", "8", "8", "8"),
    stringsAsFactors = FALSE
  )
  data <- data.frame(
    y = c(
      10.2,
      9.7,
      10.8,
      9.9,
      11.5,
      11.0,
      12.4,
      10.9,
      12.1,
      11.8,
      12.9,
      12.7
    ),
    x = rep(c(0, 1), 6),
    id = ids,
    stringsAsFactors = FALSE
  )
  Ainv <- matrix(
    c(
      1.5,
      0,
      0.5,
      0,
      -1,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      1.5,
      0,
      0.5,
      0,
      -1,
      0,
      0,
      0,
      0,
      0,
      0,
      0.5,
      0,
      1.5,
      0,
      -1,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0.5,
      0,
      2.5,
      0,
      -1,
      0,
      1,
      -1,
      -1,
      0,
      0,
      -1,
      0,
      -1,
      0,
      2.8333333333333335,
      0.5,
      -1,
      -0.6666666666666666,
      0,
      0,
      0,
      0,
      0,
      -1,
      0,
      -1,
      0.5,
      3.5,
      -1,
      1,
      0,
      0,
      -1,
      -1,
      0,
      0,
      0,
      0,
      -1,
      -1,
      2,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      1,
      -0.6666666666666666,
      1,
      0,
      3.3333333333333335,
      -1,
      -1,
      -1,
      -1,
      0,
      0,
      0,
      -1,
      0,
      0,
      0,
      -1,
      2,
      0,
      0,
      0,
      0,
      0,
      0,
      -1,
      0,
      0,
      0,
      -1,
      0,
      2,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      -1,
      0,
      -1,
      0,
      0,
      2,
      0,
      0,
      0,
      0,
      0,
      0,
      -1,
      0,
      -1,
      0,
      0,
      0,
      2
    ),
    nrow = 12,
    byrow = TRUE,
    dimnames = list(ids, ids)
  )

  list(
    name = "mrode_style_supplied_variance_outputs",
    description = paste(
      "Twelve-animal Mrode-style supplied-variance fixture for Ainv,",
      "Henderson MME outputs, ML/REML likelihoods, PEV, reliability,",
      "accuracy, and heritability."
    ),
    source = paste(
      "Mirrors the sibling HSquared.jl Phase 1 Mrode-style",
      "supplied-variance validation fixture."
    ),
    formula = y ~ x + animal(1 | id, pedigree = pedigree),
    data = data,
    pedigree = pedigree,
    sigma_a2 = 1.4,
    sigma_e2 = 0.9,
    expected = list(
      ids = ids,
      Ainv = Ainv,
      fixed_effects = c(
        "(Intercept)" = 11.317393070236822,
        x = -1.0063726022361354
      ),
      breeding_values = data.frame(
        id = ids,
        value = c(
          -0.5021061319436008,
          -0.11525433671959359,
          -0.13688874063925227,
          0.20787891031852276,
          0.1355094468877566,
          0.7022497098504573,
          0.7092602946040143,
          0.8873101719197427,
          0.6504124611509037,
          0.9594504746292134,
          1.1394542485192605,
          1.4922422619975693
        ),
        stringsAsFactors = FALSE
      ),
      fitted = data.frame(
        .fitted = c(
          10.815286938293221,
          10.195766131281093,
          11.18050432959757,
          10.518899378319208,
          11.452902517124578,
          11.013270177851144,
          12.026653364840836,
          11.198330639920428,
          11.967805531387725,
          11.270470942629899,
          12.456847318756083,
          11.803262729998256
        )
      ),
      prediction_error_variance = data.frame(
        id = ids,
        value = c(
          0.7995095200390083,
          0.7610908283670357,
          0.7995095200390082,
          0.7311017691367848,
          0.9080482784527696,
          0.7881609637019633,
          0.9046409311296871,
          0.7330130895894281,
          0.8217902414402057,
          0.8599627801418696,
          0.8487180314157624,
          0.8907596807820399
        ),
        stringsAsFactors = FALSE
      ),
      reliability = data.frame(
        id = ids,
        value = c(
          0.4289217714007082,
          0.4563636940235458,
          0.4289217714007084,
          0.47778445061658215,
          0.35139408681945017,
          0.43702788307002605,
          0.35382790633593764,
          0.4764192217218369,
          0.4130069703998529,
          0.385740871327236,
          0.3937728347030267,
          0.36374308515568576
        ),
        stringsAsFactors = FALSE
      ),
      heritability = 0.6086956521739131,
      ml_loglik = -18.181909573827813,
      reml_loglik = -16.973441618108648
    )
  )
}

# Published external canon: Mrode (2014), "Linear Models for the Prediction of
# Animal Breeding Values", 3rd ed., Example 3.1 (p.39) -- the textbook
# univariate animal model for pre-weaning gain (WWG). Unlike the
# `hs_mrode_supplied_variance_validation_fixture()` above (whose EBVs are
# regenerated by the package solver, so a wrong solver would still pass), this
# fixture pins the PUBLISHED textbook EBV digits. The inputs (8-animal pedigree,
# 5 records on animals 4-8, sex fixed effect coded 1 = male / 2 = female,
# sigma_a2 = 20 and sigma_e2 = 40 so the variance ratio alpha = 2) and the
# published solutions were confirmed against three independent citable sources
# (the masuday BLUPF90 tutorial citing Mrode 2014 p.39, the austin-putz Mrode
# chapter-3 R reproduction, and the Bioconductor `GeneticsPed` `Mrode3.1`
# dataset) and independently re-solved from the stated inputs (the published
# EBVs reproduce to ~5e-9). `Ainv` is built here by the tabular
# numerator-relationship method in pure base R (no `nadiv`) so the anchor is
# CI-runnable without Julia or optional packages; it is an INPUT to the solver,
# not the quantity under test.
hs_mrode_example_3_1_fixture <- function() {
  ids <- as.character(1:8)
  pedigree <- data.frame(
    id = ids,
    sire = c(NA, NA, NA, "1", "3", "1", "4", "3"),
    dam = c(NA, NA, NA, NA, "2", "2", "5", "6"),
    stringsAsFactors = FALSE
  )
  data <- data.frame(
    WWG = c(4.5, 2.9, 3.9, 3.5, 5.0),
    sex = c("male", "female", "female", "male", "male"),
    id = as.character(4:8),
    stringsAsFactors = FALSE
  )

  # Tabular numerator-relationship matrix A (parents precede offspring), then
  # Ainv = solve(A). Standard Henderson/Quaas recursion; 0 = unknown parent.
  sire_idx <- match(pedigree$sire, ids)
  dam_idx <- match(pedigree$dam, ids)
  sire_idx[is.na(sire_idx)] <- 0L
  dam_idx[is.na(dam_idx)] <- 0L
  n <- length(ids)
  A <- matrix(0, n, n)
  for (i in seq_len(n)) {
    s <- sire_idx[i]
    d <- dam_idx[i]
    A[i, i] <- 1 + if (s > 0L && d > 0L) 0.5 * A[s, d] else 0
    for (j in seq_len(i - 1L)) {
      aij <- 0
      if (s > 0L) {
        aij <- aij + 0.5 * A[j, s]
      }
      if (d > 0L) {
        aij <- aij + 0.5 * A[j, d]
      }
      A[i, j] <- aij
      A[j, i] <- aij
    }
  }
  Ainv <- solve(A)
  dimnames(Ainv) <- list(ids, ids)

  list(
    name = "mrode_example_3_1_published",
    description = paste(
      "Mrode (2014) Example 3.1 (p.39): published-EBV external canon for the",
      "univariate Gaussian animal model (sigma_a2 = 20, sigma_e2 = 40)."
    ),
    formula = WWG ~ sex + animal(1 | id, pedigree = pedigree),
    data = data,
    pedigree = pedigree,
    sigma_a2 = 20,
    sigma_e2 = 40,
    expected = list(
      ids = ids,
      Ainv = Ainv,
      # Published EBVs (Mrode 2014, 3rd ed., p.39), keyed to animals 1-8.
      breeding_values = stats::setNames(
        c(
          0.09844458,
          -0.01877010,
          -0.04108420,
          -0.00866312,
          -0.18573210,
          0.17687209,
          -0.24945855,
          0.18261469
        ),
        ids
      ),
      # The two sex solutions are parameterization-dependent (the sex block is
      # rank-deficient; the textbook uses a generalized inverse). The invariant
      # published quantity is the male - female contrast.
      sex_contrast_male_minus_female = 0.95407223,
      source = paste(
        "Mrode (2014) Linear Models for the Prediction of Animal Breeding",
        "Values, 3rd ed., Example 3.1 (p.39); inputs and solutions confirmed",
        "against the masuday BLUPF90 tutorial, the austin-putz reproduction,",
        "and Bioconductor GeneticsPed Mrode3.1, and independently re-solved."
      )
    )
  )
}

hs_mrode_example_5_1_multitrait_fixture <- function() {
  ids <- as.character(1:8)
  pedigree <- data.frame(
    id = ids,
    sire = c(NA, NA, NA, "1", "3", "1", "4", "3"),
    dam = c(NA, NA, NA, NA, "2", "2", "5", "6"),
    stringsAsFactors = FALSE
  )
  data <- data.frame(
    animal = as.character(4:8),
    sex = factor(c("1", "2", "2", "1", "1"), levels = c("1", "2")),
    WWG = c(4.5, 2.9, 3.9, 3.5, 5.0),
    PWG = c(6.8, 5.0, 6.8, 6.0, 7.5),
    stringsAsFactors = FALSE
  )

  sire_idx <- match(pedigree$sire, ids)
  dam_idx <- match(pedigree$dam, ids)
  sire_idx[is.na(sire_idx)] <- 0L
  dam_idx[is.na(dam_idx)] <- 0L
  n <- length(ids)
  A <- matrix(0, n, n)
  for (i in seq_len(n)) {
    s <- sire_idx[i]
    d <- dam_idx[i]
    A[i, i] <- 1 + if (s > 0L && d > 0L) 0.5 * A[s, d] else 0
    for (j in seq_len(i - 1L)) {
      aij <- 0
      if (s > 0L) {
        aij <- aij + 0.5 * A[j, s]
      }
      if (d > 0L) {
        aij <- aij + 0.5 * A[j, d]
      }
      A[i, j] <- aij
      A[j, i] <- aij
    }
  }
  Ainv <- solve(A)
  dimnames(Ainv) <- list(ids, ids)

  list(
    name = "mrode_example_5_1_multitrait_published",
    description = paste(
      "Mrode Example 5.1: published multiple-trait animal-model BLUP/MME",
      "anchor for pre-weaning and post-weaning gain at supplied G0/R0."
    ),
    data = data,
    pedigree = pedigree,
    traits = c("WWG", "PWG"),
    ids = ids,
    G0 = matrix(c(20, 18, 18, 40), nrow = 2, byrow = TRUE),
    R0 = matrix(c(40, 11, 11, 30), nrow = 2, byrow = TRUE),
    expected = list(
      Ainv = Ainv,
      fixed_effects = c(
        "trait1.sex1" = 4.3609,
        "trait2.sex1" = 6.7999,
        "trait1.sex2" = 3.3973,
        "trait2.sex2" = 5.8803
      ),
      breeding_values = data.frame(
        id = rep(ids, each = 2L),
        trait = rep(c("WWG", "PWG"), times = length(ids)),
        value = c(
          0.15092,
          0.27960,
          -0.015393,
          -0.0076101,
          -0.078392,
          -0.17034,
          -0.010239,
          -0.012671,
          -0.27033,
          -0.47783,
          0.27581,
          0.51724,
          -0.31612,
          -0.47898,
          0.24376,
          0.39196
        ),
        stringsAsFactors = FALSE
      ),
      source = paste(
        "Mrode Example 5.1 as reproduced in the LUKE Multiple trait animal",
        "model note (Mrode 1996, pp. 80-84) and Masuda's BLUPF90 tutorial",
        "for Mrode (2014) Chapter 5 Example 5.1. This is a supplied-G0/R0",
        "BLUP/MME target, not variance-component estimation."
      )
    )
  )
}

# Mrode (2014) Example 3.2 (p.48): the SIRE model on the same WWG data as 3.1.
# Random sires 1, 3, 4 are related via a sire numerator-relationship matrix
# (sire 4's sire = sire 1, matching the 3.1 pedigree); sigma_s2 = 5, sigma_e2 =
# 55 (alpha = 11). This extends the published-EBV external canon to a SECOND
# model class (the sire model) at near-zero risk: it reuses the tabular A-method
# and the reference Henderson solver, and is CI-runnable in pure base R.
#
# Inputs + published solutions confirmed against the masuday BLUPF90 tutorial
# (Mrode 3.2, p.48) and the austin-putz chapter-3 reproduction, and
# independently re-solved (~1e-7). Unlike 3.1, the sex block here is full rank
# (two means under a no-intercept design), so both sex means are published and
# directly reproducible.
hs_mrode_example_3_2_sire_fixture <- function() {
  data <- data.frame(
    WWG = c(4.5, 2.9, 3.9, 3.5, 5.0),
    sex = c("male", "female", "female", "male", "male"),
    sire = c("1", "3", "1", "4", "3"),
    stringsAsFactors = FALSE
  )
  sire_ids <- c("1", "3", "4")

  # Sire numerator-relationship matrix on {1, 3, 4} (sire 1, sire 3 founders;
  # sire 4's sire = sire 1, dam unknown). Same tabular recursion as 3.1.
  sire_ped_sire <- c(0L, 0L, 1L) # index into sire_ids; 0 = unknown
  sire_ped_dam <- c(0L, 0L, 0L)
  n <- length(sire_ids)
  A <- matrix(0, n, n)
  for (i in seq_len(n)) {
    s <- sire_ped_sire[i]
    d <- sire_ped_dam[i]
    A[i, i] <- 1 + if (s > 0L && d > 0L) 0.5 * A[s, d] else 0
    for (j in seq_len(i - 1L)) {
      aij <- 0
      if (s > 0L) {
        aij <- aij + 0.5 * A[j, s]
      }
      if (d > 0L) {
        aij <- aij + 0.5 * A[j, d]
      }
      A[i, j] <- aij
      A[j, i] <- aij
    }
  }
  Ainv <- solve(A)
  dimnames(Ainv) <- list(sire_ids, sire_ids)

  list(
    name = "mrode_example_3_2_sire_published",
    description = paste(
      "Mrode (2014) Example 3.2 (p.48): published-solution external canon for",
      "the SIRE model (related sires; sigma_s2 = 5, sigma_e2 = 55, alpha = 11)."
    ),
    data = data,
    sire_ids = sire_ids,
    sigma_s2 = 5,
    sigma_e2 = 55,
    expected = list(
      ids = sire_ids,
      Ainv = Ainv,
      # Published sire solutions (Mrode 2014, 3rd ed., p.48), keyed to sires.
      sire_solutions = stats::setNames(
        c(0.02200220, 0.01402640, -0.04304180),
        sire_ids
      ),
      # Published fixed-effect (sex) means under the no-intercept design.
      sex_male = 4.33567107,
      sex_female = 3.38198579,
      sex_contrast_male_minus_female = 4.33567107 - 3.38198579,
      source = paste(
        "Mrode (2014) Linear Models for the Prediction of Animal Breeding",
        "Values, 3rd ed., Example 3.2 (p.48); inputs and solutions confirmed",
        "against the masuday BLUPF90 tutorial and the austin-putz chapter-3",
        "reproduction, and independently re-solved (~1e-7)."
      )
    )
  )
}

hs_solve_henderson_mme_reference <- function(
  y,
  X,
  Z,
  Ainv,
  sigma_a2,
  sigma_e2,
  ids
) {
  if (!is.numeric(y)) {
    stop("`y` must be numeric.", call. = FALSE)
  }
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
  if (!inherits(Z, "Matrix")) {
    Z <- Matrix::Matrix(Z, sparse = TRUE)
  }
  if (!is.matrix(Ainv)) {
    Ainv <- as.matrix(Ainv)
  }
  if (!is.numeric(sigma_a2) || length(sigma_a2) != 1L || sigma_a2 <= 0) {
    stop("`sigma_a2` must be a positive number.", call. = FALSE)
  }
  if (!is.numeric(sigma_e2) || length(sigma_e2) != 1L || sigma_e2 <= 0) {
    stop("`sigma_e2` must be a positive number.", call. = FALSE)
  }

  residual_precision <- 1 / sigma_e2
  relationship_precision <- 1 / sigma_a2
  X <- Matrix::Matrix(X, sparse = TRUE)
  Ainv <- Matrix::Matrix(Ainv, sparse = TRUE)
  lhs <- rbind(
    cbind(
      residual_precision * Matrix::crossprod(X),
      residual_precision * Matrix::crossprod(X, Z)
    ),
    cbind(
      residual_precision * Matrix::crossprod(Z, X),
      residual_precision * Matrix::crossprod(Z) + relationship_precision * Ainv
    )
  )
  rhs <- c(
    as.numeric(residual_precision * Matrix::crossprod(X, y)),
    as.numeric(residual_precision * Matrix::crossprod(Z, y))
  )
  solution <- as.numeric(solve(as.matrix(lhs), rhs))
  nfixed <- ncol(X)
  fixed <- solution[seq_len(nfixed)]
  names(fixed) <- colnames(X)
  animal <- solution[-seq_len(nfixed)]

  list(
    fixed_effects = fixed,
    breeding_values = data.frame(
      id = ids,
      value = animal,
      stringsAsFactors = FALSE
    ),
    fitted = data.frame(
      .fitted = as.numeric(as.matrix(X) %*% fixed + as.matrix(Z) %*% animal)
    ),
    prediction_error_variance = data.frame(
      id = ids,
      value = unname(
        hs_henderson_pev_reference(X, Z, Ainv, sigma_a2, sigma_e2)
      ),
      stringsAsFactors = FALSE
    ),
    reliability = data.frame(
      id = ids,
      value = unname(
        hs_henderson_reliability_reference(
          X,
          Z,
          Ainv,
          sigma_a2,
          sigma_e2
        )
      ),
      stringsAsFactors = FALSE
    ),
    heritability = sigma_a2 / (sigma_a2 + sigma_e2)
  )
}

hs_solve_multivariate_henderson_mme_reference <- function(
  Y,
  sex,
  animal,
  ids,
  Ainv,
  G0,
  R0,
  traits = colnames(Y)
) {
  Y <- as.matrix(Y)
  if (!is.numeric(Y)) {
    stop("`Y` must be numeric.", call. = FALSE)
  }
  if (nrow(Y) != length(sex) || nrow(Y) != length(animal)) {
    stop(
      "`Y`, `sex`, and `animal` must have the same record count.",
      call. = FALSE
    )
  }
  if (!is.matrix(Ainv)) {
    Ainv <- as.matrix(Ainv)
  }
  if (!is.matrix(G0)) {
    G0 <- as.matrix(G0)
  }
  if (!is.matrix(R0)) {
    R0 <- as.matrix(R0)
  }
  n_traits <- ncol(Y)
  n_records <- nrow(Y)
  n_animals <- length(ids)
  if (!identical(dim(G0), c(n_traits, n_traits))) {
    stop("`G0` must be a square trait covariance matrix.", call. = FALSE)
  }
  if (!identical(dim(R0), c(n_traits, n_traits))) {
    stop("`R0` must be a square residual covariance matrix.", call. = FALSE)
  }

  sex <- factor(sex)
  sex_levels <- levels(sex)
  y <- as.numeric(t(Y))
  X <- matrix(0, n_records * n_traits, length(sex_levels) * n_traits)
  Z <- matrix(0, n_records * n_traits, n_animals * n_traits)

  fixed_names <- character(length(sex_levels) * n_traits)
  for (level_idx in seq_along(sex_levels)) {
    for (trait_idx in seq_len(n_traits)) {
      col <- (level_idx - 1L) * n_traits + trait_idx
      fixed_names[[col]] <- paste0(
        "trait",
        trait_idx,
        ".sex",
        sex_levels[[level_idx]]
      )
    }
  }
  colnames(X) <- fixed_names

  animal_match <- match(animal, ids)
  if (anyNA(animal_match)) {
    stop("All `animal` values must appear in `ids`.", call. = FALSE)
  }
  for (record_idx in seq_len(n_records)) {
    sex_idx <- match(as.character(sex[[record_idx]]), sex_levels)
    animal_idx <- animal_match[[record_idx]]
    for (trait_idx in seq_len(n_traits)) {
      row <- (record_idx - 1L) * n_traits + trait_idx
      X[row, (sex_idx - 1L) * n_traits + trait_idx] <- 1
      Z[row, (animal_idx - 1L) * n_traits + trait_idx] <- 1
    }
  }

  residual_precision <- kronecker(diag(n_records), solve(R0))
  relationship_precision <- kronecker(Ainv, solve(G0))
  lhs <- rbind(
    cbind(
      crossprod(X, residual_precision %*% X),
      crossprod(X, residual_precision %*% Z)
    ),
    cbind(
      crossprod(Z, residual_precision %*% X),
      crossprod(Z, residual_precision %*% Z) + relationship_precision
    )
  )
  rhs <- c(
    crossprod(X, residual_precision %*% y),
    crossprod(Z, residual_precision %*% y)
  )
  solution <- as.numeric(solve(lhs, rhs))
  nfixed <- ncol(X)
  fixed <- solution[seq_len(nfixed)]
  names(fixed) <- colnames(X)
  animal_effects <- matrix(
    solution[-seq_len(nfixed)],
    ncol = n_traits,
    byrow = TRUE,
    dimnames = list(ids, traits)
  )

  list(
    fixed_effects = fixed,
    breeding_values = data.frame(
      id = rep(ids, each = n_traits),
      trait = rep(traits, times = n_animals),
      value = as.numeric(t(animal_effects)),
      stringsAsFactors = FALSE
    )
  )
}

hs_henderson_pev_reference <- function(X, Z, Ainv, sigma_a2, sigma_e2) {
  inverse <- hs_henderson_mme_inverse_reference(
    X,
    Z,
    Ainv,
    sigma_a2,
    sigma_e2
  )
  nfixed <- ncol(as.matrix(X))
  diag(inverse[-seq_len(nfixed), -seq_len(nfixed), drop = FALSE])
}

hs_henderson_reliability_reference <- function(
  X,
  Z,
  Ainv,
  sigma_a2,
  sigma_e2
) {
  pev <- hs_henderson_pev_reference(X, Z, Ainv, sigma_a2, sigma_e2)
  relationship <- solve(as.matrix(Ainv))
  1 - pev / (sigma_a2 * diag(relationship))
}

hs_henderson_mme_inverse_reference <- function(
  X,
  Z,
  Ainv,
  sigma_a2,
  sigma_e2
) {
  residual_precision <- 1 / sigma_e2
  relationship_precision <- 1 / sigma_a2
  X <- Matrix::Matrix(X, sparse = TRUE)
  Z <- Matrix::Matrix(Z, sparse = TRUE)
  Ainv <- Matrix::Matrix(Ainv, sparse = TRUE)
  lhs <- rbind(
    cbind(
      residual_precision * Matrix::crossprod(X),
      residual_precision * Matrix::crossprod(X, Z)
    ),
    cbind(
      residual_precision * Matrix::crossprod(Z, X),
      residual_precision * Matrix::crossprod(Z) + relationship_precision * Ainv
    )
  )
  solve(as.matrix(lhs))
}

hs_gaussian_loglik_reference <- function(
  y,
  X,
  Z,
  Ainv,
  sigma_a2,
  sigma_e2,
  method = c("REML", "ML")
) {
  method <- match.arg(method)
  y <- as.numeric(y)
  X <- as.matrix(X)
  Z <- as.matrix(Z)
  A <- solve(as.matrix(Ainv))
  V <- sigma_a2 * Z %*% A %*% t(Z) + sigma_e2 * diag(length(y))
  n <- length(y)
  p <- ncol(X)

  # Near the h2 -> 1 boundary (residual variance ~ 0), an overshooting optimizer
  # can drive V to be numerically singular / non-positive-definite. A raw chol()
  # error there would propagate out of optim() and abort the whole estimate, so
  # treat any failed factorization or solve as an inadmissible point: return
  # loglik = -Inf (the minimized negative log-likelihood becomes +Inf), which
  # makes Nelder-Mead retreat from the bad region. Well-conditioned inputs are
  # unaffected and reproduce the exact previous values.
  evaluate <- function() {
    cholV <- chol(V)
    Vinv_y <- backsolve(cholV, forwardsolve(t(cholV), y))
    Vinv_X <- backsolve(cholV, forwardsolve(t(cholV), X))
    XtVinvX <- crossprod(X, Vinv_X)
    beta <- solve(XtVinvX, crossprod(X, Vinv_y))
    residual <- y - X %*% beta
    quad <- drop(crossprod(
      residual,
      backsolve(
        cholV,
        forwardsolve(t(cholV), residual)
      )
    ))
    logdetV <- 2 * sum(log(diag(cholV)))

    if (identical(method, "ML")) {
      loglik <- -0.5 * (n * log(2 * pi) + logdetV + quad)
    } else {
      cholXtVinvX <- chol(XtVinvX)
      logdetXtVinvX <- 2 * sum(log(diag(cholXtVinvX)))
      loglik <- -0.5 * ((n - p) * log(2 * pi) + logdetV + logdetXtVinvX + quad)
    }

    list(loglik = loglik, beta = as.numeric(beta))
  }

  tryCatch(
    evaluate(),
    error = function(e) list(loglik = -Inf, beta = rep(NA_real_, p))
  )
}

# Independent pure-R REML/ML variance-component optimizer used only as a
# validation reference: it maximizes the same dense Gaussian objective as
# hs_gaussian_loglik_reference() over log-variances, with no Julia involvement.
hs_reml_estimate_reference <- function(
  y,
  X,
  Z,
  Ainv,
  method = c("REML", "ML"),
  initial = c(sigma_a2 = 1, sigma_e2 = 1)
) {
  method <- match.arg(method)
  objective <- function(log_theta) {
    -hs_gaussian_loglik_reference(
      y,
      X,
      Z,
      Ainv,
      exp(log_theta[[1]]),
      exp(log_theta[[2]]),
      method = method
    )$loglik
  }
  opt <- tryCatch(
    stats::optim(
      log(as.numeric(initial[c("sigma_a2", "sigma_e2")])),
      objective,
      method = "Nelder-Mead",
      control = list(reltol = 1e-10, maxit = 1000L)
    ),
    # Never propagate a raw chol()/solve() error: if the optimization cannot be
    # evaluated at all, report a non-zero convergence code with NA estimates so
    # callers can detect failure instead of aborting.
    error = function(e) {
      list(par = c(NA_real_, NA_real_), value = Inf, convergence = 99L)
    }
  )
  estimate <- exp(opt$par)
  names(estimate) <- c("sigma_a2", "sigma_e2")
  loglik <- -opt$value
  convergence <- opt$convergence
  # If the optimum landed on an inadmissible (non-PD V) point, the objective is
  # +Inf there; surface that as a non-zero convergence code rather than a
  # spuriously "converged" non-finite estimate.
  if (!is.finite(loglik) || any(!is.finite(estimate))) {
    convergence <- if (identical(convergence, 0L)) 99L else convergence
  }
  list(
    estimate = estimate,
    loglik = loglik,
    convergence = convergence,
    method = method
  )
}

# Published REML variance components for the gryphon birth-weight univariate
# animal model (BWT ~ 1 + animal), used to anchor an external published-estimate
# recovery check of hsquared's pure-R REML reference. Source: Wilson et al.
# (2010) "An ecologist's guide to the animal model", J. Anim. Ecol. 79:13-26, as
# reproduced in the Wild Animal Models tutorial. The gryphon population is a
# teaching/simulated dataset (shipped in the CRAN package `enhancer`); the
# maintainer should confirm the headline numbers against the paper before any
# promotion. This is an external-anchor cross-check of the R reference optimizer,
# NOT the production fit path and NOT the twin-owned V1-MRODE-FIT gate row.
hs_gryphon_published_reml <- function() {
  out <- c(sigma_a2 = 3.3954, sigma_e2 = 3.8286, h2 = 0.470)
  attr(out, "source") <-
    "Wilson et al. (2010) J. Anim. Ecol. 79:13-26 (gryphon BWT animal model)"
  out
}

# Deterministic replicated animal-model dataset used only as an
# external-comparator input (e.g. pedigreemm). Three records per animal on the
# twelve-animal Mrode pedigree, so the design is non-degenerate for general REML
# fitters. The response was simulated once with set.seed(1) and is embedded
# verbatim; variance components are estimated, not supplied, and no
# data-generating recovery is claimed.
hs_replicated_animal_comparator_fixture <- function() {
  base <- hs_mrode_supplied_variance_validation_fixture()
  pedigree <- base$pedigree
  lab <- as.character(pedigree$id)
  y <- c(
    8.75230561,
    8.6588463,
    10.49847711,
    10.13870972,
    11.66745306,
    11.12747953,
    9.98559258,
    11.25827271,
    10.08334876,
    12.3774171,
    13.16984579,
    9.60592911,
    10.12178172,
    10.94582723,
    9.34616046,
    8.83855092,
    11.33115325,
    10.72724487,
    11.60897358,
    11.6475063,
    10.63796564,
    10.33658089,
    10.51332638,
    9.97539137,
    11.00568231,
    12.84065887,
    12.49999764,
    11.54006692,
    12.11236757,
    10.52352949,
    12.11579871,
    13.47549853,
    10.73007964,
    9.91801026,
    12.49008738,
    11.39403834
  )
  data <- data.frame(
    y = y,
    x = rep(c(0, 1, 0), times = length(lab)),
    id = rep(lab, each = 3L),
    stringsAsFactors = FALSE
  )
  list(
    name = "replicated_animal_reml_comparator",
    description = paste(
      "Deterministic replicated animal-model dataset (three records per animal",
      "on the twelve-animal Mrode pedigree) used only as an external-comparator",
      "input for REML variance-component estimation. Variance components are",
      "estimated, not supplied, and no data-generating recovery is claimed."
    ),
    source = paste(
      "Response simulated once with set.seed(1) from the animal model on the",
      "Mrode pedigree (sigma_a2 = sigma_e2 = 1) and embedded verbatim for",
      "reproducibility."
    ),
    formula = y ~ x + animal(1 | id, pedigree = pedigree),
    data = data,
    pedigree = pedigree,
    expected = list(ids = lab, Ainv = base$expected$Ainv)
  )
}
