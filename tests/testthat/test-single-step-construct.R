# Opt-in, experimental single-step H^-1 CONSTRUCTION bridge
# (docs/design/25-single-step-construction-bridge.md). The parser / payload /
# genotyped_rows-alignment tests run without Julia; the live tests are
# skip-guarded on the local HSquared.jl bridge. The independent correctness
# guards are the shuffled-marker REORDER test (§6.3) and the
# differs-from-pedigree-model anchor -- not a self-referential equivalence.

test_that("single_step construction parses + aligns genotyped_rows", {
  ped <- data.frame(
    id = c("c1", "s", "d", "c2", "c3"),
    sire = c("s", NA, NA, "s", "s"),
    dam = c("d", NA, NA, "d", "d"),
    stringsAsFactors = FALSE
  )
  dat <- data.frame(
    y = c(1.1, 2.2, 3.0, 2.5, 1.8),
    id = c("c1", "s", "d", "c2", "c3"),
    stringsAsFactors = FALSE
  )
  # markers for a genotyped subset (c1/c2/c3), rows deliberately out of order
  m <- matrix(
    c(0, 1, 2, 2, 1, 0, 1, 1, 1, 0, 2, 1),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("c3", "c1", "c2"), paste0("m", 1:4))
  )
  spec <- hsquared:::hs_build_model_spec(
    y ~ single_step(1 | id, pedigree = ped, markers = m),
    data = dat,
    family = stats::gaussian(),
    REML = TRUE
  )
  ss <- spec$random$single_step
  expect_equal(ss$source, "construct")
  # topological pedigree order: founders before offspring
  expect_equal(ss$ids, c("s", "d", "c1", "c2", "c3"))
  expect_equal(ss$genotyped_rows, c(3L, 4L, 5L))
  # marker rows reordered to the genotyped animals' pedigree-row order, and the
  # genotyped_rows index exactly those animals (the §3 alignment invariant).
  expect_equal(rownames(ss$markers), c("c1", "c2", "c3"))
  expect_equal(ss$ids[ss$genotyped_rows], rownames(ss$markers))

  payload <- hsquared:::hs_build_bridge_payload(spec)
  expect_equal(payload$relationship_source, "construct")
  expect_equal(payload$genotyped_rows, c(3L, 4L, 5L))
  expect_equal(payload$pedigree$id, c("s", "d", "c1", "c2", "c3"))
  expect_equal(dim(payload$markers), c(3L, 4L))
  expect_match(
    payload$metadata$julia_fit_target,
    "fit_single_step_reml",
    fixed = TRUE
  )
})

test_that("single_step construction aligns NON-contiguous genotyped animals", {
  # genotype a founder (s) and a mid animal (c1), scrambled marker rows; the
  # genotyped_rows are non-contiguous (1, 3) and must still index the right ids.
  ped <- data.frame(
    id = c("c1", "s", "d", "c2", "c3"),
    sire = c("s", NA, NA, "s", "s"),
    dam = c("d", NA, NA, "d", "d"),
    stringsAsFactors = FALSE
  )
  dat <- data.frame(y = c(1, 2, 3, 4, 5), id = ped$id)
  m <- matrix(
    c(2, 0, 1, 0, 1, 2),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("c1", "s"), c("m1", "m2", "m3"))
  )
  spec <- hsquared:::hs_build_model_spec(
    y ~ single_step(1 | id, pedigree = ped, markers = m),
    data = dat,
    family = stats::gaussian(),
    REML = TRUE
  )
  ss <- spec$random$single_step
  # ped order s,d,c1,c2,c3 -> s is row 1, c1 is row 3
  expect_equal(ss$genotyped_rows, c(1L, 3L))
  expect_equal(rownames(ss$markers), c("s", "c1"))
  expect_equal(ss$ids[ss$genotyped_rows], rownames(ss$markers))
})

test_that("single_step construction accepts ungenotyped phenotyped animals", {
  # the single-step point: observed (phenotyped) need NOT be a subset of markers.
  ped <- data.frame(
    id = c("c1", "s", "d", "c2"),
    sire = c("s", NA, NA, "s"),
    dam = c("d", NA, NA, "d"),
    stringsAsFactors = FALSE
  )
  dat <- data.frame(y = c(1, 2, 3, 4), id = c("c1", "s", "d", "c2"))
  # only c1, c2 genotyped; s, d phenotyped but ungenotyped
  m <- matrix(
    c(0, 1, 2, 1),
    nrow = 2,
    dimnames = list(c("c1", "c2"), c("m1", "m2"))
  )
  spec <- hsquared:::hs_build_model_spec(
    y ~ single_step(1 | id, pedigree = ped, markers = m),
    data = dat,
    family = stats::gaussian(),
    REML = TRUE
  )
  ss <- spec$random$single_step
  expect_equal(ss$ids, c("s", "d", "c1", "c2"))
  expect_equal(ss$genotyped_rows, c(3L, 4L)) # c1, c2 only
  # ungenotyped phenotyped animals (s, d) are not in genotyped_rows but stay in the
  # pedigree (they get GEBVs via the pedigree relationship).
  expect_false(any(c(1L, 2L) %in% ss$genotyped_rows))
})

test_that("single_step H^Gamma payload aligns group, Gamma, and genotyped_rows", {
  ped <- data.frame(
    id = c("c1", "s", "d", "c2", "c3"),
    sire = c("s", NA, NA, "s", "s"),
    dam = c("d", NA, NA, "d", "d"),
    stringsAsFactors = FALSE
  )
  dat <- data.frame(
    y = c(1.1, 2.2, 3.0, 2.5, 1.8),
    id = c("c1", "s", "d", "c2", "c3"),
    stringsAsFactors = FALSE
  )
  m <- matrix(
    c(0, 1, 2, 2, 1, 0, 1, 1, 1, 0, 2, 1),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("c3", "c1", "c2"), paste0("m", 1:4))
  )
  mf_group <- c(c1 = "", s = "base_s", d = "base_d", c2 = "", c3 = "")
  Gamma <- matrix(
    c(0.30, 0.05, 0.05, 0.25),
    nrow = 2,
    dimnames = list(c("base_d", "base_s"), c("base_d", "base_s"))
  )

  spec <- hsquared:::hs_build_model_spec(
    y ~ single_step(
      1 | id,
      pedigree = ped,
      markers = m,
      group = mf_group,
      Gamma = Gamma,
      ridge = 0.01
    ),
    data = dat,
    family = stats::gaussian(),
    REML = TRUE
  )
  ss <- spec$random$single_step
  expect_equal(ss$source, "metafounder_construct")
  expect_equal(ss$relationship, "metafounder_single_step")
  expect_equal(ss$ids, c("s", "d", "c1", "c2", "c3"))
  expect_equal(unname(ss$group_of), c("base_s", "base_d", "", "", ""))
  expect_equal(ss$gamma_labels, c("base_s", "base_d"))
  expect_equal(ss$Gamma, matrix(c(0.25, 0.05, 0.05, 0.30), nrow = 2))
  expect_equal(ss$genotyped_rows, c(3L, 4L, 5L))
  expect_equal(ss$ids[ss$genotyped_rows], rownames(ss$markers))
  expect_match(
    spec$bridge$target,
    "fit_metafounder_single_step_reml",
    fixed = TRUE
  )

  payload <- hsquared:::hs_build_bridge_payload(spec)
  expect_equal(payload$relationship_source, "metafounder_single_step")
  expect_equal(payload$relationship, "metafounder_single_step")
  expect_equal(payload$group_of, ss$group_of)
  expect_equal(payload$Gamma, ss$Gamma)
  expect_equal(payload$gamma_labels, c("base_s", "base_d"))
  expect_equal(payload$metadata$gamma_source, "supplied")
  expect_equal(payload$metadata$n_genotyped, 3L)
  expect_match(
    payload$metadata$julia_fit_target,
    "fit_metafounder_single_step_reml",
    fixed = TRUE
  )
})

test_that("single_step H^Gamma payload rejects malformed group and Gamma", {
  ped <- data.frame(
    id = c("c1", "s", "d"),
    sire = c("s", NA, NA),
    dam = c("d", NA, NA),
    stringsAsFactors = FALSE
  )
  dat <- data.frame(y = c(1, 2, 3), id = c("c1", "s", "d"))
  m <- matrix(c(0, 1), nrow = 1, dimnames = list("c1", c("m1", "m2")))
  mf_group <- c(c1 = "", s = "base_s", d = "base_d")
  Gamma <- diag(2)
  dimnames(Gamma) <- list(c("base_s", "base_d"), c("base_s", "base_d"))

  expect_error(
    hsquared:::hs_build_model_spec(
      y ~ single_step(
        1 | id,
        pedigree = ped,
        markers = m,
        group = mf_group
      ),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "requires both `group` and `Gamma`",
    fixed = TRUE
  )
  expect_error(
    hsquared:::hs_build_model_spec(
      y ~ single_step(
        1 | id,
        pedigree = ped,
        markers = m,
        group = c(c1 = "", s = "base_s"),
        Gamma = Gamma
      ),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "missing pedigree id",
    fixed = TRUE
  )
  expect_error(
    hsquared:::hs_build_model_spec(
      y ~ single_step(
        1 | id,
        pedigree = ped,
        markers = m,
        group = c(c1 = "", s = "base_s", d = ""),
        Gamma = Gamma
      ),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "non-missing labels",
    fixed = TRUE
  )
  bad_gamma <- matrix(c(1, 0, 0, 1, 0, 0), nrow = 2)
  expect_error(
    hsquared:::hs_build_model_spec(
      y ~ single_step(
        1 | id,
        pedigree = ped,
        markers = m,
        group = mf_group,
        Gamma = bad_gamma
      ),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "Gamma` must be square",
    fixed = TRUE
  )
  bad_names <- Gamma
  dimnames(bad_names) <- list(c("base_s", "ghost"), c("base_s", "ghost"))
  expect_error(
    hsquared:::hs_build_model_spec(
      y ~ single_step(
        1 | id,
        pedigree = ped,
        markers = m,
        group = mf_group,
        Gamma = bad_names
      ),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "dimnames must match",
    fixed = TRUE
  )
  bad_psd <- matrix(c(1, 2, 2, 1), nrow = 2)
  expect_error(
    hsquared:::hs_build_model_spec(
      y ~ single_step(
        1 | id,
        pedigree = ped,
        markers = m,
        group = mf_group,
        Gamma = bad_psd
      ),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "positive semidefinite",
    fixed = TRUE
  )
})

test_that("single_step construction rejects malformed inputs", {
  ped <- data.frame(
    id = c("c1", "s", "d"),
    sire = c("s", NA, NA),
    dam = c("d", NA, NA),
    stringsAsFactors = FALSE
  )
  dat <- data.frame(y = c(1, 2, 3), id = c("c1", "s", "d"))
  m_ok <- matrix(c(0, 1), nrow = 1, dimnames = list("c1", c("m1", "m2")))
  m_bad <- matrix(
    c(0, 1, 2, 1),
    nrow = 2,
    dimnames = list(c("c1", "ghost"), c("m1", "m2"))
  )
  # a genotyped id absent from the pedigree
  expect_error(
    hsquared:::hs_build_model_spec(
      y ~ single_step(1 | id, pedigree = ped, markers = m_bad),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "not in the pedigree",
    fixed = TRUE
  )
  # pedigree without markers
  expect_error(
    hsquared:::hs_build_model_spec(
      y ~ single_step(1 | id, pedigree = ped),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "requires",
    fixed = TRUE
  )
  # markers without a pedigree -> directing error (not "planned, not implemented")
  expect_error(
    hsquared:::hs_build_model_spec(
      y ~ single_step(1 | id, markers = m_ok),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "needs a pedigree",
    fixed = TRUE
  )
  # both Hinv and pedigree+markers -> "choose one"
  hinv <- diag(3)
  dimnames(hinv) <- list(c("c1", "s", "d"), c("c1", "s", "d"))
  expect_error(
    hsquared:::hs_build_model_spec(
      y ~ single_step(1 | id, pedigree = ped, markers = m_ok, Hinv = hinv),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "not both",
    fixed = TRUE
  )
})

# --- hs_data() bundle shorthand: single_step(1 | id) resolves ped + markers ---

test_that("single_step(1 | id) resolves pedigree + markers from an hs_data bundle", {
  ped <- data.frame(
    id = c("c1", "s", "d", "c2", "c3"),
    sire = c("s", NA, NA, "s", "s"),
    dam = c("d", NA, NA, "d", "d"),
    stringsAsFactors = FALSE
  )
  dat <- data.frame(y = c(1.1, 2.2, 3.0, 2.5, 1.8), id = ped$id)
  m <- matrix(
    c(0, 1, 2, 2, 1, 0, 1, 1, 1, 0, 2, 1),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("c3", "c1", "c2"), paste0("m", 1:4))
  )
  bundle <- hs_data(phenotypes = dat, pedigree = ped, genotypes = m)

  spec_bundle <- hsquared:::hs_build_model_spec(
    y ~ single_step(1 | id),
    data = bundle,
    family = stats::gaussian(),
    REML = TRUE
  )
  spec_explicit <- hsquared:::hs_build_model_spec(
    y ~ single_step(1 | id, pedigree = ped, markers = m),
    data = dat,
    family = stats::gaussian(),
    REML = TRUE
  )
  ssb <- spec_bundle$random$single_step
  sse <- spec_explicit$random$single_step
  # the shorthand parses to the identical construction spec as the explicit call
  expect_equal(ssb$source, "construct")
  expect_equal(ssb$ids, sse$ids)
  expect_equal(ssb$genotyped_rows, sse$genotyped_rows)
  expect_equal(ssb$markers, sse$markers)
})

test_that("explicit single_step args override the hs_data bundle", {
  ped <- data.frame(
    id = c("c1", "s", "d", "c2", "c3"),
    sire = c("s", NA, NA, "s", "s"),
    dam = c("d", NA, NA, "d", "d"),
    stringsAsFactors = FALSE
  )
  dat <- data.frame(y = c(1.1, 2.2, 3.0, 2.5, 1.8), id = ped$id)
  m_bundle <- matrix(
    c(0, 1, 2, 2, 1, 0, 1, 1, 1, 0, 2, 1),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("c3", "c1", "c2"), paste0("m", 1:4))
  )
  # a different, smaller genotyped set passed explicitly
  m_explicit <- matrix(
    c(0, 1, 2, 2),
    nrow = 1,
    dimnames = list("c1", paste0("m", 1:4))
  )
  bundle <- hs_data(phenotypes = dat, pedigree = ped, genotypes = m_bundle)

  spec <- hsquared:::hs_build_model_spec(
    y ~ single_step(1 | id, markers = m_explicit),
    data = bundle,
    family = stats::gaussian(),
    REML = TRUE
  )
  ss <- spec$random$single_step
  # the explicit markers win: only c1 is genotyped (one row), not the bundle's 3
  expect_equal(rownames(ss$markers), "c1")
  expect_equal(length(ss$genotyped_rows), 1L)
})

test_that("hs_single_step_bundle_markers coerces genotype representations", {
  m <- matrix(
    c(0, 1, 2, 1),
    nrow = 2,
    dimnames = list(c("a", "b"), c("m1", "m2"))
  )
  expect_identical(hsquared:::hs_single_step_bundle_markers(m, "id"), m)

  df_id <- data.frame(
    id = c("a", "b"),
    m1 = c(0, 2),
    m2 = c(1, 1),
    stringsAsFactors = FALSE
  )
  out <- hsquared:::hs_single_step_bundle_markers(df_id, "id")
  expect_equal(rownames(out), c("a", "b"))
  expect_equal(colnames(out), c("m1", "m2"))
  expect_true(is.numeric(out))

  df_rn <- data.frame(m1 = c(0, 2), m2 = c(1, 1))
  rownames(df_rn) <- c("a", "b")
  expect_equal(
    rownames(hsquared:::hs_single_step_bundle_markers(df_rn, "id")),
    c("a", "b")
  )

  # a data frame with neither an id column nor explicit row names is rejected
  df_bad <- data.frame(m1 = c(0, 2), m2 = c(1, 1))
  expect_error(
    hsquared:::hs_single_step_bundle_markers(df_bad, "id"),
    "row names"
  )
})

test_that("the bundle shorthand threads a non-default id + data-frame genotypes", {
  # end-to-end exercise of the new mechanism: id = data$id is threaded into the
  # model-data context and used to coerce a DATA-FRAME genotype bundle (not just
  # the matrix fast path) under a non-default id column ("animal").
  ped <- data.frame(
    id = c("c1", "s", "d", "c2", "c3"),
    sire = c("s", NA, NA, "s", "s"),
    dam = c("d", NA, NA, "d", "d"),
    stringsAsFactors = FALSE
  )
  dat <- data.frame(y = c(1.1, 2.2, 3.0, 2.5, 1.8), animal = ped$id)
  geno_df <- data.frame(
    animal = c("c3", "c1", "c2"),
    m1 = c(0, 1, 1),
    m2 = c(1, 1, 0),
    m3 = c(2, 0, 2),
    m4 = c(2, 1, 1),
    stringsAsFactors = FALSE
  )
  bundle <- hs_data(
    phenotypes = dat,
    pedigree = ped,
    genotypes = geno_df,
    id = "animal"
  )
  spec <- hsquared:::hs_build_model_spec(
    y ~ single_step(1 | animal),
    data = bundle,
    family = stats::gaussian(),
    REML = TRUE
  )
  ss <- spec$random$single_step
  expect_equal(ss$source, "construct")
  expect_true(is.numeric(ss$markers))
  expect_equal(sort(rownames(ss$markers)), c("c1", "c2", "c3"))
  # genotyped_rows index exactly the genotyped animals' pedigree rows
  expect_equal(ss$ids[ss$genotyped_rows], rownames(ss$markers))
})

test_that("single_step construction errors direct to the on-ramps", {
  dat <- data.frame(y = c(1, 2, 3), id = c("c1", "s", "d"))
  # bare single_step(1 | id) on a plain data frame: the Hinv error now also points
  # at the construction on-ramps (explicit args or an hs_data() bundle).
  expect_error(
    hsquared:::hs_build_model_spec(
      y ~ single_step(1 | id),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "construct one"
  )
  # a partial bundle (pedigree only, no genotypes) also directs to construction
  ped <- data.frame(
    id = c("c1", "s", "d"),
    sire = c("s", NA, NA),
    dam = c("d", NA, NA),
    stringsAsFactors = FALSE
  )
  bundle <- hs_data(phenotypes = dat, pedigree = ped)
  expect_error(
    hsquared:::hs_build_model_spec(
      y ~ single_step(1 | id),
      data = bundle,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "construct one"
  )
})

test_that("metafounder_single_step target requires the H^Gamma formula branch", {
  expect_equal(
    hsquared:::hs_validate_julia_target("metafounder_single_step"),
    "metafounder_single_step"
  )

  ped <- data.frame(
    id = c("c1", "s", "d"),
    sire = c("s", NA, NA),
    dam = c("d", NA, NA),
    stringsAsFactors = FALSE
  )
  dat <- data.frame(y = c(1, 2, 3), id = c("c1", "s", "d"))
  m <- matrix(c(0, 1), nrow = 1, dimnames = list("c1", c("m1", "m2")))
  mf_group <- c(c1 = "", s = "base_s", d = "base_d")
  Gamma <- diag(2)
  dimnames(Gamma) <- list(c("base_s", "base_d"), c("base_s", "base_d"))

  expect_error(
    hsquared(
      y ~ single_step(1 | id, pedigree = ped, markers = m),
      data = dat,
      control = hs_control(
        engine = "julia",
        engine_control = list(target = "metafounder_single_step")
      )
    ),
    "requires `single_step",
    fixed = TRUE
  )

  expect_error(
    hsquared(
      y ~ single_step(
        1 | id,
        pedigree = ped,
        markers = m,
        group = mf_group,
        Gamma = Gamma
      ),
      data = dat,
      control = hs_control(
        engine = "julia",
        engine_control = list(target = "single_step_construct")
      )
    ),
    "not the ordinary",
    fixed = TRUE
  )
})

# --- live tests (skip-guarded on the local HSquared.jl bridge) ---------------

hs_ss_construct_fit <- function(
  ped,
  dat,
  markers,
  ridge = 0.01,
  iterations = 100L
) {
  hsquared(
    y ~ single_step(1 | id, pedigree = ped, markers = markers, ridge = ridge),
    data = dat,
    control = hs_control(
      engine = "julia",
      engine_control = list(
        target = "single_step_construct",
        iterations = iterations
      )
    )
  )
}

hs_mf_group <- function(ped, label = "base") {
  out <- rep("", nrow(ped))
  out[is.na(ped$sire) | is.na(ped$dam)] <- label
  names(out) <- as.character(ped$id)
  out
}

hs_mf_single_step_fit <- function(
  ped,
  dat,
  markers,
  group,
  Gamma,
  ridge = 0.01,
  iterations = 100L
) {
  hsquared(
    y ~ single_step(
      1 | id,
      pedigree = ped,
      markers = markers,
      group = group,
      Gamma = Gamma,
      ridge = ridge
    ),
    data = dat,
    control = hs_control(
      engine = "julia",
      engine_control = list(
        target = "metafounder_single_step",
        iterations = iterations
      )
    )
  )
}

test_that("single_step construction is invariant to marker row order [live]", {
  # THE alignment guard (docs/design/25 §6.3): shuffling the marker rows (same
  # underlying genotypes) must give an identical fit -- a missing/wrong reorder
  # would place G at the wrong H^-1 rows and change the fit.
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for the live single-step test."
  )
  hsquared:::hs_julia_setup(hsquared:::hs_default_julia_project())

  ped <- hs_sim_pedigree(n_founder = 10, n_per_gen = 20, n_gen = 2, seed = 21)
  dat <- hs_sim_genedrop_phenotypes(
    ped,
    sigma_a2 = 0.4,
    sigma_e2 = 0.6,
    seed = 21
  )
  geno <- utils::tail(ped$id, 26L)
  set.seed(21)
  markers <- matrix(
    stats::rbinom(length(geno) * 60L, 2L, 0.3),
    nrow = length(geno),
    dimnames = list(geno, paste0("snp", seq_len(60L)))
  )
  perm <- c(seq(2L, length(geno), by = 2L), seq(1L, length(geno), by = 2L))
  markers_shuf <- markers[perm, , drop = FALSE]

  fit_a <- hs_ss_construct_fit(ped, dat, markers)
  fit_b <- hs_ss_construct_fit(ped, dat, markers_shuf)

  expect_equal(
    variance_components(fit_a)$estimate,
    variance_components(fit_b)$estimate,
    tolerance = 1e-8
  )
  ba <- breeding_values(fit_a)
  bb <- breeding_values(fit_b)
  expect_equal(ba$value[order(ba$id)], bb$value[order(bb$id)], tolerance = 1e-8)
})

test_that("metafounder single-step Gamma = 0 reduces to ordinary construction [live]", {
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for the live metafounder single-step test."
  )
  hsquared:::hs_julia_setup(hsquared:::hs_default_julia_project())

  ped <- hs_sim_pedigree(n_founder = 10, n_per_gen = 20, n_gen = 2, seed = 25)
  dat <- hs_sim_genedrop_phenotypes(
    ped,
    sigma_a2 = 0.4,
    sigma_e2 = 0.6,
    seed = 25
  )
  geno <- utils::tail(ped$id, 30L)
  set.seed(25)
  markers <- matrix(
    stats::rbinom(length(geno) * 60L, 2L, 0.3),
    nrow = length(geno),
    dimnames = list(geno, paste0("snp", seq_len(60L)))
  )
  group <- hs_mf_group(ped)
  Gamma0 <- matrix(0, nrow = 1, dimnames = list("base", "base"))

  ordinary <- hs_ss_construct_fit(
    ped,
    dat,
    markers,
    ridge = 0.05,
    iterations = 500L
  )
  hgamma0 <- hs_mf_single_step_fit(
    ped,
    dat,
    markers,
    group = group,
    Gamma = Gamma0,
    ridge = 0.05,
    iterations = 500L
  )

  expect_equal(hgamma0$spec$target, "metafounder_single_step")
  expect_equal(
    fit_diagnostics(hgamma0)$value[
      fit_diagnostics(hgamma0)$metric == "gamma_source"
    ],
    "supplied"
  )
  expect_equal(gamma_matrix(hgamma0), Gamma0)
  expect_equal(
    variance_components(hgamma0)$estimate,
    variance_components(ordinary)$estimate,
    tolerance = 1e-6
  )
  expect_equal(hgamma0$result$converged, ordinary$result$converged)
  expect_equal(hgamma0$result$loglik, ordinary$result$loglik, tolerance = 1e-6)
  bv_hg <- breeding_values(hgamma0)
  bv_ord <- breeding_values(ordinary)
  expect_equal(bv_hg$id[order(bv_hg$id)], bv_ord$id[order(bv_ord$id)])
  expect_equal(
    bv_hg$value[order(bv_hg$id)],
    bv_ord$value[order(bv_ord$id)],
    tolerance = 1e-6
  )
})

test_that("metafounder single-step nonzero Gamma changes predictions [live]", {
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for the live metafounder single-step test."
  )
  hsquared:::hs_julia_setup(hsquared:::hs_default_julia_project())

  ped <- hs_sim_pedigree(n_founder = 10, n_per_gen = 20, n_gen = 2, seed = 26)
  dat <- hs_sim_genedrop_phenotypes(
    ped,
    sigma_a2 = 0.4,
    sigma_e2 = 0.6,
    seed = 26
  )
  geno <- utils::tail(ped$id, 30L)
  set.seed(26)
  markers <- matrix(
    stats::rbinom(length(geno) * 60L, 2L, 0.3),
    nrow = length(geno),
    dimnames = list(geno, paste0("snp", seq_len(60L)))
  )
  group <- hs_mf_group(ped)
  Gamma0 <- matrix(0, nrow = 1, dimnames = list("base", "base"))
  Gamma <- matrix(0.25, nrow = 1, dimnames = list("base", "base"))

  hgamma0 <- hs_mf_single_step_fit(
    ped,
    dat,
    markers,
    group = group,
    Gamma = Gamma0,
    ridge = 0.05,
    iterations = 500L
  )
  hgamma <- hs_mf_single_step_fit(
    ped,
    dat,
    markers,
    group = group,
    Gamma = Gamma,
    ridge = 0.05,
    iterations = 500L
  )

  bv0 <- breeding_values(hgamma0)
  bvg <- breeding_values(hgamma)
  diag <- fit_diagnostics(hgamma)
  pev <- prediction_error_variance(hgamma)
  rel <- reliability(hgamma)
  acc <- accuracy(hgamma)

  expect_equal(gamma_matrix(hgamma0), Gamma0)
  expect_equal(gamma_matrix(hgamma), Gamma)
  expect_equal(
    metafounder_groups(hgamma),
    data.frame(
      id = as.character(ped$id),
      metafounder_group = ifelse(nzchar(group[ped$id]), group[ped$id], NA),
      is_metafounder = nzchar(group[ped$id]),
      stringsAsFactors = FALSE
    )
  )
  expect_setequal(bvg$id, as.character(ped$id))
  expect_equal(nrow(bvg), nrow(ped))
  expect_true(all(is.finite(bvg$value)))
  expect_equal(
    diag$value[diag$metric == "target"],
    "metafounder_single_step"
  )
  expect_equal(
    diag$value[diag$metric == "variance_components_source"],
    "estimated_metafounder_single_step_ai_reml"
  )
  expect_equal(diag$value[diag$metric == "gamma_source"], "supplied")
  expect_equal(pev$id, as.character(ped$id))
  expect_equal(rel$id, as.character(ped$id))
  expect_equal(acc$id, as.character(ped$id))
  expect_true(all(is.finite(pev$value)))
  expect_true(all(rel$value >= 0 & rel$value <= 1))
  expect_equal(acc$value, sqrt(rel$value), tolerance = 1e-12)
  merged <- merge(bv0, bvg, by = "id", suffixes = c("_zero", "_gamma"))
  expect_equal(nrow(merged), nrow(ped))
  expect_gt(mean(abs(merged$value_zero - merged$value_gamma)), 1e-6)
})

test_that("single_step construction labels + covers all pedigree animals [live]", {
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for the live single-step test."
  )
  hsquared:::hs_julia_setup(hsquared:::hs_default_julia_project())

  ped <- hs_sim_pedigree(n_founder = 10, n_per_gen = 20, n_gen = 2, seed = 22)
  dat <- hs_sim_genedrop_phenotypes(
    ped,
    sigma_a2 = 0.4,
    sigma_e2 = 0.6,
    seed = 22
  )
  geno <- utils::tail(ped$id, 26L)
  set.seed(22)
  markers <- matrix(
    stats::rbinom(length(geno) * 60L, 2L, 0.3),
    nrow = length(geno),
    dimnames = list(geno, paste0("snp", seq_len(60L)))
  )
  fit_c <- hs_ss_construct_fit(ped, dat, markers)
  bv_c <- breeding_values(fit_c)

  # GEBVs are labelled by the real pedigree ids (not positional integers) and
  # cover ALL pedigree animals, genotyped and ungenotyped.
  expect_setequal(bv_c$id, as.character(ped$id))
  expect_equal(nrow(bv_c), nrow(ped))
  expect_true(all(is.finite(bv_c$value)))
  ungeno <- setdiff(as.character(ped$id), geno)
  expect_true(length(ungeno) > 0L)
  expect_true(all(is.finite(bv_c$value[bv_c$id %in% ungeno])))
})

test_that("single_step construction uses the genomic info (differs from the pedigree model) [live]", {
  # Independent correctness anchor: the single-step GEBVs must TRACK the pedigree
  # animal-model GEBVs (same signal) but DIFFER from them (the genomic info is
  # actually used) -- not a self-referential check.
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for the live single-step test."
  )
  hsquared:::hs_julia_setup(hsquared:::hs_default_julia_project())

  ped <- hs_sim_pedigree(n_founder = 10, n_per_gen = 20, n_gen = 2, seed = 23)
  dat <- hs_sim_genedrop_phenotypes(
    ped,
    sigma_a2 = 0.4,
    sigma_e2 = 0.6,
    seed = 23
  )
  geno <- utils::tail(ped$id, 26L)
  set.seed(23)
  markers <- matrix(
    stats::rbinom(length(geno) * 60L, 2L, 0.3),
    nrow = length(geno),
    dimnames = list(geno, paste0("snp", seq_len(60L)))
  )
  fit_c <- hs_ss_construct_fit(ped, dat, markers)
  fit_ped <- hsquared(
    y ~ animal(1 | id, pedigree = ped),
    data = dat,
    family = stats::gaussian(),
    REML = TRUE
  )
  m <- merge(
    breeding_values(fit_c),
    breeding_values(fit_ped),
    by = "id",
    suffixes = c("_ss", "_ped")
  )
  expect_equal(nrow(m), nrow(ped))
  # tracks the same signal ...
  expect_gt(stats::cor(m$value_ss, m$value_ped), 0.5)
  # ... but the genomic information genuinely changes the predictions.
  expect_gt(mean(abs(m$value_ss - m$value_ped)), 1e-6)
})

test_that("single_step construction ridge fits a rank-deficient G [live]", {
  # singular-G path (docs/design/25 §6.4): more genotyped animals than markers ->
  # raw G is rank-deficient; a positive ridge makes H^-1 PD and the fit succeed.
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for the live single-step test."
  )
  hsquared:::hs_julia_setup(hsquared:::hs_default_julia_project())

  ped <- hs_sim_pedigree(n_founder = 10, n_per_gen = 20, n_gen = 2, seed = 24)
  dat <- hs_sim_genedrop_phenotypes(
    ped,
    sigma_a2 = 0.4,
    sigma_e2 = 0.6,
    seed = 24
  )
  geno <- utils::tail(ped$id, 30L)
  set.seed(24)
  markers <- matrix(
    stats::rbinom(length(geno) * 8L, 2L, 0.3), # 8 markers < 30 genotyped -> singular
    nrow = length(geno),
    dimnames = list(geno, paste0("snp", seq_len(8L)))
  )
  fit_c <- hs_ss_construct_fit(ped, dat, markers, ridge = 0.05)
  expect_true(fit_c$result$converged)
  expect_true(all(variance_components(fit_c)$estimate > 0))
  expect_true(
    heritability(fit_c)$estimate > 0 && heritability(fit_c)$estimate < 1
  )
})

test_that("the single_step(1 | id) bundle shorthand fits like the explicit call [live]", {
  hs_skip_live_julia()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for the live single-step test."
  )
  hsquared:::hs_julia_setup(hsquared:::hs_default_julia_project())

  # mirror the proven rank-deficient setup (8 markers < genotyped, ridge = 0.05)
  # so the genotyped genomic block is positive definite for both fits.
  ped <- hs_sim_pedigree(n_founder = 10, n_per_gen = 20, n_gen = 2, seed = 24)
  dat <- hs_sim_genedrop_phenotypes(
    ped,
    sigma_a2 = 0.4,
    sigma_e2 = 0.6,
    seed = 24
  )
  geno <- utils::tail(ped$id, 30L)
  set.seed(24)
  markers <- matrix(
    stats::rbinom(length(geno) * 8L, 2L, 0.3),
    nrow = length(geno),
    dimnames = list(geno, paste0("snp", seq_len(8L)))
  )

  # explicit call vs the hs_data() bundle shorthand: identical fit
  fit_explicit <- hs_ss_construct_fit(ped, dat, markers, ridge = 0.05)
  bundle <- hs_data(phenotypes = dat, pedigree = ped, genotypes = markers)
  fit_bundle <- hsquared(
    y ~ single_step(1 | id, ridge = 0.05),
    data = bundle,
    control = hs_control(
      engine = "julia",
      engine_control = list(target = "single_step_construct")
    )
  )

  expect_equal(
    variance_components(fit_bundle)$estimate,
    variance_components(fit_explicit)$estimate,
    tolerance = 1e-8
  )
  be <- breeding_values(fit_explicit)
  bb <- breeding_values(fit_bundle)
  expect_equal(bb$id[order(bb$id)], be$id[order(be$id)])
  expect_equal(
    bb$value[order(bb$id)],
    be$value[order(be$id)],
    tolerance = 1e-8
  )
})
