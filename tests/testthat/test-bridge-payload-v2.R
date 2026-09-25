# Tests for the payload-v2 random_effects block structure added by P0.4.
# Spec: docs/design/21-payload-v2-multiblock-schema.md (frozen schema).
# These tests ONLY exercise hs_build_bridge_payload() (the emitter). They do NOT
# call the Julia bridge. All legacy v0.1 top-level fields must be byte-identical
# to what the function produced before P0.4.

# ---- shared test fixtures -------------------------------------------------- #

make_ped <- function() {
  data.frame(
    id   = c("a", "b", "c", "d"),
    sire = c(NA,  NA,  "a", "a"),
    dam  = c(NA,  NA,  "b", "c"),
    stringsAsFactors = FALSE
  )
}

make_dat_single <- function() {
  data.frame(
    y   = c(1, 2, 3),
    sex = c("f", "m", "f"),
    id  = c("a", "c", "d"),
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------- #
# Case (a): v0.1 single-pedigree animal model
# ---------------------------------------------------------------------------- #

test_that("payload-v2 single animal model: payload_version = 2L", {
  ped <- make_ped()
  dat <- make_dat_single()
  spec <- hsquared:::hs_build_model_spec(
    y ~ sex + animal(1 | id, pedigree = ped),
    data = dat,
    family = stats::gaussian(),
    REML = TRUE
  )
  payload <- hsquared:::hs_build_bridge_payload(spec)
  expect_equal(payload$payload_version, 2L)
})

test_that("payload-v2 single animal model: exactly one pedigree block", {
  ped <- make_ped()
  dat <- make_dat_single()
  spec <- hsquared:::hs_build_model_spec(
    y ~ sex + animal(1 | id, pedigree = ped),
    data = dat,
    family = stats::gaussian(),
    REML = TRUE
  )
  payload <- hsquared:::hs_build_bridge_payload(spec)
  re <- payload$random_effects

  expect_length(re, 1L)
  expect_equal(re[[1L]]$name, "animal")
  expect_equal(re[[1L]]$type, "pedigree")
  expect_equal(re[[1L]]$relmat_status, "build_in_julia")
  expect_null(re[[1L]]$relmat_inverse)
  # pedigree sub-list must carry the same fields as the top-level pedigree
  expect_equal(re[[1L]]$pedigree$id,             payload$pedigree$id)
  expect_equal(re[[1L]]$pedigree$sire,           payload$pedigree$sire)
  expect_equal(re[[1L]]$pedigree$dam,            payload$pedigree$dam)
  expect_equal(re[[1L]]$pedigree$sire_index,     payload$pedigree$sire_index)
  expect_equal(re[[1L]]$pedigree$dam_index,      payload$pedigree$dam_index)
  expect_equal(re[[1L]]$pedigree$original_order, payload$pedigree$original_order)
  # ids must equal the top-level ids
  expect_equal(re[[1L]]$ids, payload$ids)
  # Z must be the same sparse matrix as the top-level Z
  expect_equal(as.matrix(re[[1L]]$Z), as.matrix(payload$Z))
})

test_that("payload-v2 single animal model: all legacy v0.1 fields unchanged", {
  # Snapshot the key fields that existed before P0.4 and verify they are
  # byte-identical now that payload_version + random_effects have been added.
  ped <- make_ped()
  dat <- make_dat_single()
  spec <- hsquared:::hs_build_model_spec(
    y ~ sex + animal(1 | id, pedigree = ped),
    data = dat,
    family = stats::gaussian(),
    REML = TRUE
  )
  payload <- hsquared:::hs_build_bridge_payload(spec)

  # class
  expect_s3_class(payload, "hs_bridge_payload")

  # response
  expect_equal(payload$y, dat$y)
  expect_null(payload$Y)

  # design matrix (intercept + sexm = 2 columns)
  expect_equal(dim(payload$X), c(3L, 2L))
  expect_null(colnames(payload$X))
  expect_equal(payload$metadata$fixed_colnames, c("(Intercept)", "sexm"))

  # primary incidence
  expect_s4_class(payload$Z, "dgCMatrix")
  expect_equal(dim(payload$Z), c(3L, 4L))
  expect_equal(colnames(payload$Z), c("a", "b", "c", "d"))
  expect_equal(
    unname(as.matrix(payload$Z)),
    rbind(c(1, 0, 0, 0), c(0, 0, 1, 0), c(0, 0, 0, 1))
  )

  # no second effect
  expect_null(payload$Z2)
  expect_null(payload$effect2)

  # Ainv always NULL
  expect_null(payload$Ainv)

  # method / family
  expect_equal(payload$method, "REML")
  expect_equal(payload$family, "gaussian")

  # ids + pedigree
  expect_equal(payload$ids, c("a", "b", "c", "d"))
  expect_equal(payload$pedigree$sire_index, c(0L, 0L, 1L, 1L))
  expect_equal(payload$pedigree$dam_index,  c(0L, 0L, 2L, 3L))

  # metadata
  expect_equal(payload$metadata$observed_ids,       c("a", "c", "d"))
  expect_equal(payload$metadata$observed_id_index,  c(1L, 3L, 4L))
  expect_equal(payload$metadata$ainv_status,         "build_in_julia")
  expect_match(payload$metadata$julia_spec_target,  "animal_model_spec")
})

# ---------------------------------------------------------------------------- #
# Case (b): animal + common_env()  → pedigree block + iid block
# ---------------------------------------------------------------------------- #

test_that("payload-v2 animal + common_env: two blocks (pedigree + iid)", {
  ped <- make_ped()
  dat <- data.frame(
    y      = c(1, 2, 3, 4),
    id     = c("a", "c", "d", "b"),
    litter = c("L1", "L1", "L2", "L2"),
    stringsAsFactors = FALSE
  )
  spec <- hsquared:::hs_build_model_spec(
    y ~ animal(1 | id, pedigree = ped) + common_env(1 | litter),
    data = dat,
    family = stats::gaussian(),
    REML = TRUE
  )
  payload <- hsquared:::hs_build_bridge_payload(spec)
  re <- payload$random_effects

  expect_length(re, 2L)

  # block 1
  expect_equal(re[[1L]]$name,          "animal")
  expect_equal(re[[1L]]$type,          "pedigree")
  expect_equal(re[[1L]]$relmat_status, "build_in_julia")

  # block 2
  expect_equal(re[[2L]]$name,          "common_env")
  expect_equal(re[[2L]]$type,          "iid")
  expect_equal(re[[2L]]$relmat_status, "identity")
  expect_null(re[[2L]]$relmat_inverse)
  expect_null(re[[2L]]$pedigree)

  # Z2 incidence dims: 4 records × 2 litter levels
  expect_s4_class(re[[2L]]$Z, "dgCMatrix")
  expect_equal(dim(re[[2L]]$Z), c(4L, 2L))
  # ids are the litter levels
  expect_equal(sort(re[[2L]]$ids), sort(c("L1", "L2")))

  # legacy fields preserved
  expect_s3_class(payload, "hs_bridge_payload")
  expect_equal(payload$payload_version, 2L)
  expect_s4_class(payload$Z2, "dgCMatrix")
  expect_equal(
    payload$effect2$type,         "common_env"
  )
  expect_equal(payload$effect2$relationship, "identity")
})

# ---------------------------------------------------------------------------- #
# Case (c): animal + permanent()  → pedigree block + iid block
# ---------------------------------------------------------------------------- #

test_that("payload-v2 animal + permanent: two blocks (pedigree + iid)", {
  ped <- make_ped()
  # repeatability: repeated records per individual
  dat <- data.frame(
    y  = c(1, 2, 3, 4, 5, 6),
    id = c("a", "a", "b", "b", "c", "c"),
    stringsAsFactors = FALSE
  )
  spec <- hsquared:::hs_build_model_spec(
    y ~ animal(1 | id, pedigree = ped) + permanent(1 | id),
    data = dat,
    family = stats::gaussian(),
    REML = TRUE
  )
  payload <- hsquared:::hs_build_bridge_payload(spec)
  re <- payload$random_effects

  expect_length(re, 2L)

  # block 1 — animal pedigree
  expect_equal(re[[1L]]$name,          "animal")
  expect_equal(re[[1L]]$type,          "pedigree")
  expect_equal(re[[1L]]$relmat_status, "build_in_julia")

  # block 2 — permanent environment (iid, shares animal Z)
  expect_equal(re[[2L]]$name,          "permanent")
  expect_equal(re[[2L]]$type,          "iid")
  expect_equal(re[[2L]]$relmat_status, "identity")
  expect_null(re[[2L]]$relmat_inverse)
  expect_null(re[[2L]]$pedigree)

  # permanent block shares the same Z as the animal block (same dims)
  expect_s4_class(re[[2L]]$Z, "dgCMatrix")
  expect_equal(dim(re[[2L]]$Z), dim(re[[1L]]$Z))
  # ids are the observed individual ids (unique values from dat$id)
  expect_true(length(re[[2L]]$ids) >= 1L)

  # legacy: Z2 / effect2 are still NULL for permanent (the repeatability bridge
  # uses Z directly — no second incidence matrix)
  expect_null(payload$Z2)
  expect_null(payload$effect2)
  expect_equal(payload$payload_version, 2L)
  expect_s3_class(payload, "hs_bridge_payload")
})

# ---------------------------------------------------------------------------- #
# Case (d): animal + maternal_genetic() → two pedigree blocks
# ---------------------------------------------------------------------------- #

test_that("payload-v2 animal + maternal_genetic: two pedigree blocks", {
  ped <- make_ped()
  dat <- data.frame(
    y   = c(1, 2, 3),
    id  = c("c", "d", "c"),
    dam = c("b", "c", "c"),
    stringsAsFactors = FALSE
  )
  spec <- hsquared:::hs_build_model_spec(
    y ~ animal(1 | id, pedigree = ped) + maternal_genetic(1 | dam),
    data = dat,
    family = stats::gaussian(),
    REML = TRUE
  )
  payload <- hsquared:::hs_build_bridge_payload(spec)
  re <- payload$random_effects

  expect_length(re, 2L)

  # block 1 — direct animal (pedigree)
  expect_equal(re[[1L]]$name,          "animal")
  expect_equal(re[[1L]]$type,          "pedigree")
  expect_equal(re[[1L]]$relmat_status, "build_in_julia")

  # block 2 — maternal genetic (independent pedigree, NOT correlated — schema §3)
  expect_equal(re[[2L]]$name,          "maternal")
  expect_equal(re[[2L]]$type,          "pedigree")
  expect_equal(re[[2L]]$relmat_status, "build_in_julia")
  expect_null(re[[2L]]$relmat_inverse)
  # maternal block carries pedigree rows (shares Ainv construction)
  expect_false(is.null(re[[2L]]$pedigree))
  expect_equal(re[[2L]]$pedigree$id, payload$pedigree$id)

  # Z2 incidence dims: 3 records × q2 pedigree levels
  expect_s4_class(re[[2L]]$Z, "dgCMatrix")
  expect_equal(nrow(re[[2L]]$Z), nrow(dat))

  expect_equal(payload$payload_version, 2L)
})

test_that("payload-v2 preserves the frozen bridge schema", {
  ped <- make_ped()
  dat <- make_dat_single()
  spec <- hsquared:::hs_build_model_spec(
    y ~ sex + animal(1 | id, pedigree = ped),
    data = dat,
    family = stats::gaussian(),
    REML = TRUE
  )
  payload <- hsquared:::hs_build_bridge_payload(spec)

  # These names are the R-to-Julia contract, not an incidental list shape.
  # Keep this Julia-free guard close to the emitter so a new public route
  # cannot silently bypass the declared payload schema.
  expect_setequal(
    names(payload),
    c(
      "payload_version", "y", "Y", "X", "Z", "Z2", "effect2",
      "random_effects", "random_regression", "Ainv", "group_of",
      "Gamma", "gamma_labels", "relationship_source", "method", "family",
      "n_trials", "ids", "pedigree", "metadata"
    )
  )
  expect_equal(payload$payload_version, 2L)
  expect_setequal(
    names(payload$metadata),
    c(
      "response", "response_type", "trait_names", "random_regression",
      "fixed_colnames", "animal_id_column", "observed_ids",
      "observed_id_index", "relationship", "gamma_source", "fixed_terms",
      "contrasts", "ainv_status", "ainv_target", "julia_spec_target",
      "julia_fit_target"
    )
  )
  expect_setequal(
    names(payload$random_effects[[1L]]),
    c("name", "type", "Z", "relmat_inverse", "relmat_status", "pedigree", "ids")
  )
  expect_equal(payload$random_effects[[1L]]$type, "pedigree")
  expect_equal(payload$random_effects[[1L]]$relmat_status, "build_in_julia")
  expect_match(payload$metadata$julia_spec_target, "animal_model_spec", fixed = TRUE)
  expect_match(payload$metadata$julia_fit_target, "fit_animal_model", fixed = TRUE)
})
