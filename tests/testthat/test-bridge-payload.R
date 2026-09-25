test_that("bridge payload exposes the v0.1 Julia contract", {
  ped <- data.frame(
    id = c("a", "b", "c", "d"),
    sire = c(NA, NA, "a", "a"),
    dam = c(NA, NA, "b", "c")
  )
  dat <- data.frame(
    y = c(1, 2, 3),
    sex = c("f", "m", "f"),
    age = c(1, 2, 3),
    id = c("a", "c", "d")
  )

  spec <- hsquared:::hs_build_model_spec(
    y ~ sex + age + animal(1 | id, pedigree = ped),
    data = dat,
    family = stats::gaussian(),
    REML = TRUE
  )
  payload <- hsquared:::hs_build_bridge_payload(spec)

  expect_s3_class(payload, "hs_bridge_payload")
  expect_equal(payload$y, dat$y)
  expect_equal(dim(payload$X), c(nrow(dat), 3L))
  expect_null(colnames(payload$X))
  expect_equal(payload$metadata$fixed_colnames, c("(Intercept)", "sexm", "age"))
  expect_s4_class(payload$Z, "dgCMatrix")
  expect_equal(dim(payload$Z), c(3L, 4L))
  expect_equal(colnames(payload$Z), c("a", "b", "c", "d"))
  expect_equal(
    unname(as.matrix(payload$Z)),
    rbind(
      c(1, 0, 0, 0),
      c(0, 0, 1, 0),
      c(0, 0, 0, 1)
    )
  )
  expect_null(payload$Ainv)
  expect_equal(payload$method, "REML")
  expect_equal(payload$family, "gaussian")
  expect_equal(payload$ids, c("a", "b", "c", "d"))
  expect_equal(payload$pedigree$sire_index, c(0L, 0L, 1L, 1L))
  expect_equal(payload$pedigree$dam_index, c(0L, 0L, 2L, 3L))
  expect_equal(payload$metadata$observed_ids, c("a", "c", "d"))
  expect_equal(payload$metadata$observed_id_index, c(1L, 3L, 4L))
  expect_equal(payload$metadata$ainv_status, "build_in_julia")
  expect_match(payload$metadata$julia_spec_target, "animal_model_spec")
})

test_that("bridge payload normalizes pedigree order for Julia Ainv construction", {
  ped <- data.frame(
    id = c("d", "c", "b", "a"),
    sire = c("a", "a", NA, NA),
    dam = c("c", "b", NA, NA)
  )
  dat <- data.frame(
    y = c(1, 2, 3),
    id = c("d", "a", "c")
  )

  spec <- hsquared:::hs_build_model_spec(
    y ~ animal(1 | id, pedigree = ped),
    data = dat,
    family = stats::gaussian(),
    REML = FALSE
  )
  payload <- hsquared:::hs_build_bridge_payload(spec)

  expect_equal(payload$ids, c("a", "b", "c", "d"))
  expect_equal(payload$pedigree$id, c("a", "b", "c", "d"))
  expect_equal(payload$pedigree$sire_index, c(0L, 0L, 1L, 1L))
  expect_equal(payload$pedigree$dam_index, c(0L, 0L, 2L, 3L))
  expect_equal(payload$metadata$observed_id_index, c(4L, 1L, 3L))
  expect_equal(payload$method, "ML")
})

test_that("bridge payload carries animal-only metafounder fields", {
  ped <- data.frame(
    id = c("calf", "sire", "dam"),
    sire = c("sire", NA, NA),
    dam = c("dam", NA, NA)
  )
  dat <- data.frame(y = c(1, 2.5, 4), id = c("sire", "dam", "calf"))
  mf_group <- c(sire = "base", dam = "base", calf = "")
  Gamma <- matrix(0.25, nrow = 1, dimnames = list("base", "base"))

  spec <- hsquared:::hs_build_model_spec(
    y ~ metafounder(1 | id, pedigree = ped, group = mf_group, Gamma = Gamma),
    data = dat,
    family = stats::gaussian(),
    REML = TRUE
  )
  payload <- hsquared:::hs_build_bridge_payload(spec)

  expect_s3_class(payload, "hs_bridge_payload")
  expect_equal(payload$relationship_source, "metafounder")
  expect_equal(payload$ids, c("sire", "dam", "calf"))
  expect_equal(unname(payload$group_of), c("base", "base", ""))
  expect_equal(payload$Gamma, matrix(0.25, nrow = 1))
  expect_equal(payload$gamma_labels, "base")
  expect_equal(payload$metadata$gamma_source, "supplied")
  expect_match(payload$metadata$julia_spec_target, "metafounder_animal_model")
  expect_equal(
    unname(as.matrix(payload$Z)),
    rbind(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1))
  )
})

test_that("pedigree cycles are rejected before bridge payload construction", {
  ped <- data.frame(
    id = c("a", "b"),
    sire = c("b", "a"),
    dam = c(NA, NA)
  )
  dat <- data.frame(y = c(1, 2), id = c("a", "b"))

  expect_error(
    hsquared:::hs_build_model_spec(
      y ~ animal(1 | id, pedigree = ped),
      data = dat,
      family = stats::gaussian(),
      REML = TRUE
    ),
    "parent-offspring cycle",
    fixed = TRUE
  )
})

test_that("hsquared validate-only message describes the Julia fit target", {
  ped <- data.frame(
    id = c("a", "b"),
    sire = c(NA, NA),
    dam = c(NA, NA)
  )
  dat <- data.frame(y = c(1, 2), id = c("a", "b"))

  expect_message(
    spec <- hsquared(
      y ~ animal(1 | id, pedigree = ped),
      data = dat,
      family = stats::gaussian(),
      control = hs_control(engine = "validate")
    ),
    "HSquared.fit_animal_model",
    fixed = TRUE
  )
  # The validate path also returns the validated spec (a list) invisibly.
  expect_type(spec, "list")
  expect_match(spec$bridge$target, "fit_animal_model", fixed = TRUE)
})
