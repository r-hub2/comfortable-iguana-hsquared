test_that("Mrode9 validation fixture exposes a nadiv sparse Ainv comparator", {
  testthat::skip_if_not_installed("nadiv")

  fixture <- hsquared:::hs_mrode9_pedigree_validation_fixture()

  expect_equal(fixture$name, "mrode9_nadiv_pedigree")
  expect_equal(nrow(fixture$pedigree), 12L)
  expect_equal(names(fixture$pedigree), c("id", "sire", "dam"))
  expect_s4_class(fixture$expected$Ainv, "dgCMatrix")
  expect_equal(dim(fixture$expected$Ainv), c(12L, 12L))
  expect_equal(rownames(fixture$expected$Ainv), fixture$pedigree$id)
  expect_equal(colnames(fixture$expected$Ainv), fixture$pedigree$id)
})

test_that("Mrode9 nadiv Ainv matches Julia pedigree_inverse when available", {
  testthat::skip_if_not_installed("nadiv")
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(),
    "JuliaCall, Julia, and local HSquared.jl are required for live Mrode9 Ainv validation."
  )

  fixture <- hsquared:::hs_mrode9_pedigree_validation_fixture()

  hs_skip_live_julia()
  hsquared:::hs_julia_setup(hsquared:::hs_default_julia_project())
  JuliaCall::julia_assign("hsq_mrode_id", fixture$pedigree$id)
  JuliaCall::julia_assign(
    "hsq_mrode_sire",
    hsquared:::hs_parent_for_julia(fixture$pedigree$sire)
  )
  JuliaCall::julia_assign(
    "hsq_mrode_dam",
    hsquared:::hs_parent_for_julia(fixture$pedigree$dam)
  )
  JuliaCall::julia_command(paste(
    "hsq_mrode_ped = HSquared.normalize_pedigree(",
    "hsq_mrode_id, hsq_mrode_sire, hsq_mrode_dam);",
    "hsq_mrode_Ainv = Matrix(HSquared.pedigree_inverse(hsq_mrode_ped));",
    "hsq_mrode_ids = hsq_mrode_ped.ids;"
  ))

  observed <- JuliaCall::julia_eval("hsq_mrode_Ainv")
  ids <- JuliaCall::julia_eval("hsq_mrode_ids")
  dimnames(observed) <- list(ids, ids)
  expected <- as.matrix(fixture$expected$Ainv[ids, ids])

  expect_equal(observed, expected, tolerance = 1e-10)
})
