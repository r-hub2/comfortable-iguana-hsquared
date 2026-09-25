# C1-ext / H1+H3 R PATH_ONLY pointer (Julia #294 twin).
# Fast fixture/contract test: no Julia, no coverage numbers, no flip.
# Optional extra smoke (same path, still not a claim):
#   HSQUARED_C1EXT_SMOKE=1 Rscript sim/phase1_interval_coverage_ext.R
# sim/ is .Rbuildignore'd, so the helper is source-tree only (same pattern as
# test-v07-genomic-boundary-oracle.R and data-raw comparator skips).

ext_helper <- testthat::test_path(
  "..",
  "..",
  "sim",
  "phase1_interval_coverage_ext.R"
)
if (file.exists(ext_helper)) {
  source(normalizePath(ext_helper, mustWork = TRUE), local = TRUE)
}

test_that("C1-ext H1/H3 R pointer freezes Julia #294 campaign contract", {
  skip_if_not(
    file.exists(ext_helper),
    "sim/phase1_interval_coverage_ext.R is not in the build tarball"
  )
  expect_equal(EXT_INTERPRETABLE_FRACTION, 0.9)
  expect_equal(EXT_CONFIRM_REPS_TARGET, 2000L)
  expect_equal(EXT_SEED_STRIDE, 40009L)
  expect_equal(EXT_PROMOTABLE_LEVEL, 0.95)
  expect_setequal(
    EXT_CAMPAIGNS,
    c("h1_two", "h1_multi", "h1_t", "h3_rg", "h3_ram")
  )

  expect_true(all(
    SYMBOLIC_ALIGNMENT$role[SYMBOLIC_ALIGNMENT$estimand == "t"] ==
      "characterization_only"
  ))
  expect_true(all(
    SYMBOLIC_ALIGNMENT$role[SYMBOLIC_ALIGNMENT$estimand == "ratio1"] ==
      "covered_pillar_bank"
  ))
  expect_true(all(
    SYMBOLIC_ALIGNMENT$role[
      SYMBOLIC_ALIGNMENT$estimand %in% c("r_g", "r_am")
    ] ==
      "covered_pillar_bank"
  ))
  expect_true(all(SYMBOLIC_ALIGNMENT$campaign %in% EXT_CAMPAIGNS))
})

test_that("C1-ext R parser is PATH_ONLY smoke and rejects claim modes", {
  skip_if_not(
    file.exists(ext_helper),
    "sim/phase1_interval_coverage_ext.R is not in the build tarball"
  )
  cfg <- hs_c1ext_parse_args(character())
  expect_equal(cfg$mode, "smoke")
  expect_equal(cfg$campaigns, EXT_CAMPAIGNS)
  expect_match(cfg$output, "c1ext-r-smoke.tsv", fixed = TRUE)

  cfg2 <- hs_c1ext_parse_args(c("--mode=smoke", "--campaigns=h1_two,h3_rg"))
  expect_equal(cfg2$campaigns, c("h1_two", "h3_rg"))

  expect_error(
    hs_c1ext_parse_args("--campaigns=h1_nbinom"),
    "Unknown campaigns"
  )
  expect_error(hs_c1ext_parse_args("--mode=promote"), "Unsupported --mode")
  expect_error(
    hs_c1ext_parse_args("--mode=confirm"),
    "not available on the R twin"
  )
  expect_error(
    hs_c1ext_parse_args("--mode=screen"),
    "not available on the R twin"
  )
})

test_that("C1-ext smoke TSV is PATH_ONLY and never claim-eligible", {
  skip_if_not(
    file.exists(ext_helper),
    "sim/phase1_interval_coverage_ext.R is not in the build tarball"
  )
  path <- tempfile(fileext = ".tsv")
  on.exit(unlink(path), add = TRUE)
  written <- hs_c1ext_write_smoke_tsv(path)
  expect_true(file.exists(written))

  rows <- utils::read.delim(written, stringsAsFactors = FALSE)
  expect_true(all(
    c("claim_eligible", "gate", "public_covered_count") %in% names(rows)
  ))
  expect_true(all(rows$gate == "PATH_ONLY"))
  expect_true(all(rows$claim_eligible == FALSE))
  expect_true(all(rows$public_covered_count == 7L))
  expect_true(all(rows$experimental_version == "0.9.0"))
  expect_true(all(rows$campaign %in% EXT_CAMPAIGNS))
  expect_equal(sum(rows$role == "characterization_only"), 1L)
})

test_that("R interval claim surface stays experimental / blocked for H1/H3", {
  expect_true(exists("heritability_interval", mode = "function"))
  expect_true(exists("repeatability_interval", mode = "function"))
  expect_true(exists("common_env_proportion_interval", mode = "function"))
  expect_true(exists("maternal_proportion_interval", mode = "function"))
  expect_true(exists("genetic_correlation", mode = "function"))
  expect_false(exists(
    "genetic_correlation_interval",
    mode = "function",
    inherits = TRUE
  ))

  desc <- read.dcf(system.file("DESCRIPTION", package = "hsquared"))
  expect_equal(unname(desc[, "Version"]), "0.9.0")
  expect_match(
    unname(desc[, "Description"]),
    "public covered count is 7",
    fixed = TRUE
  )

  fit <- hsquared:::hs_new_fit(
    call = quote(hsquared(y ~ animal(1 | id, pedigree = ped), data = dat)),
    spec = list(method = "REML", family = list(family = "gaussian")),
    payload = list(y = seq_len(10)),
    result = list(heritability = data.frame(term = "animal", estimate = 0.4))
  )
  expect_error(stats::confint(fit), "Validated confidence intervals")
  expect_error(
    stats::vcov(fit),
    "variance-covariance matrix is not implemented"
  )
  expect_error(stats::profile(fit), "Profile-likelihood intervals")
})
