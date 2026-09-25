# A16 phase 1 — Tier 0 bridge CI contracts (fixture-first, no live Julia).
# Plan: ~/local-scratch/h2-b4-bridge-plan.md § A16 Tier 0
# Registry: docs/dev-log/dashboard/bridge-ci-tier0.md

bridge_tier0_repo_root <- function() {
  normalizePath(file.path(testthat::test_path(), "..", ".."), mustWork = TRUE)
}

# `docs/` is .Rbuildignore'd, so the dashboard registry is absent under
# `R CMD check`, which installs from the tarball. These are source-tree
# contract tests.
bridge_tier0_skip_without_docs <- function() {
  testthat::skip_if_not(
    dir.exists(file.path(
      bridge_tier0_repo_root(),
      "docs",
      "dev-log",
      "dashboard"
    )),
    "docs/dev-log/dashboard is not in the build tarball"
  )
}

bridge_tier0_read_parity_tsv <- function() {
  path <- file.path(
    bridge_tier0_repo_root(),
    "docs",
    "dev-log",
    "dashboard",
    "bridge-parity-smoke-status.tsv"
  )
  lines <- readLines(path, warn = FALSE)
  utils::read.delim(
    text = paste(lines, collapse = "\n"),
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

bridge_tier0_payload_fingerprint <- function(payload) {
  re_names <- vapply(
    payload$random_effects,
    function(blk) blk$name,
    character(1L)
  )
  re_types <- vapply(
    payload$random_effects,
    function(blk) blk$type,
    character(1L)
  )
  paste(
    payload$payload_version,
    length(payload$y),
    nrow(payload$X),
    ncol(payload$X),
    paste(re_names, collapse = ","),
    paste(re_types, collapse = ","),
    sep = "|"
  )
}

test_that("Tier 0 registry documents Julia-free parity smoke rows", {
  bridge_tier0_skip_without_docs()
  parity <- bridge_tier0_read_parity_tsv()
  tier0_ids <- c(
    "smoke_payload_v2_emitter",
    "smoke_bridge_payload_v01",
    "smoke_gryphon_r_reference"
  )

  rows <- parity[parity$smoke_id %in% tier0_ids, , drop = FALSE]
  expect_equal(nrow(rows), length(tier0_ids))
  expect_equal(sort(rows$smoke_id), sort(tier0_ids))
  expect_true(all(rows$julia_path == "none"))
  expect_true(all(rows$test_status == "covered"))

  repo_root <- bridge_tier0_repo_root()
  for (ref in rows$evidence_url) {
    expect_true(
      file.exists(file.path(repo_root, ref)),
      info = paste("missing:", ref)
    )
  }
})

test_that("Tier 0 canonical test files exist on disk", {
  repo_root <- bridge_tier0_repo_root()
  tier0_tests <- c(
    "tests/testthat/test-bridge-payload-v2.R",
    "tests/testthat/test-bridge-payload.R",
    "tests/testthat/test-julia-bridge.R",
    "tests/testthat/test-validation-fixtures.R",
    "tests/testthat/test-bridge-dashboard-contracts.R",
    "tests/testthat/test-bridge-ci-tier0-contracts.R"
  )
  for (path in tier0_tests) {
    expect_true(file.exists(file.path(repo_root, path)), info = path)
  }

  bridge_tier0_skip_without_docs()
  tier0_doc <- file.path(
    repo_root,
    "docs",
    "dev-log",
    "dashboard",
    "bridge-ci-tier0.md"
  )
  expect_true(file.exists(tier0_doc))
  doc_lines <- readLines(tier0_doc, warn = FALSE)
  expect_true(any(grepl("h2-b4-bridge-plan", doc_lines, fixed = TRUE)))
})

test_that("Tier 0 payload v2 emitter fingerprints are stable (in-memory)", {
  ped <- data.frame(
    id = c("a", "b", "c", "d"),
    sire = c(NA, NA, "a", "a"),
    dam = c(NA, NA, "b", "c"),
    stringsAsFactors = FALSE
  )
  dat <- data.frame(
    y = c(1.0, 2.0, 3.0),
    sex = c("f", "m", "f"),
    id = c("a", "c", "d"),
    stringsAsFactors = FALSE
  )
  spec <- hsquared:::hs_build_model_spec(
    y ~ sex + animal(1 | id, pedigree = ped),
    data = dat,
    family = stats::gaussian(),
    REML = TRUE
  )
  payload <- hsquared:::hs_build_bridge_payload(spec)

  expect_equal(payload$payload_version, 2L)
  expect_equal(
    bridge_tier0_payload_fingerprint(payload),
    "2|3|3|2|animal|pedigree"
  )

  dat_b <- data.frame(
    y = c(14.0, 13.0, 12.1, 8.9),
    id = c("a", "b", "c", "d"),
    litter = c("L1", "L1", "L2", "L2"),
    stringsAsFactors = FALSE
  )
  spec_b <- hsquared:::hs_build_model_spec(
    y ~ animal(1 | id, pedigree = ped) + common_env(1 | litter),
    data = dat_b,
    family = stats::gaussian(),
    REML = TRUE
  )
  payload_b <- hsquared:::hs_build_bridge_payload(spec_b)
  expect_equal(
    bridge_tier0_payload_fingerprint(payload_b),
    "2|4|4|1|animal,common_env|pedigree,iid"
  )
})

test_that("Tier 0 emit_payload_v2_fixtures.R script is present and syntactically valid", {
  repo_root <- bridge_tier0_repo_root()
  script <- file.path(
    repo_root,
    "tests",
    "fixtures",
    "emit_payload_v2_fixtures.R"
  )
  expect_true(file.exists(script))
  parsed <- parse(file = script)
  expect_true(length(parsed) > 0L)
})
