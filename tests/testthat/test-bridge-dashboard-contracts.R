# Bridge dashboard contract tests (A14 phase 1).
# Pattern donor: drmTMB test-structured-re-conversion-contracts.R (dashboard TSV reads).

bridge_dashboard_path <- function(file) {
  candidates <- c(
    file.path(testthat::test_path(), "..", "..", "docs", "dev-log", "dashboard", file),
    file.path(getwd(), "docs", "dev-log", "dashboard", file)
  )
  found <- candidates[file.exists(candidates)]
  if (length(found) == 0L) {
    testthat::skip(paste("dashboard source file not available:", file))
  }
  found[[1L]]
}

bridge_read_dashboard_tsv <- function(file) {
  path <- bridge_dashboard_path(file)
  lines <- readLines(path, warn = FALSE)
  if (length(lines) < 2L) {
    return(data.frame())
  }
  utils::read.delim(
    text = paste(lines, collapse = "\n"),
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

test_that("bridge dashboard TSV column contracts match the A14 schema", {
  schema <- bridge_read_dashboard_tsv("bridge-payload-schema.tsv")
  parity <- bridge_read_dashboard_tsv("bridge-parity-smoke-status.tsv")
  boundary <- bridge_read_dashboard_tsv("bridge-boundary.tsv")

  expect_named(
    schema,
    c(
      "schema_id", "route", "r_target", "julia_dispatch", "estimator",
      "payload_fields", "relinv_builder", "r_bridge_status", "julia_status",
      "claim_boundary", "evidence_url", "next_gate"
    )
  )
  expect_named(
    parity,
    c(
      "smoke_id", "schema_id", "model_cell", "r_path", "julia_path",
      "parity_target", "tolerance_rule", "parity_status", "test_status",
      "bridge_status", "evidence_url", "claim_boundary", "next_gate"
    )
  )
  expect_named(
    boundary,
    c(
      "boundary_id", "target", "smoke_status", "parity_required",
      "bridge_status", "boundary_doc_status", "evidence_url", "claim_boundary",
      "next_gate"
    )
  )

  expect_true(nrow(schema) >= 5L)
  expect_true(nrow(parity) >= 5L)
  expect_true(nrow(boundary) >= 5L)
})

test_that("bridge dashboard rows carry claim_boundary and link schema_ids to tests", {
  schema <- bridge_read_dashboard_tsv("bridge-payload-schema.tsv")
  parity <- bridge_read_dashboard_tsv("bridge-parity-smoke-status.tsv")
  boundary <- bridge_read_dashboard_tsv("bridge-boundary.tsv")

  expect_true(all(nzchar(schema$claim_boundary)))
  expect_true(all(nzchar(parity$claim_boundary)))
  expect_true(all(nzchar(boundary$claim_boundary)))

  expect_true(all(parity$schema_id %in% schema$schema_id))

  public_covered <- c(
    "v02_two_effect",
    "v02_multi_effect",
    "v02_direct_maternal"
  )
  expect_true(all(public_covered %in% schema$schema_id))
  expect_true(all(schema$r_bridge_status[schema$schema_id %in% public_covered] == "covered"))

  evidence_files <- unique(c(
    schema$evidence_url,
    parity$evidence_url,
    boundary$evidence_url[boundary$boundary_doc_status == "documented"]
  ))
  evidence_files <- evidence_files[!grepl("^HSquared\\.jl/", evidence_files)]
  repo_root <- normalizePath(file.path(testthat::test_path(), "..", ".."), mustWork = TRUE)
  for (ref in evidence_files) {
    path <- file.path(repo_root, ref)
    expect_true(file.exists(path), info = paste("missing evidence:", ref))
  }
})

test_that("boundary_doc_status vocabulary cannot be read as bridge coverage", {
  boundary <- bridge_read_dashboard_tsv("bridge-boundary.tsv")

  doc_statuses <- c("documented", "partial", "planned")
  expect_true(all(boundary$boundary_doc_status %in% doc_statuses))
  # The two columns must not share the words that mean "covered", or a reader
  # scanning either one will count a documented boundary as a covered surface.
  expect_false("covered" %in% boundary$boundary_doc_status)
  expect_false("documented" %in% boundary$bridge_status)

  no_smoke <- boundary[boundary$smoke_status == "no_smoke", , drop = FALSE]
  expect_true(nrow(no_smoke) > 0L)
  expect_false(any(no_smoke$bridge_status == "covered"))
})

test_that("Tier 0 parity smoke rows stay Julia-free with covered test_status", {
  parity <- bridge_read_dashboard_tsv("bridge-parity-smoke-status.tsv")
  tier0_ids <- c(
    "smoke_payload_v2_emitter",
    "smoke_bridge_payload_v01",
    "smoke_gryphon_r_reference"
  )
  rows <- parity[parity$smoke_id %in% tier0_ids, , drop = FALSE]
  expect_equal(nrow(rows), length(tier0_ids))
  expect_true(all(rows$julia_path == "none"))
  expect_true(all(rows$test_status == "covered"))
})

test_that("bridge dashboard validator passes on the committed ledgers", {
  repo_root <- normalizePath(file.path(testthat::test_path(), "..", ".."), mustWork = TRUE)
  validator <- file.path(repo_root, "tools", "validate-bridge-dashboard.py")
  skip_if_not(file.exists(validator), "validate-bridge-dashboard.py not present")

  out_file <- tempfile(fileext = ".txt")
  on.exit(unlink(out_file), add = TRUE)
  status <- system2(
    "python3",
    args = c(validator),
    stdout = out_file,
    stderr = out_file
  )
  result <- readLines(out_file, warn = FALSE)
  expect_equal(status, 0L)
  expect_match(paste(result, collapse = "\n"), "bridge_dashboard_ok")
})
