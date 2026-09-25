test_that("real-data validation manifest defines 3-tier ladder with claim boundaries (A13)", {
  hs_skip_without_realdata_manifest()
  manifest <- hs_read_realdata_manifest()
  expect_identical(manifest$schema_version, 1L)
  expect_identical(manifest$lane, "hsquared")
  expect_match(
    manifest$claim_boundary,
    "does not promote any capability to covered",
    ignore.case = TRUE
  )
  expect_match(
    manifest$claim_boundary,
    "does not claim field-empirical",
    ignore.case = TRUE
  )

  expect_true(!is.null(manifest$tier_summary))
  expect_match(
    manifest$tier_summary$tier_1,
    "Not empirical validation",
    fixed = FALSE
  )
  expect_match(
    manifest$tier_summary$tier_2,
    "Not field empirical",
    fixed = FALSE
  )
  expect_match(manifest$tier_summary$tier_4, "NOT STARTED", fixed = FALSE)

  ids <- hs_realdata_arc_ids(manifest)
  expect_length(unique(ids), length(ids))

  tiers_present <- sort(unique(vapply(manifest$arc, `[[`, integer(1), "tier")))
  expect_identical(tiers_present, c(1L, 2L, 3L, 4L))

  tier_counts <- table(vapply(manifest$arc, `[[`, integer(1), "tier"))
  expect_true(tier_counts[["1"]] >= 5L)
  expect_true(tier_counts[["2"]] >= 4L)
  expect_true(tier_counts[["3"]] >= 3L)
  expect_identical(as.integer(tier_counts[["4"]]), 1L)

  allowed_data_origins <- c(
    "synthetic",
    "teaching-simulated",
    "published-textbook",
    "generated-comparator",
    "field-empirical"
  )
  allowed_darwin <- c("pending", "signed", "blocked")
  allowed_ci <- c("in-CI", "skip-guarded", "opt-in")

  for (arc in manifest$arc) {
    expect_true(nzchar(arc$id))
    expect_true(arc$data_origin %in% allowed_data_origins)
    expect_true(arc$darwin_review %in% allowed_darwin)
    expect_true(arc$ci_policy %in% allowed_ci)
    expect_true(nzchar(arc$claim_boundary))
    if (arc$tier < 4L) {
      expect_match(
        arc$claim_boundary,
        "not|no|NOT|only|open|Internal|Test-of-test|FROZEN",
        ignore.case = TRUE
      )
    }

    if (arc$tier < 4L) {
      expect_false(arc$data_origin == "field-empirical")
    }
  }

  gryphon <- manifest$arc[[match("gryphon_bwt_reml", ids)]]
  expect_identical(gryphon$tier, 2L)
  expect_identical(gryphon$data_origin, "teaching-simulated")
  expect_match(gryphon$claim_boundary, "Not field empirical", fixed = FALSE)

  placeholder <- manifest$arc[[match("field_empirical_placeholder", ids)]]
  expect_identical(placeholder$tier, 4L)
  expect_identical(placeholder$darwin_review, "blocked")
  expect_match(placeholder$claim_boundary, "NOT STARTED", fixed = FALSE)
})

test_that("Darwin review stub is pending with checklist items (A13)", {
  hs_skip_without_realdata_manifest()
  manifest <- hs_read_realdata_manifest()
  review <- manifest$darwin_review

  expect_identical(review$status, "pending")
  expect_identical(review$signed_date, "")
  expect_identical(review$reviewer, "")
  expect_true(length(review$questions) >= 3L)
  expect_true(length(review$checklist) >= 3L)

  checklist_statuses <- vapply(review$checklist, `[[`, character(1), "status")
  expect_true(all(checklist_statuses == "pending"))

  gryphon_q <- vapply(
    review$questions,
    function(q) "gryphon_bwt_reml" %in% q$arc_ids,
    logical(1)
  )
  expect_true(any(gryphon_q))
})
