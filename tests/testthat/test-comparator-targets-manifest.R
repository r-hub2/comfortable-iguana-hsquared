test_that("comparator target manifest indexes mirrored R fixtures (#49 coordination)", {
  manifest <- hs_read_comparator_manifest()
  expect_identical(manifest$schema_version, 1L)
  expect_identical(manifest$lane, "hsquared")
  expect_match(manifest$claim_boundary, "not add external comparator evidence", fixed = FALSE)

  targets <- manifest$target
  ids <- vapply(targets, `[[`, character(1), "id")
  expect_length(unique(ids), length(ids))
  expect_setequal(
    ids,
    c(
      "animal_model_fitted_target",
      "sire_model_fitted_target",
      "phase4_multitrait_parity",
      "genomic_gblup_snpblup_target",
      "marker_scan_parity",
      "structured_covariance_parity",
      "non_gaussian_parity"
    )
  )

  allowed_evidence <- c(
    "julia_target",
    "julia_target_r_consumed",
    "julia_target_external_one_leg",
    "bridge_payload_fixture"
  )

  pkg_root <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)

  for (target in targets) {
    expect_true(length(target$capability_rows) > 0L)
    expect_true(target$evidence_type %in% allowed_evidence)
    expect_match(target$boundary, "no|not", ignore.case = TRUE)
    expect_true(nzchar(trimws(target$required_comparator)))

    if (isTRUE(target$r_mirror)) {
      fixture_dir <- file.path(pkg_root, target$r_fixture_path)
      expect_true(dir.exists(fixture_dir))
      for (file in target$required_files) {
        expect_true(file.exists(file.path(fixture_dir, file)))
      }
    } else {
      expect_false(target$r_mirror)
      expect_identical(target$r_fixture_path, "")
    }
  }

  sire <- targets[[match("sire_model_fitted_target", ids)]]
  expect_identical(sire$r_mirror, FALSE)
  expect_identical(sire$issue, 16L)

  multivariate <- targets[[match("phase4_multitrait_parity", ids)]]
  expect_match(multivariate$external_status, "sommer")
  expect_match(multivariate$boundary, "second independent")

  genomic <- targets[[match("genomic_gblup_snpblup_target", ids)]]
  expect_identical(genomic$evidence_type, "julia_target_r_consumed")
  expect_match(genomic$external_status, "PR #84")

  nongaussian <- targets[[match("non_gaussian_parity", ids)]]
  expect_match(nongaussian$boundary, "no per-record varying-trial R activation")
})

test_that("frozen comparator fixture SHA256 pins match canonical CSV content", {
  pins <- utils::read.csv(
    hs_comparator_fixture_shas_path(),
    stringsAsFactors = FALSE
  )
  pkg_root <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
  manifest <- hs_read_comparator_manifest()
  mirrored_ids <- vapply(
    manifest$target,
    function(t) if (isTRUE(t$r_mirror)) t$id else NA_character_,
    character(1)
  )
  mirrored_ids <- mirrored_ids[!is.na(mirrored_ids)]

  for (id in mirrored_ids) {
    target <- manifest$target[[match(id, vapply(manifest$target, `[[`, "", "id"))]]
    fixture_dir <- file.path(pkg_root, target$r_fixture_path)
    rows <- pins[pins$target == id, , drop = FALSE]
    expect_true(nrow(rows) > 0L)
    for (i in seq_len(nrow(rows))) {
      path <- file.path(fixture_dir, rows$file[i])
      expect_true(file.exists(path))
      expect_identical(hs_sha256_canonical_csv(path), rows$sha256[i])
    }
  }
})

test_that("canonical CSV hashes ignore CRLF checkout conversion", {
  lf <- tempfile(fileext = ".csv")
  crlf <- tempfile(fileext = ".csv")
  on.exit(unlink(c(lf, crlf)), add = TRUE)
  writeBin(charToRaw("a,b\n1,2\n"), lf)
  writeBin(charToRaw("a,b\r\n1,2\r\n"), crlf)

  expect_identical(hs_sha256_canonical_csv(lf), hs_sha256_canonical_csv(crlf))
})

test_that("SHA-256 output parser accepts the Windows certutil form", {
  expected <- "8bc191bf660c0c26d0e51f0f243502fe601c648f3c97a2047fa2c2b62bc3884f"
  certutil_output <- c(
    "SHA256 hash of file C:\\fixtures\\sealed.csv:",
    "8b c1 91 bf 66 0c 0c 26 d0 e5 1f 0f 24 35 02 fe 60 1c 64 8f 3c 97 a2 04 7f a2 c2 b6 2b c3 88 4f",
    "CertUtil: -hashfile command completed successfully."
  )

  expect_identical(hs_extract_sha256(certutil_output), expected)
})

test_that("fixture SHA-256 helper falls through a failed backend", {
  expected <- "8bc191bf660c0c26d0e51f0f243502fe601c648f3c97a2047fa2c2b62bc3884f"
  find_command <- function(command) {
    if (command %in% c("shasum", "sha256sum")) command else ""
  }
  run_command <- function(command, ...) {
    if (identical(command, "shasum")) {
      return(structure("Git shasum failed", status = 29L))
    }
    paste(expected, "sealed.csv")
  }

  expect_identical(
    hs_sha256_file("sealed.csv", find_command = find_command, run_command = run_command),
    expected
  )
})

test_that("fixture SHA-256 helper uses the Windows system certutil path", {
  expected <- "8bc191bf660c0c26d0e51f0f243502fe601c648f3c97a2047fa2c2b62bc3884f"
  commands <- character()
  run_command <- function(command, ...) {
    commands <<- c(commands, command)
    c(
      "SHA256 hash of file C:\\fixtures\\sealed.csv:",
      "8b c1 91 bf 66 0c 0c 26 d0 e5 1f 0f 24 35 02 fe 60 1c 64 8f 3c 97 a2 04 7f a2 c2 b6 2b c3 88 4f",
      "CertUtil: -hashfile command completed successfully."
    )
  }

  expect_identical(
    hs_sha256_file(
      "sealed.csv",
      find_command = function(...) "",
      run_command = run_command,
      windows_directory = "C:/Windows"
    ),
    expected
  )
  expect_identical(commands, "C:/Windows/System32/certutil.exe")
})
