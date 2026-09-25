checker <- system.file(
  "tools", "v07_genomic_recovery_recompute.R", package = "hsquared"
)
if (!nzchar(checker)) {
  checker <- testthat::test_path(
    "..", "..", "inst", "tools", "v07_genomic_recovery_recompute.R"
  )
}
source(checker, local = TRUE)

v07_test_hash <- function(letter) paste(rep(letter, 64L), collapse = "")

v07_test_inputs <- function(tier = "confirm") {
  seeds <- 1:5
  raw <- data.frame(
    tier = tier, cell_id = "n120_m600_r020", seed = seeds, n = 120, m = 600,
    truth_sigma_g2 = 0.2, truth_sigma_e2 = 0.8, truth_ratio = 0.2,
    ridge = 0.01,
    estimate_sigma_g2 = c(0.19, 0.20, 0.21, 0.20, NA),
    estimate_sigma_e2 = c(0.79, 0.80, 0.81, 0.80, NA),
    estimate_ratio = c(0.194, 0.200, 0.206, 0.200, NA),
    converged = c(TRUE, TRUE, TRUE, TRUE, FALSE), iterations = c(8, 8, 9, 8, -1),
    objective = c(1, 2, 3, 4, NA), gradient_norm = c(1e-8, 1e-8, 1e-8, 1e-8, NA),
    runtime_seconds = 1, peak_rss_mb = 10,
    marker_hash = vapply(seeds, function(i) v07_test_hash(as.character(i)), character(1)),
    id_hash = v07_test_hash("a"),
    kernel_hash = vapply(seeds, function(i) v07_test_hash(letters[i]), character(1)),
    error_class = c(rep("none", 4), "fit_not_converged"),
    stringsAsFactors = FALSE
  )
  manifest <- raw[v07_manifest_columns[v07_manifest_columns %in% names(raw)]]
  manifest$ridge <- 0.01
  manifest$regime <- "marker_rich_n120"
  manifest <- manifest[v07_manifest_columns]
  list(raw = raw, manifest = manifest)
}

test_that("raw reader accepts real one-row TSV data frames", {
  root <- withr::local_tempdir()
  raw_dir <- file.path(root, "raw", "confirm", "n120_m600_r020")
  dir.create(raw_dir, recursive = TRUE)
  x <- v07_test_inputs("confirm")
  for (i in 1:2) {
    utils::write.table(
      x$raw[i, , drop = FALSE], file.path(raw_dir, paste0(i, ".tsv")),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
  }
  read_back <- v07_read_raw(root, "confirm")
  expect_equal(nrow(read_back), 2L)
  expect_identical(names(read_back), v07_result_columns)

  utils::write.table(
    x$raw[1:2, , drop = FALSE], file.path(raw_dir, "bad.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  expect_error(v07_read_raw(root, "confirm"), "exactly one row")
})

test_that("independent R summary reproduces hand-calculated statistics", {
  x <- v07_test_inputs("confirm")
  out <- v07_recompute_summary(x$raw, x$manifest, "confirm", "n120_m600_r020")

  expect_equal(out$n_attempted, rep(5, 3))
  expect_equal(out$n_converged, rep(4, 3))
  expect_equal(out$mean_estimate, c(0.20, 0.80, 0.20), tolerance = 1e-14)
  expect_equal(out$bias, c(0, 0, 0), tolerance = 1e-14)
  expect_equal(out$mcse, c(sd(c(.19, .20, .21, .20)), sd(c(.79, .80, .81, .80)), sd(c(.194, .2, .206, .2))) / 2)
  expect_equal(out$convergence_rate, rep(0.8, 3))
  expect_match(out$failure_classes, "fit_not_converged=1;none=4", fixed = TRUE)
  expect_identical(out$cell_status, rep("FAIL", 3))
})

test_that("R versus Julia comparison uses a 1e-10 numeric gate", {
  x <- v07_test_inputs("confirm")
  julia <- v07_recompute_summary(x$raw, x$manifest, "confirm", "n120_m600_r020")
  expect_invisible(v07_compare_summary(julia, julia, 1e-10))
  within <- julia
  within$bias[[1]] <- within$bias[[1]] + 0.9e-10
  expect_invisible(v07_compare_summary(julia, within, 1e-10))
  outside <- julia
  outside$bias[[1]] <- outside$bias[[1]] + 1.1e-10
  expect_error(v07_compare_summary(julia, outside, 1e-10), "summary mismatch in bias")
})

test_that("pilot sizing uses the preregistered upper SD bound", {
  x <- v07_test_inputs("pilot")
  # Make every row eligible so the numeric anchor has df = 4.
  x$raw$converged <- TRUE
  x$raw$error_class <- "none"
  x$raw$estimate_sigma_g2[[5]] <- 0.22
  x$raw$estimate_sigma_e2[[5]] <- 0.82
  x$raw$estimate_ratio[[5]] <- 0.212
  out <- v07_recompute_summary(x$raw, x$manifest, "pilot", "n120_m600_r020")

  values <- c(0.19, 0.20, 0.21, 0.20, 0.22)
  raw_sd <- sd(values)
  upper_sd <- raw_sd * sqrt(4 / qchisq(0.05, df = 4))
  expected_n <- ceiling((qnorm(0.975) * upper_sd / (0.01 / 2))^2)
  naive_n <- ceiling((qnorm(0.975) * raw_sd / (0.01 / 2))^2)

  expect_gt(upper_sd, raw_sd)
  expect_gt(expected_n, naive_n)
  expect_equal(out$pilot_sd_upper[out$target == "sigma_g2"], upper_sd)
  expect_equal(out$required_n_raw[out$target == "sigma_g2"], expected_n)
})

test_that("every preregistered mutation turns at least one gate red", {
  x <- v07_test_inputs("confirm")
  expected <- "n120_m600_r020"
  julia <- v07_recompute_summary(x$raw, x$manifest, "confirm", expected)

  estimate <- x$raw
  estimate$estimate_ratio[[1]] <- estimate$estimate_ratio[[1]] + 0.1
  recomputed <- v07_recompute_summary(estimate, x$manifest, "confirm", expected)
  expect_error(v07_compare_summary(recomputed, julia), "summary mismatch")

  truth <- x$raw
  truth$truth_ratio[[1]] <- 0.3
  expect_error(v07_validate_inputs(truth, x$manifest, "confirm", expected), "raw/manifest mismatch")

  duplicate <- rbind(x$raw, x$raw[1, ])
  expect_error(v07_validate_inputs(duplicate, x$manifest, "confirm", expected), "duplicate")

  removed_failed <- x$raw[x$raw$converged, ]
  expect_error(v07_validate_inputs(removed_failed, x$manifest, "confirm", expected), "exactly equal")

  cell_label <- x$manifest
  cell_label$cell_id[[1]] <- "wrong_cell"
  expect_error(v07_validate_inputs(x$raw, cell_label, "confirm", expected), "cell labels")

  ridge <- x$manifest
  ridge$ridge[[1]] <- 0.02
  expect_error(v07_validate_inputs(x$raw, ridge, "confirm", expected), "frozen ridge")

  raw_ridge <- x$raw
  raw_ridge$ridge[[1]] <- 0.02
  expect_error(v07_validate_inputs(raw_ridge, x$manifest, "confirm", expected), "raw/manifest mismatch")

  marker_hash <- x$raw
  marker_hash$marker_hash[[1]] <- "mutated"
  expect_error(v07_validate_inputs(marker_hash, x$manifest, "confirm", expected), "SHA-256")

  id_order <- x$raw
  id_order$id_hash[[1]] <- v07_test_hash("b")
  expect_error(v07_validate_inputs(id_order, x$manifest, "confirm", expected), "ID-order")

  pilot <- x$manifest
  pilot$tier <- "pilot"
  confirm <- x$manifest
  confirm$tier <- "confirm"
  expect_error(v07_validate_disjoint_seeds(pilot, confirm), "overlap")
})

test_that("raw checksum sealing detects a valid-looking hash mutation", {
  skip_if(!v07_has_sha256_command())
  root <- withr::local_tempdir()
  dir.create(file.path(root, "raw", "confirm"), recursive = TRUE)
  dir.create(file.path(root, "raw", "pilot"), recursive = TRUE)
  x <- v07_test_inputs("confirm")
  confirm_path <- file.path(root, "raw", "confirm", "1.tsv")
  pilot_path <- file.path(root, "raw", "pilot", "1.tsv")
  utils::write.table(x$raw[1, ], confirm_path, sep = "\t", quote = FALSE, row.names = FALSE)
  pilot_row <- x$raw[1, ]
  pilot_row$tier <- "pilot"
  utils::write.table(pilot_row, pilot_path, sep = "\t", quote = FALSE, row.names = FALSE)

  v07_write_raw_lock(root, "pilot")
  v07_write_raw_lock(root, "confirm")
  expect_true(file.exists(file.path(root, "pilot_raw_sha256.tsv")))
  expect_true(file.exists(file.path(root, "confirm_raw_sha256.tsv")))
  expect_invisible(v07_verify_raw_lock(root, "pilot"))
  expect_invisible(v07_verify_raw_lock(root, "confirm"))
  expect_error(v07_write_raw_lock(root, "pilot"), "refusing to reseal")
  expect_error(v07_write_raw_lock(root, "confirm"), "refusing to reseal")

  changed <- x$raw[1, ]
  changed$marker_hash <- v07_test_hash("f")
  utils::write.table(changed, confirm_path, sep = "\t", quote = FALSE, row.names = FALSE)
  expect_error(v07_verify_raw_lock(root, "confirm"), "checksum lock")
  expect_invisible(v07_verify_raw_lock(root, "pilot"))

  extra_path <- file.path(root, "raw", "pilot", "extra.tsv")
  utils::write.table(pilot_row, extra_path, sep = "\t", quote = FALSE, row.names = FALSE)
  expect_error(v07_verify_raw_lock(root, "pilot"), "file set differs")
  unlink(extra_path)
  unlink(pilot_path)
  expect_error(v07_verify_raw_lock(root, "pilot"), "file set differs")
})

test_that("raw-lock SHA-256 parser accepts the Windows certutil form", {
  expected <- "8bc191bf660c0c26d0e51f0f243502fe601c648f3c97a2047fa2c2b62bc3884f"
  certutil_output <- c(
    "SHA256 hash of file C:\\raw\\confirm\\1.tsv:",
    "8b c1 91 bf 66 0c 0c 26 d0 e5 1f 0f 24 35 02 fe 60 1c 64 8f 3c 97 a2 04 7f a2 c2 b6 2b c3 88 4f",
    "CertUtil: -hashfile command completed successfully."
  )

  expect_identical(v07_extract_sha256(certutil_output), expected)
})
