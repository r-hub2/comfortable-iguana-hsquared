local_comparator_script <- function(...) {
  candidates <- c(
    testthat::test_path(
      "..",
      "..",
      "inst",
      "comparator-scripts",
      ...
    ),
    system.file("comparator-scripts", ..., package = "hsquared")
  )
  candidates <- candidates[nzchar(candidates) & file.exists(candidates)]
  testthat::skip_if(
    length(candidates) == 0L,
    paste("Comparator script is not available:", file.path(...))
  )
  normalizePath(candidates[[1]], mustWork = TRUE)
}

local_package_root <- function() {
  normalizePath(testthat::test_path("..", ".."), mustWork = FALSE)
}

run_comparator_script <- function(script, args = character()) {
  rscript <- file.path(R.home("bin"), "Rscript")
  testthat::skip_if_not(
    file.exists(rscript),
    "Rscript executable is not available."
  )

  output <- system2(
    rscript,
    args = shQuote(c(script, args)),
    stdout = TRUE,
    stderr = TRUE,
    env = paste0("HSQUARED_REPO=", shQuote(local_package_root()))
  )
  status <- attr(output, "status")
  if (is.null(status)) {
    status <- 0L
  }
  list(status = status, output = paste(output, collapse = "\n"))
}

test_that("MCMCglmm agreement probe is reproducible and fenced", {
  script <- testthat::test_path(
    "..",
    "..",
    "data-raw",
    "multivariate-mcmcglmm-agreement-study.R"
  )
  testthat::skip_if_not(
    file.exists(script),
    "data-raw MCMCglmm agreement script is not available in built-package checks."
  )
  script <- normalizePath(
    testthat::test_path(
      "..",
      "..",
      "data-raw",
      "multivariate-mcmcglmm-agreement-study.R"
    ),
    mustWork = TRUE
  )
  lines <- readLines(script, warn = FALSE)
  text <- paste(lines, collapse = "\n")

  expect_match(text, "MCMCglmm::MCMCglmm", fixed = TRUE)
  expect_match(text, "set.seed(20260621)", fixed = TRUE)
  expect_match(text, "nitt = 50000", fixed = TRUE)
  expect_match(text, "burnin = 10000", fixed = TRUE)
  expect_match(text, "thin = 40", fixed = TRUE)
  expect_match(text, "targets_inside_95_HPD", fixed = TRUE)
  expect_match(text, "not a same-estimand REML comparator", fixed = TRUE)
  expect_match(text, "must not promote V4-MV-REML", fixed = TRUE)
})

test_that("manual ASReml comparator script has a dry-run mode", {
  script <- local_comparator_script("asreml", "multivariate-animal.R")

  run <- run_comparator_script(script)

  expect_identical(run$status, 0L)
  expect_match(run$output, "Prepared ASReml-R candidate fixture", fixed = TRUE)
  expect_match(run$output, "records: 160", fixed = TRUE)
  expect_match(run$output, "traits: trait1, trait2", fixed = TRUE)
  expect_match(run$output, "Dry run only", fixed = TRUE)
})

test_that("manual BLUPF90 comparator script has dry-run and write modes", {
  skip_if_not_installed("withr")
  script <- local_comparator_script("blupf90", "prepare-multivariate-animal.R")

  dry_run <- run_comparator_script(script)
  expect_identical(dry_run$status, 0L)
  expect_match(
    dry_run$output,
    "Prepared BLUPF90-family candidate fixture",
    fixed = TRUE
  )
  expect_match(dry_run$output, "data rows: 80", fixed = TRUE)
  expect_match(dry_run$output, "pedigree rows: 20", fixed = TRUE)
  expect_match(dry_run$output, "Dry run only", fixed = TRUE)

  write_dir <- withr::local_tempdir(pattern = "hsquared-blupf90-")
  write_run <- run_comparator_script(
    script,
    paste0("--write=", write_dir)
  )
  expect_identical(write_run$status, 0L)
  expect_match(
    write_run$output,
    "Wrote BLUPF90-family comparator files",
    fixed = TRUE
  )

  data_file <- file.path(write_dir, "multivariate-animal.dat")
  ped_file <- file.path(write_dir, "multivariate-animal.ped")
  renf90_file <- file.path(write_dir, "multivariate-animal.renf90")
  par_file <- file.path(write_dir, "multivariate-animal.par")
  readme_file <- file.path(write_dir, "README.txt")

  expect_true(file.exists(data_file))
  expect_true(file.exists(ped_file))
  expect_true(file.exists(renf90_file))
  expect_true(file.exists(par_file))
  expect_true(file.exists(readme_file))

  dat <- utils::read.table(data_file)
  ped <- utils::read.table(ped_file)
  expect_equal(nrow(dat), 80L)
  expect_equal(ncol(dat), 5L)
  expect_equal(nrow(ped), 20L)
  expect_equal(ncol(ped), 3L)

  renf90 <- readLines(renf90_file, warn = FALSE)
  par <- readLines(par_file, warn = FALSE)
  expect_false(any(grepl("__DATAFILE__", renf90, fixed = TRUE)))
  expect_false(any(grepl("__PEDFILE__", renf90, fixed = TRUE)))
  expect_false(any(grepl("__DATAFILE__", par, fixed = TRUE)))
  expect_false(any(grepl("__PEDFILE__", par, fixed = TRUE)))
  expect_true(any(grepl("multivariate-animal.dat", renf90, fixed = TRUE)))
  expect_true(any(grepl("multivariate-animal.ped", par, fixed = TRUE)))
})

test_that("BLUPF90 multivariate summary ingester validates synthetic reports", {
  skip_if_not_installed("withr")
  path <- withr::local_tempfile(fileext = ".csv")
  summary <- data.frame(
    quantity = c(
      "G[1,1]",
      "G[1,2]",
      "G[2,2]",
      "R[1,1]",
      "R[1,2]",
      "R[2,2]",
      "h2 trait 1",
      "h2 trait 2",
      "EBV correlation trait 1",
      "EBV correlation trait 2"
    ),
    target = c(1, 0.2, 1.5, 0.7, 0.1, 0.8, 0.58, 0.65, NA, NA),
    estimate = c(
      1.0001,
      0.2001,
      1.4999,
      0.7001,
      0.1001,
      0.8001,
      0.5801,
      0.6499,
      0.9995,
      0.9994
    ),
    difference = c(rep(0.0001, 8), NA, NA),
    tolerance = c(rep(0.001, 8), NA, NA),
    verdict = rep("pass", 10),
    stringsAsFactors = FALSE
  )
  utils::write.csv(summary, path, row.names = FALSE)

  parsed <- hs_read_blupf90_multivariate_summary(path)
  expect_s3_class(parsed, "hs_blupf90_multivariate_summary")
  expect_equal(parsed$quantity[[1]], "G[1,1]")
  expect_true(is.numeric(parsed$estimate))

  verdict <- hs_validate_blupf90_multivariate_summary(parsed)
  expect_s3_class(verdict, "hs_blupf90_multivariate_summary_validation")
  expect_true(verdict$ok)
  expect_equal(verdict$n_failed, 0L)
  expect_equal(verdict$n_review, 0L)
})

test_that("BLUPF90 multivariate summary ingester rejects incomplete reports", {
  skip_if_not_installed("withr")
  path <- withr::local_tempfile(fileext = ".csv")
  missing_verdict <- data.frame(
    quantity = "G[1,1]",
    target = 1,
    estimate = 1,
    difference = 0,
    tolerance = 0.001
  )
  utils::write.csv(missing_verdict, path, row.names = FALSE)

  expect_error(
    hs_read_blupf90_multivariate_summary(path),
    "missing required columns"
  )
})

test_that("BLUPF90 multivariate summary validation keeps blockers explicit", {
  skip_if_not_installed("withr")
  path <- withr::local_tempfile(fileext = ".csv")
  summary <- data.frame(
    quantity = c(
      "G[1,2]",
      "G[2,2]",
      "R[1,1]",
      "R[1,2]",
      "R[2,2]",
      "h2 trait 1",
      "h2 trait 2"
    ),
    target = 1,
    estimate = 1,
    difference = 0,
    tolerance = 0.001,
    verdict = c(rep("pass", 6), "review")
  )
  utils::write.csv(summary, path, row.names = FALSE)

  parsed <- hs_read_blupf90_multivariate_summary(path)
  verdict <- hs_validate_blupf90_multivariate_summary(parsed)

  expect_false(verdict$ok)
  expect_true("G[1,1]" %in% verdict$missing_quantities)
  expect_true("h2 trait 2" %in% verdict$review_quantities)
})
