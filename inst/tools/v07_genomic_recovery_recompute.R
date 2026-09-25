#!/usr/bin/env Rscript

# Independent R recomputation gate for the v0.7 genomic public-activation
# recovery campaign. This file deliberately uses base R only and does not call
# hsquared or the Julia summary implementation.

v07_result_columns <- c(
  "tier", "cell_id", "seed", "n", "m", "truth_sigma_g2",
  "truth_sigma_e2", "truth_ratio", "ridge", "estimate_sigma_g2",
  "estimate_sigma_e2", "estimate_ratio", "converged", "iterations",
  "objective", "gradient_norm", "runtime_seconds", "peak_rss_mb",
  "marker_hash", "id_hash", "kernel_hash", "error_class"
)

v07_manifest_columns <- c(
  "tier", "cell_id", "seed", "n", "m", "truth_sigma_g2",
  "truth_sigma_e2", "truth_ratio", "ridge", "regime"
)

v07_summary_columns <- c(
  "tier", "cell_id", "n_expected", "n_attempted", "n_converged",
  "n_bias_rows", "convergence_rate", "wilson_lower", "wilson_upper",
  "target", "truth", "mean_estimate", "bias", "mcse",
  "pilot_sd_upper", "bias_ci_lower", "bias_ci_upper", "margin", "target_pass",
  "required_n_raw", "required_n", "cell_status", "failure_classes"
)

v07_expected_cells <- as.vector(outer(
  c("n120_m600", "n300_m150", "n300_m1000"),
  c("r020", "r050", "r080"), paste, sep = "_"
))

v07_abort <- function(...) stop(sprintf(...), call. = FALSE)

v07_read_tsv <- function(path, columns) {
  if (!file.exists(path)) v07_abort("missing file: %s", path)
  out <- utils::read.delim(
    path, header = TRUE, sep = "\t", quote = "", comment.char = "",
    stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("NA", "NaN")
  )
  if (!identical(names(out), columns)) {
    v07_abort("schema drift in %s", path)
  }
  out
}

v07_read_raw <- function(out_dir, tier) {
  root <- file.path(out_dir, "raw", tier)
  if (!dir.exists(root)) v07_abort("missing raw result directory: %s", root)
  paths <- sort(list.files(root, pattern = "\\.tsv$", recursive = TRUE, full.names = TRUE))
  if (!length(paths)) v07_abort("no raw result files under %s", root)
  rows <- lapply(paths, v07_read_tsv, columns = v07_result_columns)
  if (any(vapply(rows, nrow, integer(1)) != 1L)) {
    v07_abort("each raw result file must contain exactly one row")
  }
  out <- do.call(rbind, rows)
  attr(out, "raw_paths") <- normalizePath(paths, winslash = "/", mustWork = TRUE)
  rownames(out) <- NULL
  out
}

v07_as_numeric <- function(x, field) {
  out <- suppressWarnings(as.numeric(x))
  bad <- is.na(out) & !is.na(x)
  if (any(bad)) v07_abort("field %s contains nonnumeric values", field)
  out
}

v07_normalize_raw <- function(raw) {
  for (field in c(
    "seed", "n", "m", "truth_sigma_g2", "truth_sigma_e2", "truth_ratio", "ridge",
    "estimate_sigma_g2", "estimate_sigma_e2", "estimate_ratio", "iterations",
    "objective", "gradient_norm", "runtime_seconds", "peak_rss_mb"
  )) raw[[field]] <- v07_as_numeric(raw[[field]], field)
  conv <- tolower(as.character(raw$converged))
  if (any(!conv %in% c("true", "false"))) v07_abort("converged must be true or false")
  raw$converged <- conv == "true"
  raw
}

v07_normalize_manifest <- function(manifest) {
  for (field in c(
    "seed", "n", "m", "truth_sigma_g2", "truth_sigma_e2",
    "truth_ratio", "ridge"
  )) manifest[[field]] <- v07_as_numeric(manifest[[field]], field)
  manifest
}

v07_keys <- function(x) paste(x$tier, x$cell_id, format(x$seed, scientific = FALSE), sep = "\r")

v07_assert_unique_keys <- function(x, label) {
  keys <- v07_keys(x)
  if (anyDuplicated(keys)) v07_abort("duplicate tier/cell/seed key in %s", label)
  invisible(TRUE)
}

v07_assert_sha256 <- function(x, field) {
  if (any(!grepl("^[[:xdigit:]]{64}$", x))) {
    v07_abort("%s must contain lowercase or uppercase SHA-256 hex digests", field)
  }
}

v07_validate_inputs <- function(raw, manifest, tier,
                                expected_cells = v07_expected_cells,
                                tolerance = 1e-10) {
  raw <- v07_normalize_raw(raw)
  manifest <- v07_normalize_manifest(manifest)
  v07_assert_unique_keys(raw, "raw results")
  v07_assert_unique_keys(manifest, "manifest")

  if (any(raw$tier != tier) || any(manifest$tier != tier)) {
    v07_abort("tier mismatch: expected %s", tier)
  }
  if (!setequal(unique(manifest$cell_id), expected_cells)) {
    v07_abort("manifest cell labels differ from the frozen cell set")
  }
  if (!identical(sort(v07_keys(raw)), sort(v07_keys(manifest)))) {
    v07_abort("raw attempted seeds do not exactly equal manifest seeds")
  }
  if (any(abs(manifest$ridge - 0.01) > tolerance)) {
    v07_abort("manifest ridge differs from frozen ridge 0.01")
  }

  joined <- match(v07_keys(raw), v07_keys(manifest))
  for (field in c("n", "m", "truth_sigma_g2", "truth_sigma_e2", "truth_ratio", "ridge")) {
    delta <- abs(raw[[field]] - manifest[[field]][joined])
    if (any(!is.finite(delta)) || any(delta > tolerance)) {
      v07_abort("raw/manifest mismatch in %s", field)
    }
  }

  for (field in c("marker_hash", "id_hash", "kernel_hash")) {
    v07_assert_sha256(as.character(raw[[field]]), field)
  }
  # IDs are deterministically id1,...,idn. Therefore every seed with the same n
  # must have exactly the same order fingerprint, even though marker/kernel
  # fingerprints legitimately vary by seed.
  id_hashes_by_n <- split(raw$id_hash, raw$n)
  if (any(vapply(id_hashes_by_n, function(x) length(unique(x)) != 1L, logical(1)))) {
    v07_abort("ID-order fingerprints disagree within an n regime")
  }
  raw
}

v07_validate_disjoint_seeds <- function(pilot_manifest, confirm_manifest) {
  pilot_manifest <- v07_normalize_manifest(pilot_manifest)
  confirm_manifest <- v07_normalize_manifest(confirm_manifest)
  scientific_key <- function(x) {
    paste(x$cell_id, format(x$seed, scientific = FALSE), sep = "\r")
  }
  overlap <- intersect(scientific_key(pilot_manifest), scientific_key(confirm_manifest))
  if (length(overlap)) v07_abort("pilot and confirmation seed blocks overlap")
  invisible(TRUE)
}

v07_wilson <- function(k, n, z = stats::qnorm(0.975)) {
  if (!n) return(c(lower = NA_real_, upper = NA_real_))
  phat <- k / n
  den <- 1 + z^2 / n
  center <- (phat + z^2 / (2 * n)) / den
  half <- z * sqrt(phat * (1 - phat) / n + z^2 / (4 * n^2)) / den
  c(lower = center - half, upper = center + half)
}

v07_failure_classes <- function(x) {
  tab <- table(as.character(x))
  paste(paste0(names(tab), "=", as.integer(tab)), collapse = ";")
}

v07_recompute_summary <- function(raw, manifest, tier,
                                  expected_cells = v07_expected_cells) {
  raw <- v07_validate_inputs(raw, manifest, tier, expected_cells)
  manifest <- v07_normalize_manifest(manifest)
  z <- stats::qnorm(0.975)
  output <- list()

  for (cell in expected_cells) {
    cr <- raw[raw$cell_id == cell, , drop = FALSE]
    cm <- manifest[manifest$cell_id == cell, , drop = FALSE]
    finite <- is.finite(cr$estimate_sigma_g2) &
      is.finite(cr$estimate_sigma_e2) & is.finite(cr$estimate_ratio)
    eligible <- cr$converged & finite
    n_attempted <- nrow(cr)
    n_converged <- sum(eligible)
    rate <- n_converged / n_attempted
    wilson <- v07_wilson(n_converged, n_attempted, z)
    truths <- unique(cm[c("truth_sigma_g2", "truth_sigma_e2", "truth_ratio")])
    if (nrow(truths) != 1L) v07_abort("truth varies within cell %s", cell)
    targets <- list(
      c("sigma_g2", "estimate_sigma_g2", "truth_sigma_g2"),
      c("sigma_e2", "estimate_sigma_e2", "truth_sigma_e2"),
      c("ratio", "estimate_ratio", "truth_ratio")
    )
    stats_rows <- lapply(targets, function(spec) {
      truth <- truths[[spec[[3L]]]]
      margin <- if (spec[[1L]] == "ratio") 0.02 else 0.05 * truth
      values <- cr[[spec[[2L]]]][eligible]
      mean_estimate <- if (length(values)) mean(values) else NA_real_
      bias <- mean_estimate - truth
      sd_estimate <- if (length(values) > 1L) stats::sd(values) else NA_real_
      mcse <- sd_estimate / sqrt(length(values))
      if (tier == "pilot" && is.finite(sd_estimate)) {
        # Frozen preregistration: protect confirmation sizing against downward
        # noise in the 48-seed pilot SD. This is the one-sided 95% upper
        # confidence bound for a Normal standard deviation.
        df <- length(values) - 1L
        sd_upper <- sd_estimate * sqrt(df / stats::qchisq(0.05, df = df))
        required_raw <- ceiling((z * sd_upper / (margin / 2))^2)
        ci <- c(NA_real_, NA_real_)
        target_pass <- FALSE
      } else if (tier == "confirm" && is.finite(mcse)) {
        sd_upper <- NA_real_
        critical <- stats::qt(0.975, df = length(values) - 1L)
        ci <- bias + c(-1, 1) * critical * mcse
        target_pass <- ci[[1L]] > -margin && ci[[2L]] < margin
        required_raw <- 0
      } else {
        sd_upper <- NA_real_
        required_raw <- 0
        ci <- c(NA_real_, NA_real_)
        target_pass <- FALSE
      }
      data.frame(
        target = spec[[1L]], truth = truth, mean_estimate = mean_estimate,
        bias = bias, mcse = mcse, bias_ci_lower = ci[[1L]],
        pilot_sd_upper = sd_upper, bias_ci_upper = ci[[2L]], margin = margin,
        target_pass = target_pass, required_n_raw = required_raw,
        stringsAsFactors = FALSE
      )
    })
    stats_rows <- do.call(rbind, stats_rows)
    raw_max <- max(stats_rows$required_n_raw)
    required_n <- min(max(200, raw_max), 2000)
    cell_status <- if (tier == "pilot" && rate < 0.95) {
      "STOP_LOW_PILOT_CONVERGENCE"
    } else if (tier == "pilot" && max(200, raw_max) > 2000) {
      "PRECISION_BLOCKER"
    } else if (tier == "pilot") {
      "CONFIRMATION_ELIGIBLE"
    } else if (all(stats_rows$target_pass) && rate >= 0.95 && wilson[[1L]] >= 0.90) {
      "PASS"
    } else {
      "FAIL"
    }
    common <- data.frame(
      tier = tier, cell_id = cell, n_expected = nrow(cm),
      n_attempted = n_attempted, n_converged = n_converged,
      n_bias_rows = n_converged, convergence_rate = rate,
      wilson_lower = wilson[[1L]], wilson_upper = wilson[[2L]],
      stringsAsFactors = FALSE
    )
    tail <- data.frame(
      required_n = required_n, cell_status = cell_status,
      failure_classes = v07_failure_classes(cr$error_class),
      stringsAsFactors = FALSE
    )
    output[[cell]] <- cbind(common[rep(1L, 3L), ], stats_rows, tail[rep(1L, 3L), ])
  }
  out <- do.call(rbind, output)
  rownames(out) <- NULL
  out[v07_summary_columns]
}

v07_equal_numeric <- function(x, y, tolerance) {
  same_missing <- is.na(x) == is.na(y)
  same_value <- is.na(x) | (is.finite(x) & is.finite(y) & abs(x - y) <= tolerance)
  all(same_missing & same_value)
}

v07_compare_summary <- function(recomputed, julia_summary, tolerance = 1e-10) {
  if (!identical(names(julia_summary), v07_summary_columns)) {
    v07_abort("Julia summary schema drift")
  }
  key <- function(x) paste(x$tier, x$cell_id, x$target, sep = "\r")
  if (anyDuplicated(key(recomputed)) || anyDuplicated(key(julia_summary))) {
    v07_abort("duplicate tier/cell/target row in summary")
  }
  if (!identical(sort(key(recomputed)), sort(key(julia_summary)))) {
    v07_abort("R and Julia summary keys differ")
  }
  julia_summary <- julia_summary[match(key(recomputed), key(julia_summary)), , drop = FALSE]
  numeric_fields <- c(
    "n_expected", "n_attempted", "n_converged", "n_bias_rows",
    "convergence_rate", "wilson_lower", "wilson_upper", "truth",
    "mean_estimate", "bias", "mcse", "bias_ci_lower", "bias_ci_upper",
    "pilot_sd_upper", "margin", "required_n_raw", "required_n"
  )
  logical_fields <- "target_pass"
  character_fields <- c("tier", "cell_id", "target", "cell_status", "failure_classes")
  for (field in numeric_fields) {
    lhs <- v07_as_numeric(recomputed[[field]], field)
    rhs <- v07_as_numeric(julia_summary[[field]], field)
    if (!v07_equal_numeric(lhs, rhs, tolerance)) v07_abort("summary mismatch in %s", field)
  }
  for (field in logical_fields) {
    lhs <- tolower(as.character(recomputed[[field]]))
    rhs <- tolower(as.character(julia_summary[[field]]))
    if (!identical(lhs, rhs)) v07_abort("summary mismatch in %s", field)
  }
  for (field in character_fields) {
    if (!identical(as.character(recomputed[[field]]), as.character(julia_summary[[field]]))) {
      v07_abort("summary mismatch in %s", field)
    }
  }
  invisible(TRUE)
}

# A SHA-256 command is used because base R exposes MD5 but not SHA-256. The
# semantic gates remain authoritative if no usable platform command is present.
v07_extract_sha256 <- function(output) {
  output <- tolower(trimws(output))
  unix_candidates <- sub("^([0-9a-f]{64})[[:space:]].*$", "\\1", output)
  unix_candidates <- unix_candidates[grepl("^[0-9a-f]{64}$", unix_candidates)]
  certutil_candidates <- gsub("[[:space:]]", "", output)
  certutil_candidates <- certutil_candidates[
    grepl("^[0-9a-f]{64}$", certutil_candidates)
  ]
  candidates <- unique(c(unix_candidates, certutil_candidates))

  if (length(candidates) != 1L) {
    v07_abort("SHA-256 command did not return exactly one digest")
  }

  candidates[[1L]]
}

v07_sha256 <- function(path) {
  commands <- list(
    list(command = "shasum", args = c("-a", "256", path)),
    list(command = "sha256sum", args = path),
    list(command = "certutil", args = c("-hashfile", path, "SHA256"))
  )

  for (candidate in commands) {
    if (!nzchar(Sys.which(candidate$command))) next
    output <- tryCatch(
      suppressWarnings(system2(candidate$command, candidate$args, stdout = TRUE, stderr = TRUE)),
      error = function(...) NULL
    )
    if (is.null(output) || (!is.null(attr(output, "status")) && attr(output, "status") != 0L)) next
    digest <- tryCatch(v07_extract_sha256(output), error = function(...) NULL)
    if (!is.null(digest)) return(digest)
  }

  v07_abort("no usable SHA-256 command is available (tried shasum, sha256sum, certutil)")
}

v07_has_sha256_command <- function() {
  probe <- tempfile("hsquared-v07-sha256-")
  on.exit(unlink(probe), add = TRUE)
  writeLines("SHA-256 backend probe", probe, useBytes = TRUE)
  !inherits(try(v07_sha256(probe), silent = TRUE), "try-error")
}

v07_raw_lock_path <- function(out_dir, tier) {
  if (!tier %in% c("pilot", "confirm")) v07_abort("tier must be pilot or confirm")
  file.path(out_dir, paste0(tier, "_raw_sha256.tsv"))
}

v07_tier_raw_files <- function(out_dir, tier) {
  root <- file.path(out_dir, "raw", tier)
  sort(list.files(root, pattern = "\\.tsv$", recursive = TRUE, full.names = TRUE))
}

v07_write_raw_lock <- function(out_dir, tier, path = v07_raw_lock_path(out_dir, tier)) {
  if (file.exists(path)) {
    v07_abort("raw campaign checksum lock already exists; refusing to reseal")
  }
  files <- v07_tier_raw_files(out_dir, tier)
  if (!length(files)) v07_abort("no raw %s TSV files to seal", tier)
  root <- paste0(normalizePath(out_dir, winslash = "/", mustWork = TRUE), "/")
  rel <- sub(paste0("^", root), "", normalizePath(files, winslash = "/", mustWork = TRUE))
  lock <- data.frame(path = rel, sha256 = vapply(files, v07_sha256, character(1)), stringsAsFactors = FALSE)
  utils::write.table(lock, path, sep = "\t", row.names = FALSE, quote = FALSE)
  invisible(lock)
}

v07_verify_raw_lock <- function(out_dir, tier, path = v07_raw_lock_path(out_dir, tier)) {
  lock <- v07_read_tsv(path, c("path", "sha256"))
  root <- paste0(normalizePath(out_dir, winslash = "/", mustWork = TRUE), "/")
  current_paths <- sub(
    paste0("^", root), "",
    normalizePath(v07_tier_raw_files(out_dir, tier), winslash = "/", mustWork = TRUE)
  )
  if (!identical(sort(as.character(lock$path)), sort(current_paths))) {
    v07_abort("%s raw checksum lock file set differs from current files", tier)
  }
  current <- vapply(file.path(out_dir, lock$path), v07_sha256, character(1))
  if (!identical(unname(current), as.character(lock$sha256))) {
    v07_abort("raw campaign checksum lock does not match current files")
  }
  invisible(TRUE)
}

v07_run_gate <- function(out_dir, tier, tolerance = 1e-10, verify_lock = TRUE) {
  manifest_name <- if (tier == "pilot") "pilot_manifest.tsv" else "confirmation_manifest.tsv"
  raw <- v07_read_raw(out_dir, tier)
  manifest <- v07_read_tsv(file.path(out_dir, manifest_name), v07_manifest_columns)
  if (verify_lock) v07_verify_raw_lock(out_dir, tier)
  if (file.exists(file.path(out_dir, "pilot_manifest.tsv")) &&
      file.exists(file.path(out_dir, "confirmation_manifest.tsv"))) {
    v07_validate_disjoint_seeds(
      v07_read_tsv(file.path(out_dir, "pilot_manifest.tsv"), v07_manifest_columns),
      v07_read_tsv(file.path(out_dir, "confirmation_manifest.tsv"), v07_manifest_columns)
    )
  }
  recomputed <- v07_recompute_summary(raw, manifest, tier)
  julia <- v07_read_tsv(file.path(out_dir, paste0(tier, "_summary.tsv")), v07_summary_columns)
  v07_compare_summary(recomputed, julia, tolerance)
  invisible(recomputed)
}

v07_option <- function(args, key, default = NULL) {
  prefix <- paste0("--", key, "=")
  hit <- args[startsWith(args, prefix)]
  if (!length(hit)) return(default)
  sub(prefix, "", hit[[1L]], fixed = TRUE)
}

v07_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  out_dir <- v07_option(args, "out-dir")
  tier <- v07_option(args, "tier")
  if (is.null(out_dir) || !tier %in% c("pilot", "confirm")) {
    v07_abort("usage: Rscript tools/v07_genomic_recovery_recompute.R --out-dir=PATH --tier=pilot|confirm [--write-lock=true]")
  }
  write_lock <- identical(tolower(v07_option(args, "write-lock", "false")), "true")
  if (write_lock) v07_write_raw_lock(out_dir, tier)
  result <- v07_run_gate(out_dir, tier)
  message(sprintf("v0.7 %s independent R recomputation: PASS (%d summary rows)", tier, nrow(result)))
}

if (sys.nframe() == 0L) v07_main()
