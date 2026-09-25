# emit_payload_v2_fixtures.R
# P0.5 payload-v2 round-trip parity: R-side fixture emitter.
#
# Builds three bridge payloads via hs_build_bridge_payload(), serializes each
# to a portable JSON file that the Julia parity testset can read directly:
#   (a) v0.1 single-pedigree animal model → dispatch :animal
#   (b) animal + common_env() → dispatch :two_effect
#   (c) animal + permanent()  → dispatch :two_effect
#
# Run from the project root:
#   Rscript tests/fixtures/emit_payload_v2_fixtures.R [OUTPUT_DIR]
#
# The OUTPUT_DIR defaults to
#   ../../HSquared.jl/test/fixtures/payload_v2
# relative to this file's directory (works if both repos are siblings).
#
# CONTRACT: no covered-status change; public_covered_count unchanged (5).
# This script is purely a serialization helper; it contains no Julia calls.

# ---- setup ----------------------------------------------------------------- #

suppressMessages({
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("jsonlite required — install.packages('jsonlite')")
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Matrix required")
  }
})

# Resolve hsquared root: detect from script location or working directory.
cmd_args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("--file=", cmd_args_full, value = TRUE)
if (length(file_arg) > 0L) {
  script_path <- normalizePath(sub("--file=", "", file_arg[[1L]]), mustWork = FALSE)
  this_dir <- dirname(script_path)
  hsq_root <- normalizePath(file.path(this_dir, "../.."), mustWork = FALSE)
} else {
  hsq_root <- normalizePath(".", mustWork = FALSE)
}
if (!file.exists(file.path(hsq_root, "DESCRIPTION"))) {
  # fallback: try working directory
  hsq_root <- normalizePath(".", mustWork = FALSE)
}
suppressMessages(devtools::load_all(hsq_root, quiet = TRUE))
cat("Loaded hsquared from:", hsq_root, "\n")

# Output directory (default: sibling HSquared.jl/test/fixtures/payload_v2).
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1L) {
  out_dir <- normalizePath(args[[1L]], mustWork = FALSE)
} else {
  out_dir <- normalizePath(
    file.path(hsq_root, "..", "HSquared.jl", "test", "fixtures", "payload_v2"),
    mustWork = FALSE
  )
}
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
  cat("Created output directory:", out_dir, "\n")
}
cat("Writing fixtures to:", out_dir, "\n")

# ---- helper: serialize one payload ----------------------------------------- #

#' Serialize a sparse dgCMatrix to a list with triplet + dims.
sparse_to_triplet <- function(M) {
  M <- as(M, "CsparseMatrix")
  M <- as(M, "TsparseMatrix")
  list(
    i    = as.integer(M@i + 1L),   # 1-based row indices for Julia
    j    = as.integer(M@j + 1L),   # 1-based col indices for Julia
    v    = as.numeric(M@x),
    nrow = nrow(M),
    ncol = ncol(M)
  )
}

#' Serialize a dense matrix as flat row-major array with explicit dims.
#' This avoids jsonlite auto_unbox collapsing single-column matrices to scalars.
dense_to_rowmajor <- function(M) {
  M <- unname(as.matrix(M))
  list(
    data = as.numeric(t(M)),   # row-major: row1col1, row1col2, ..., row2col1, ...
    nrow = nrow(M),
    ncol = ncol(M)
  )
}

#' Serialize a pedigree sub-list (id, sire, dam + integer indexes).
pedigree_to_json <- function(ped_list) {
  list(
    id             = as.character(ped_list$id),
    sire           = as.character(ped_list$sire),
    dam            = as.character(ped_list$dam),
    sire_index     = as.integer(ped_list$sire_index),
    dam_index      = as.integer(ped_list$dam_index),
    original_order = as.integer(ped_list$original_order)
  )
}

#' Serialize one block from random_effects.
block_to_json <- function(blk) {
  out <- list(
    name          = blk$name,
    type          = blk$type,
    relmat_status = blk$relmat_status,
    relmat_inverse = NULL,          # always NULL: Julia builds it
    ids            = as.character(blk$ids)
  )
  # Z is always a sparse matrix in the emitter.
  out$Z <- sparse_to_triplet(blk$Z)
  # pedigree sub-list (may be NULL for iid blocks).
  if (!is.null(blk$pedigree)) {
    out$pedigree <- pedigree_to_json(blk$pedigree)
  } else {
    out$pedigree <- NULL
  }
  out
}

#' Serialize a full payload to a JSON-ready list.
payload_to_json_list <- function(payload) {
  # Build the blocks list.
  re_json <- lapply(payload$random_effects, block_to_json)

  list(
    payload_version = as.integer(payload$payload_version),
    y               = as.numeric(payload$y),
    X               = dense_to_rowmajor(payload$X),
    method          = as.character(payload$method),
    family          = as.character(payload$family),
    random_effects  = re_json,
    # --- metadata for human readability + Julia cross-check ---
    metadata = list(
      n_obs             = length(payload$y),
      n_re_blocks       = length(payload$random_effects),
      fixed_colnames    = as.character(payload$metadata$fixed_colnames),
      ainv_status       = as.character(payload$metadata$ainv_status)
    )
  )
}

# ---- shared pedigree + data ------------------------------------------------- #

ped_abcd <- data.frame(
  id   = c("a", "b", "c", "d"),
  sire = c(NA,  NA,  "a", "a"),
  dam  = c(NA,  NA,  "b", "c"),
  stringsAsFactors = FALSE
)

# ============================================================================
# Fixture (a): v0.1 single-pedigree animal model
#   Formula: y ~ sex + animal(1 | id, pedigree = ped)
#   Dispatch (Julia): :animal
# ============================================================================

dat_a <- data.frame(
  y   = c(1.0, 2.0, 3.0),
  sex = c("f", "m", "f"),
  id  = c("a", "c", "d"),
  stringsAsFactors = FALSE
)

spec_a <- hsquared:::hs_build_model_spec(
  y ~ sex + animal(1 | id, pedigree = ped_abcd),
  data   = dat_a,
  family = stats::gaussian(),
  REML   = TRUE
)
payload_a <- hsquared:::hs_build_bridge_payload(spec_a)

stopifnot(payload_a$payload_version == 2L)
stopifnot(length(payload_a$random_effects) == 1L)
stopifnot(payload_a$random_effects[[1L]]$relmat_status == "build_in_julia")

json_a <- payload_to_json_list(payload_a)
out_path_a <- file.path(out_dir, "fixture_a_single_animal.json")
jsonlite::write_json(json_a, out_path_a, auto_unbox = TRUE, null = "null", digits = 15)
cat("Written fixture (a):", out_path_a, "\n")

# ============================================================================
# Fixture (b): animal + common_env()  → pedigree block + iid block
#   Formula: y ~ animal(1 | id, pedigree = ped) + common_env(1 | litter)
#   Dispatch (Julia): :two_effect
# ============================================================================

# 4 records, 4 animals, 2 litter groups.
dat_b <- data.frame(
  y      = c(14.0, 13.0, 12.1, 8.9),
  id     = c("a",  "b",  "c",  "d"),
  litter = c("L1", "L1", "L2", "L2"),
  stringsAsFactors = FALSE
)

spec_b <- hsquared:::hs_build_model_spec(
  y ~ animal(1 | id, pedigree = ped_abcd) + common_env(1 | litter),
  data   = dat_b,
  family = stats::gaussian(),
  REML   = TRUE
)
payload_b <- hsquared:::hs_build_bridge_payload(spec_b)

stopifnot(payload_b$payload_version == 2L)
stopifnot(length(payload_b$random_effects) == 2L)
stopifnot(payload_b$random_effects[[1L]]$name == "animal")
stopifnot(payload_b$random_effects[[2L]]$name == "common_env")
stopifnot(payload_b$random_effects[[2L]]$type == "iid")
stopifnot(payload_b$random_effects[[2L]]$relmat_status == "identity")

json_b <- payload_to_json_list(payload_b)
out_path_b <- file.path(out_dir, "fixture_b_animal_common_env.json")
jsonlite::write_json(json_b, out_path_b, auto_unbox = TRUE, null = "null", digits = 15)
cat("Written fixture (b):", out_path_b, "\n")

# ============================================================================
# Fixture (c): animal + permanent()  → pedigree block + iid block
#   Formula: y ~ animal(1 | id, pedigree = ped) + permanent(1 | id)
#   Dispatch (Julia): :two_effect (permanent is iid over the observed subset)
# ============================================================================

# 6 records, 3 animals with 2 records each.
ped_abc <- data.frame(
  id   = c("a", "b", "c"),
  sire = c(NA,  NA,  "a"),
  dam  = c(NA,  NA,  "b"),
  stringsAsFactors = FALSE
)
dat_c <- data.frame(
  y  = c(10.0, 11.0, 9.0, 12.0, 8.5, 9.5),
  id = c("a",  "a",  "b", "b",  "c", "c"),
  stringsAsFactors = FALSE
)

spec_c <- hsquared:::hs_build_model_spec(
  y ~ animal(1 | id, pedigree = ped_abc) + permanent(1 | id),
  data   = dat_c,
  family = stats::gaussian(),
  REML   = TRUE
)
payload_c <- hsquared:::hs_build_bridge_payload(spec_c)

stopifnot(payload_c$payload_version == 2L)
stopifnot(length(payload_c$random_effects) == 2L)
stopifnot(payload_c$random_effects[[1L]]$name == "animal")
stopifnot(payload_c$random_effects[[2L]]$name == "permanent")
stopifnot(payload_c$random_effects[[2L]]$type == "iid")
stopifnot(payload_c$random_effects[[2L]]$relmat_status == "identity")
# permanent ids must be the observed subset (≤ pedigree size)
stopifnot(length(payload_c$random_effects[[2L]]$ids) <= length(payload_c$ids))

json_c <- payload_to_json_list(payload_c)
out_path_c <- file.path(out_dir, "fixture_c_animal_permanent.json")
jsonlite::write_json(json_c, out_path_c, auto_unbox = TRUE, null = "null", digits = 15)
cat("Written fixture (c):", out_path_c, "\n")

# ---- summary --------------------------------------------------------------- #

cat("\nAll fixtures written successfully.\n")
cat("(a) dispatch expected: :animal     |", out_path_a, "\n")
cat("(b) dispatch expected: :two_effect |", out_path_b, "\n")
cat("(c) dispatch expected: :two_effect |", out_path_c, "\n")
