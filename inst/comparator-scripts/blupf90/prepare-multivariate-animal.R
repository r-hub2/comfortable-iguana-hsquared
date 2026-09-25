#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
write_arg <- grep("^--write=", args, value = TRUE)
write_dir <- if (length(write_arg)) sub("^--write=", "", write_arg[[1]]) else ""

repo <- Sys.getenv("HSQUARED_REPO", unset = getwd())
repo <- normalizePath(repo, mustWork = FALSE)
fixture_dir <- file.path(
  repo,
  "tests",
  "testthat",
  "fixtures",
  "phase4_multitrait_parity"
)

if (!dir.exists(fixture_dir)) {
  stop(
    "Could not find the Phase 4 fixture. Run from the hsquared repo root or ",
    "set HSQUARED_REPO=/path/to/hsquared.",
    call. = FALSE
  )
}

read_fixture <- function(file) {
  utils::read.csv(
    file.path(fixture_dir, file),
    stringsAsFactors = FALSE,
    na.strings = c("", "NA")
  )
}

ped <- read_fixture("pedigree.csv")
pheno <- read_fixture("phenotypes.csv")
id_map <- stats::setNames(seq_len(nrow(ped)), ped$animal)

ped_out <- data.frame(
  animal = unname(id_map[ped$animal]),
  sire = ifelse(ped$sire == "0", 0L, unname(id_map[ped$sire])),
  dam = ifelse(ped$dam == "0", 0L, unname(id_map[ped$dam]))
)
dat_out <- data.frame(
  trait1 = pheno$trait1,
  trait2 = pheno$trait2,
  intercept = 1L,
  x = pheno$x,
  animal = unname(id_map[pheno$animal])
)

cat("Prepared BLUPF90-family candidate fixture\n")
cat("  data rows:", nrow(dat_out), "\n")
cat("  pedigree rows:", nrow(ped_out), "\n")
cat(
  "  write dir:",
  if (nzchar(write_dir)) write_dir else "<dry-run only>",
  "\n"
)

if (!nzchar(write_dir)) {
  cat(
    "Dry run only. Re-run with --write=/path/to/output to create flat files.\n"
  )
  quit(status = 0)
}

dir.create(write_dir, recursive = TRUE, showWarnings = FALSE)
data_file <- file.path(write_dir, "multivariate-animal.dat")
ped_file <- file.path(write_dir, "multivariate-animal.ped")
renf90_file <- file.path(write_dir, "multivariate-animal.renf90")
par_file <- file.path(write_dir, "multivariate-animal.par")

utils::write.table(
  dat_out,
  data_file,
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)
utils::write.table(
  ped_out,
  ped_file,
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

template_candidates <- c(
  file.path(repo, "inst", "comparator-scripts", "blupf90"),
  file.path(repo, "hsquared", "comparator-scripts", "blupf90"),
  file.path(repo, "comparator-scripts", "blupf90"),
  system.file("comparator-scripts", "blupf90", package = "hsquared")
)
template_candidates <- template_candidates[
  nzchar(template_candidates) & dir.exists(template_candidates)
]
if (!length(template_candidates)) {
  stop(
    "Could not find the BLUPF90 comparator templates. Run from the hsquared ",
    "repo root or set HSQUARED_REPO=/path/to/hsquared.",
    call. = FALSE
  )
}
template_dir <- template_candidates[[1]]
copy_template <- function(template, destination) {
  lines <- readLines(file.path(template_dir, template), warn = FALSE)
  lines <- gsub("__DATAFILE__", basename(data_file), lines, fixed = TRUE)
  lines <- gsub("__PEDFILE__", basename(ped_file), lines, fixed = TRUE)
  writeLines(lines, destination, useBytes = TRUE)
}
copy_template("multivariate-animal.renf90", renf90_file)
copy_template("multivariate-animal.par", par_file)

writeLines(
  c(
    "Manual BLUPF90-family comparator bundle.",
    "",
    "Generated from hsquared tests/testthat/fixtures/phase4_multitrait_parity.",
    "Columns in multivariate-animal.dat:",
    "1 trait1",
    "2 trait2",
    "3 intercept",
    "4 x",
    "5 animal",
    "",
    "Run only on a machine with renumf90 and airemlf90 or blupf90+ installed.",
    "Record program versions, convergence, outputs, and scale mapping before",
    "using any result as validation evidence."
  ),
  file.path(write_dir, "README.txt"),
  useBytes = TRUE
)

cat("Wrote BLUPF90-family comparator files to:", write_dir, "\n")
