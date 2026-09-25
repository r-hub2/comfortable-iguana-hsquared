#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
run_fit <- "--run" %in% args
out_arg <- grep("^--out=", args, value = TRUE)
out_file <- if (length(out_arg)) sub("^--out=", "", out_arg[[1]]) else ""

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
long <- stats::reshape(
  pheno[c("record", "animal", "x", "trait1", "trait2")],
  varying = c("trait1", "trait2"),
  v.names = "value",
  timevar = "trait",
  times = c("trait1", "trait2"),
  idvar = "record",
  direction = "long"
)
long$trait <- factor(long$trait, levels = c("trait1", "trait2"))
long$animal <- factor(long$animal, levels = ped$animal)
long$.record_index <- match(long$record, pheno$record)
long <- long[with(long, order(trait, .record_index)), ]
row.names(long) <- NULL

cat("Prepared ASReml-R candidate fixture\n")
cat("  records:", nrow(long), "\n")
cat("  traits:", paste(levels(long$trait), collapse = ", "), "\n")
cat("  animals:", length(unique(ped$animal)), "\n")
cat("  run fit:", run_fit, "\n")

if (!run_fit) {
  cat("Dry run only. Re-run with --run on a licensed ASReml-R machine.\n")
  quit(status = 0)
}

if (!requireNamespace("asreml", quietly = TRUE)) {
  stop(
    "ASReml-R is not installed or licensed in this R session.",
    call. = FALSE
  )
}

ped_asreml <- ped
ped_asreml$sire[ped_asreml$sire == "0"] <- NA
ped_asreml$dam[ped_asreml$dam == "0"] <- NA

library(asreml)

# Candidate ASReml-R 4 syntax. Review locally before treating the output as
# evidence; ASReml syntax and object slots can differ by licensed version.
ainv <- asreml::ainverse(ped_asreml)
fit <- asreml::asreml(
  fixed = value ~ trait + trait:x - 1,
  random = ~ us(trait):vm(animal, ainv),
  residual = ~ idv(record):us(trait),
  data = long,
  na.action = asreml::na.method(y = "include", x = "include")
)
fit <- asreml::update.asreml(fit)

if (nzchar(out_file)) {
  saveRDS(
    list(
      fit = fit,
      fixture_dir = fixture_dir,
      session_info = utils::sessionInfo()
    ),
    out_file
  )
  cat("Saved ASReml comparator result:", out_file, "\n")
} else {
  print(summary(fit))
}
