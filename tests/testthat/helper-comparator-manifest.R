hs_comparator_manifest_path <- function() {
  normalizePath(
    testthat::test_path("..", "fixtures", "comparator_targets.toml"),
    mustWork = TRUE
  )
}

hs_comparator_fixture_shas_path <- function() {
  normalizePath(
    testthat::test_path("..", "fixtures", "comparator_fixture_shas.csv"),
    mustWork = TRUE
  )
}

# Minimal TOML reader for the comparator manifest schema (schema_version,
# claim_boundary, and [[target]] tables with string/int/bool keys and string
# arrays). Not a general TOML parser.
hs_read_comparator_manifest <- function(path = hs_comparator_manifest_path()) {
  lines <- readLines(path, warn = FALSE)
  lines <- sub("^[[:space:]]*#.*$", "", lines)
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]

  top <- list()
  targets <- list()
  current <- NULL
  array_key <- NULL
  array_vals <- character()

  flush_array <- function() {
    if (is.null(array_key) || is.null(current)) {
      return()
    }
    current[[array_key]] <<- array_vals
    targets[[length(targets)]] <<- current
    array_key <<- NULL
    array_vals <<- character()
  }

  for (line in lines) {
    if (line == "[[target]]") {
      flush_array()
      current <- list()
      targets[[length(targets) + 1L]] <- current
      next
    }

    if (grepl("^\\[", line)) {
      next
    }

    if (grepl("= *\\[$", line)) {
      flush_array()
      array_key <- sub(" *= *\\[$", "", line)
      next
    }

    if (!is.null(array_key) && line == "]") {
      flush_array()
      next
    }

    if (!is.null(array_key)) {
      val <- trimws(line)
      val <- sub(",$", "", val)
      val <- gsub('^"|"$', "", val)
      array_vals <- c(array_vals, val)
      next
    }

    if (!grepl("=", line)) {
      next
    }

    key <- sub(" *=.*", "", line)
    val <- sub("^[^=]*= *", "", line)
    val <- gsub('^"|"$', "", val)
    val <- sub(",$", "", val)

    parsed <- if (val %in% c("true", "false")) {
      as.logical(val)
    } else if (grepl("^[0-9]+$", val)) {
      as.integer(val)
    } else {
      val
    }

    if (is.null(current)) {
      top[[key]] <- parsed
    } else {
      current[[key]] <- parsed
      targets[[length(targets)]] <- current
    }
  }

  flush_array()
  top$target <- targets
  top
}

hs_extract_sha256 <- function(output) {
  output <- tolower(trimws(output))
  unix_candidates <- sub("^([0-9a-f]{64})[[:space:]].*$", "\\1", output)
  unix_candidates <- unix_candidates[grepl("^[0-9a-f]{64}$", unix_candidates)]
  certutil_candidates <- gsub("[[:space:]]", "", output)
  certutil_candidates <- certutil_candidates[
    grepl("^[0-9a-f]{64}$", certutil_candidates)
  ]
  candidates <- unique(c(unix_candidates, certutil_candidates))

  if (length(candidates) != 1L) {
    stop("SHA-256 command did not return exactly one digest", call. = FALSE)
  }

  candidates[[1L]]
}

hs_sha256_file <- function(
  path,
  find_command = Sys.which,
  run_command = system2,
  tools_sha256 = get0("sha256sum", envir = asNamespace("tools"), inherits = FALSE),
  windows_directory = Sys.getenv("WINDIR"),
  system_root = Sys.getenv("SystemRoot")
) {
  if (!is.null(tools_sha256)) {
    digest <- tryCatch(
      unname(tools_sha256(path)),
      error = function(...) NULL
    )
    digest <- tolower(as.character(digest))
    if (length(digest) == 1L && grepl("^[0-9a-f]{64}$", digest)) return(digest)
  }

  commands <- list(
    list(command = "shasum", args = c("-a", "256", path)),
    list(command = "sha256sum", args = path),
    list(command = "certutil", args = c("-hashfile", path, "SHA256"))
  )

  windows_directories <- unique(c(windows_directory, system_root))
  windows_directories <- windows_directories[nzchar(windows_directories)]

  for (directory in windows_directories) {
    commands[[length(commands) + 1L]] <- list(
      command = file.path(directory, "System32", "certutil.exe"),
      args = c("-hashfile", path, "SHA256"),
      system_path = TRUE
    )
  }

  for (candidate in commands) {
    command <- if (isTRUE(candidate$system_path)) {
      candidate$command
    } else {
      find_command(candidate$command)
    }
    if (!nzchar(command)) next
    out <- tryCatch(
      suppressWarnings(run_command(command, candidate$args, stdout = TRUE, stderr = TRUE)),
      error = function(...) NULL
    )
    if (is.null(out) || (!is.null(attr(out, "status")) && attr(out, "status") != 0L)) next
    digest <- tryCatch(hs_extract_sha256(out), error = function(...) NULL)
    if (!is.null(digest)) return(digest)
  }

  stop("no usable SHA-256 command is available (tried shasum, sha256sum, certutil)", call. = FALSE)
}

hs_sha256_canonical_csv <- function(path) {
  bytes <- readBin(path, what = "raw", n = file.info(path)$size)
  canonical <- raw()
  i <- 1L

  while (i <= length(bytes)) {
    if (
      identical(bytes[[i]], as.raw(13L)) &&
        i < length(bytes) &&
        identical(bytes[[i + 1L]], as.raw(10L))
    ) {
      canonical <- c(canonical, as.raw(10L))
      i <- i + 2L
    } else {
      canonical <- c(canonical, bytes[[i]])
      i <- i + 1L
    }
  }

  normalized <- tempfile(fileext = ".csv")
  on.exit(unlink(normalized), add = TRUE)
  writeBin(canonical, normalized)
  hs_sha256_file(normalized)
}
