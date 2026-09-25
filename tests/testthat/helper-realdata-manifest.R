hs_realdata_manifest_path <- function() {
  file.path(
    testthat::test_path("..", ".."),
    "docs",
    "design",
    "real-data-validation-manifest.toml"
  )
}

# `docs/` is .Rbuildignore'd, so the manifest is absent under `R CMD check`,
# which installs from the tarball. These are source-tree contract tests.
hs_skip_without_realdata_manifest <- function() {
  testthat::skip_if_not(
    file.exists(hs_realdata_manifest_path()),
    "docs/design/real-data-validation-manifest.toml is not in the build tarball"
  )
}

hs_parse_toml_scalar <- function(val) {
  val <- trimws(val)
  val <- gsub('^"|"$', "", val)
  val <- sub(",$", "", val)
  if (identical(val, "[]")) {
    return(character(0))
  }
  if (grepl("^\\[.*\\]$", val)) {
    inner <- sub("^\\[(.*)\\]$", "\\1", val)
    if (!nzchar(inner)) {
      return(character(0))
    }
    parts <- strsplit(inner, ",\\s*")[[1]]
    return(vapply(
      trimws(parts),
      function(x) gsub('^"|"$', "", x),
      character(1)
    ))
  }
  if (val %in% c("true", "false")) {
    return(as.logical(val))
  }
  if (grepl("^[0-9]+$", val)) {
    return(as.integer(val))
  }
  val
}

# Minimal TOML reader for the real-data validation manifest schema.
hs_read_realdata_manifest <- function(path = hs_realdata_manifest_path()) {
  lines <- readLines(path, warn = FALSE)
  lines <- sub("^[[:space:]]*#.*$", "", lines)
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]

  result <- list(
    tier_summary = list(),
    darwin_review = list(questions = list(), checklist = list()),
    arc = list()
  )
  section <- "root"
  table <- NULL
  array_key <- NULL
  array_vals <- character()

  commit_table <- function() {
    invisible(NULL)
  }

  update_table <- function(key, value) {
    table[[key]] <<- value
    if (section == "arc") {
      result$arc[[length(result$arc)]] <<- table
    } else if (grepl("^darwin_review\\.", section)) {
      subkey <- sub("^darwin_review\\.", "", section)
      result$darwin_review[[subkey]][[length(result$darwin_review[[
        subkey
      ]])]] <<- table
    }
    invisible(NULL)
  }

  flush_array <- function() {
    if (is.null(array_key)) {
      return(invisible(NULL))
    }
    if (section == "darwin_review") {
      result$darwin_review[[array_key]] <<- array_vals
    } else {
      update_table(array_key, array_vals)
    }
    array_key <<- NULL
    array_vals <<- character()
    invisible(NULL)
  }

  for (line in lines) {
    if (line == "[[arc]]") {
      flush_array()
      commit_table()
      section <- "arc"
      table <- list()
      result$arc[[length(result$arc) + 1L]] <- table
      next
    }

    if (grepl("^\\[\\[darwin_review\\.", line)) {
      flush_array()
      commit_table()
      section <- sub("^\\[\\[(.*)\\]\\]$", "\\1", line)
      table <- list()
      subkey <- sub("^darwin_review\\.", "", section)
      result$darwin_review[[subkey]][[
        length(result$darwin_review[[subkey]]) + 1L
      ]] <- table
      next
    }

    if (line == "[tier_summary]") {
      flush_array()
      commit_table()
      section <- "tier_summary"
      table <- NULL
      next
    }

    if (line == "[darwin_review]") {
      flush_array()
      commit_table()
      section <- "darwin_review"
      table <- NULL
      next
    }

    if (grepl("^\\[", line)) {
      flush_array()
      commit_table()
      section <- "root"
      table <- NULL
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
    parsed <- hs_parse_toml_scalar(val)

    if (section == "tier_summary") {
      result$tier_summary[[key]] <- parsed
    } else if (section == "darwin_review") {
      result$darwin_review[[key]] <- parsed
    } else if (!is.null(table)) {
      update_table(key, parsed)
    } else {
      result[[key]] <- parsed
    }
  }

  flush_array()
  commit_table()
  result
}

hs_realdata_arc_ids <- function(manifest = hs_read_realdata_manifest()) {
  vapply(manifest$arc, `[[`, character(1), "id")
}

hs_realdata_arcs_by_tier <- function(
  tier,
  manifest = hs_read_realdata_manifest()
) {
  arcs <- manifest$arc
  arcs[vapply(arcs, function(a) identical(a$tier, tier), logical(1))]
}
