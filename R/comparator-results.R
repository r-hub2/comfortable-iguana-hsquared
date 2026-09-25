hs_blupf90_summary_required_columns <- function() {
  c("quantity", "target", "estimate", "difference", "tolerance", "verdict")
}

hs_blupf90_summary_required_quantities <- function() {
  c(
    "G[1,1]",
    "G[1,2]",
    "G[2,2]",
    "R[1,1]",
    "R[1,2]",
    "R[2,2]",
    "h2 trait 1",
    "h2 trait 2"
  )
}

hs_read_blupf90_multivariate_summary <- function(path) {
  if (
    !is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)
  ) {
    stop("`path` must be a single non-empty file path.", call. = FALSE)
  }
  if (!file.exists(path)) {
    stop("BLUPF90 summary file does not exist: ", path, call. = FALSE)
  }

  out <- utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = c("", "NA")
  )
  required_columns <- hs_blupf90_summary_required_columns()
  missing_columns <- setdiff(required_columns, names(out))
  if (length(missing_columns) > 0L) {
    stop(
      "BLUPF90 summary is missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  out <- out[,
    c(required_columns, setdiff(names(out), required_columns)),
    drop = FALSE
  ]
  out$quantity <- trimws(as.character(out$quantity))
  out$verdict <- tolower(trimws(as.character(out$verdict)))

  numeric_columns <- c("target", "estimate", "difference", "tolerance")
  for (column in numeric_columns) {
    out[[column]] <- suppressWarnings(as.numeric(out[[column]]))
  }

  if (any(is.na(out$quantity) | !nzchar(out$quantity))) {
    stop("BLUPF90 summary contains blank `quantity` values.", call. = FALSE)
  }
  if (any(is.na(out$verdict) | !nzchar(out$verdict))) {
    stop("BLUPF90 summary contains blank `verdict` values.", call. = FALSE)
  }

  class(out) <- c("hs_blupf90_multivariate_summary", class(out))
  out
}

hs_validate_blupf90_multivariate_summary <- function(
  x,
  required_quantities = hs_blupf90_summary_required_quantities()
) {
  if (!inherits(x, "hs_blupf90_multivariate_summary")) {
    stop(
      "`x` must be read by `hs_read_blupf90_multivariate_summary()`.",
      call. = FALSE
    )
  }
  if (!is.character(required_quantities) || any(!nzchar(required_quantities))) {
    stop(
      "`required_quantities` must be non-empty quantity labels.",
      call. = FALSE
    )
  }

  accepted_verdicts <- c("pass", "fail", "failed", "review", "unclear")
  bad_verdict_values <- setdiff(unique(x$verdict), accepted_verdicts)
  missing_quantities <- setdiff(required_quantities, x$quantity)
  failed_quantities <- x$quantity[x$verdict %in% c("fail", "failed")]
  review_quantities <- x$quantity[x$verdict %in% c("review", "unclear")]

  out <- list(
    ok = length(missing_quantities) == 0L &&
      length(failed_quantities) == 0L &&
      length(review_quantities) == 0L &&
      length(bad_verdict_values) == 0L,
    missing_quantities = missing_quantities,
    failed_quantities = failed_quantities,
    review_quantities = review_quantities,
    bad_verdict_values = bad_verdict_values,
    n_required = length(required_quantities),
    n_observed = length(unique(x$quantity)),
    n_failed = length(failed_quantities),
    n_review = length(review_quantities)
  )
  class(out) <- c("hs_blupf90_multivariate_summary_validation", class(out))
  out
}
