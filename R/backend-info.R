#' Inspect planned compute backends
#'
#' `r lifecycle::badge("experimental")`
#'
#' `backend_info()` reports which backend names are accepted by
#' [hs_control()] and whether any of them are execution-ready from the R
#' package. In the current package state, backend names are control metadata
#' only: they are selectable but not dispatched.
#'
#' @param control An object created by [hs_control()].
#'
#' @return A data frame of backend status records with class
#'   `"hs_backend_info"`.
#' @export
backend_info <- function(control = hs_control()) {
  if (!inherits(control, "hs_control")) {
    stop("`control` must be created by `hs_control()`.", call. = FALSE)
  }

  backends <- c("cpu", "threads", "cuda", "amdgpu", "metal", "oneapi")
  accelerator <- c("none", "none", "cuda", "amdgpu", "metal", "oneapi")

  out <- data.frame(
    backend = backends,
    accelerator = accelerator,
    requested = hs_backend_requested(backends, control),
    selectable = TRUE,
    execution_available = FALSE,
    status = "planned",
    note = hs_backend_note(backends),
    stringsAsFactors = FALSE
  )
  class(out) <- c("hs_backend_info", class(out))
  attr(out, "control") <- list(
    engine = control$engine,
    backend = control$backend,
    accelerator = control$accelerator,
    precision = control$precision
  )
  out
}

#' @export
print.hs_backend_info <- function(x, ...) {
  control <- attr(x, "control")
  cat("<hs_backend_info>\n")
  cat("  engine: ", control$engine, "\n", sep = "")
  cat("  requested backend: ", control$backend, "\n", sep = "")
  cat("  requested accelerator: ", control$accelerator, "\n", sep = "")
  cat("  execution: metadata only; no backend dispatch yet\n")
  print.data.frame(unclass(x), row.names = FALSE)
  invisible(x)
}

hs_backend_requested <- function(backends, control) {
  backend_match <- backends == control$backend
  accelerator_match <- backends == control$accelerator
  gpu_match <- identical(control$accelerator, "gpu") &
    backends %in% c("cuda", "amdgpu", "metal", "oneapi")
  backend_match | accelerator_match | gpu_match
}

hs_backend_note <- function(backends) {
  notes <- rep(
    "Accepted by `hs_control()`; execution dispatch is planned.",
    length(backends)
  )
  notes[backends == "cpu"] <- paste(
    "Trusted default target, but `backend = \"cpu\"` is not dispatched by R",
    "yet."
  )
  notes[backends == "threads"] <- paste(
    "Planned multi-threaded CPU control; not dispatched by R yet."
  )
  notes[backends == "metal"] <- paste(
    "Planned Apple/Mac accelerator control; not dispatched by R yet."
  )
  notes
}
