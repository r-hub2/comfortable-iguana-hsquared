#' Inspect a parsed hsquared model specification
#'
#' `r lifecycle::badge("experimental")`
#'
#' `model_spec()` validates the narrow v0.1 animal-model grammar and returns
#' the R-side model specification that would be used to build the Julia bridge
#' payload. It is a preview and diagnostics helper: it does not fit a model.
#'
#' @param formula A model formula. The current implemented grammar is
#'   `y ~ fixed + animal(1 | id, pedigree = ped)`, or
#'   `y ~ fixed + animal(1 | id)` when `data` is an [hs_data()] object with a
#'   pedigree component.
#' @param data A data frame containing model variables, or an [hs_data()]
#'   object whose `phenotypes` component contains the model variables. When
#'   `data` is an `hs_data` object, formula arguments such as
#'   `pedigree = pedigree` can refer to named components in the bundle, and
#'   `animal(1 | id)` uses the bundle pedigree by default.
#' @param family A response family. The current parser accepts only
#'   `gaussian()`.
#' @param REML Logical; whether the target Gaussian animal model would use
#'   REML.
#' @param ... Reserved for future arguments.
#'
#' @return An `"hs_model_spec"` object.
#' @export
model_spec <- function(
  formula,
  data,
  family = stats::gaussian(),
  REML = TRUE,
  ...
) {
  if (missing(formula)) {
    stop("`formula` is required.", call. = FALSE)
  }
  if (missing(data)) {
    stop("`data` is required.", call. = FALSE)
  }
  dots <- list(...)
  if (length(dots) > 0L) {
    stop("`...` is reserved for future `model_spec()` options.", call. = FALSE)
  }

  spec <- hs_build_model_spec(
    formula = formula,
    data = data,
    family = family,
    REML = REML
  )
  relinv_type <- intersect(names(spec$random), c("genomic", "single_step"))
  if (length(relinv_type) > 0L) {
    hs_abort_unsupported_syntax(
      "`model_spec()` previews the pedigree animal-model grammar only. The ",
      relinv_type[[1L]],
      " model is experimental and opt-in; fit it with ",
      "`hsquared(..., control = hs_control(engine = \"julia\", ",
      "engine_control = list(target = \"",
      relinv_type[[1L]],
      "\")))`.",
      call. = FALSE
    )
  }
  payload <- hs_build_bridge_payload(spec)
  hs_new_model_spec(spec, payload)
}

hs_new_model_spec <- function(spec, payload) {
  animal <- spec$random$animal
  pedigree <- animal$pedigree
  sire <- pedigree$parent_index$sire
  dam <- pedigree$parent_index$dam
  founders <- sire == 0L & dam == 0L

  structure(
    list(
      formula = spec$formula,
      method = spec$method,
      family = spec$family,
      response = spec$response$name,
      dimensions = list(
        observations = length(payload$y),
        fixed_columns = ncol(payload$X),
        animal_ids = length(payload$ids),
        random_design_nonzeros = length(payload$Z@x),
        pedigree_founders = sum(founders)
      ),
      fixed = list(
        terms = spec$fixed$terms,
        columns = payload$metadata$fixed_colnames,
        contrasts = spec$fixed$contrasts
      ),
      animal = list(
        term = animal$term,
        group = animal$group,
        relationship = animal$relationship,
        covariance = animal$covariance,
        pedigree_source = animal$pedigree_source,
        ids = payload$ids,
        observed_ids = payload$metadata$observed_ids,
        observed_id_index = payload$metadata$observed_id_index,
        pedigree = pedigree$data,
        parent_index = pedigree$parent_index
      ),
      julia = list(
        ainv_status = payload$metadata$ainv_status,
        ainv_target = payload$metadata$ainv_target,
        spec_target = payload$metadata$julia_spec_target,
        fit_target = payload$metadata$julia_fit_target
      ),
      payload = payload
    ),
    class = "hs_model_spec"
  )
}

#' @export
print.hs_model_spec <- function(x, ...) {
  dims <- x$dimensions
  family <- paste0(x$family$family, "(", x$family$link, ")")

  cat("<hs_model_spec>\n")
  cat("  response: ", x$response, "\n", sep = "")
  cat("  family: ", family, "\n", sep = "")
  cat("  method: ", x$method, "\n", sep = "")
  cat(
    "  observations: ",
    dims$observations,
    "; fixed columns: ",
    dims$fixed_columns,
    "; animal IDs: ",
    dims$animal_ids,
    "\n",
    sep = ""
  )
  cat(
    "  animal term: ",
    x$animal$term,
    "; group: ",
    x$animal$group,
    "\n",
    sep = ""
  )
  cat(
    "  random design: sparse Z with ",
    dims$random_design_nonzeros,
    " nonzero entries\n",
    sep = ""
  )
  cat("  Ainv: ", x$julia$ainv_status, "\n", sep = "")
  cat("  fitting: not run by `model_spec()`\n", sep = "")
  invisible(x)
}

#' @export
summary.hs_model_spec <- function(object, ...) {
  dims <- object$dimensions
  structure(
    list(
      model = data.frame(
        response = object$response,
        family = object$family$family,
        link = object$family$link,
        method = object$method,
        observations = dims$observations,
        fixed_columns = dims$fixed_columns,
        animal_ids = dims$animal_ids,
        random_design_nonzeros = dims$random_design_nonzeros,
        pedigree_founders = dims$pedigree_founders,
        stringsAsFactors = FALSE
      ),
      fixed = data.frame(
        column = object$fixed$columns,
        stringsAsFactors = FALSE
      ),
      animal = data.frame(
        term = object$animal$term,
        group = object$animal$group,
        relationship = object$animal$relationship,
        covariance = object$animal$covariance,
        pedigree_source = object$animal$pedigree_source,
        ainv_status = object$julia$ainv_status,
        stringsAsFactors = FALSE
      ),
      julia = data.frame(
        stage = c("Ainv", "model specification", "fit"),
        target = c(
          object$julia$ainv_target,
          object$julia$spec_target,
          object$julia$fit_target
        ),
        stringsAsFactors = FALSE
      )
    ),
    class = "summary_hs_model_spec"
  )
}

#' @export
print.summary_hs_model_spec <- function(x, ...) {
  cat("<summary_hs_model_spec>\n")
  print.data.frame(x$model, row.names = FALSE)
  invisible(x)
}
