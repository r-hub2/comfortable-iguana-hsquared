#' Create an hsquared data container
#'
#' `r lifecycle::badge("experimental")`
#'
#' `hs_data()` collects phenotype, pedigree, genotype, marker, expression,
#' annotation, and environment inputs into one checked container. It is a
#' lightweight data-contract object for future genomic, QTL/eQTL, and
#' multi-omics workflows. It does not fit models. The v0.1 parser can use an
#' `hs_data` object directly as `data`, reading model variables from
#' `phenotypes` and making named components such as `pedigree` available to
#' formula terms.
#'
#' @param phenotypes A data frame of phenotypic records.
#' @param pedigree Optional pedigree data frame.
#' @param genotypes Optional genotype matrix or data frame. Matrix row names or
#'   data-frame ID values identify individuals. When `markers` is supplied,
#'   genotype marker column names must match marker-map IDs exactly.
#' @param markers Optional marker map data frame. When supplied, it must contain
#'   marker ID, chromosome, and position columns. Recognized aliases include
#'   `marker`, `snp`, or `id`; `chromosome`, `chr`, or `chrom`; and
#'   `position`, `pos`, `bp`, or `base_pair`.
#' @param expression Optional expression matrix or data frame.
#' @param annotation Optional annotation data frame.
#' @param annotation_id Optional column name used to match `annotation` rows to
#'   expression feature columns. When supplied, the column must exist in
#'   `annotation`.
#' @param environment Optional environment/covariate data frame.
#' @param environment_id Optional column name used to match `environment` rows
#'   to phenotype records. When supplied, the column must exist in both
#'   `phenotypes` and `environment`.
#' @param id Name of the individual ID column in `phenotypes`.
#'
#' @return An `hs_data` object.
#'
#' @details
#' `summary(hs_data(...))` reports ID overlap diagnostics, pedigree diagnostics,
#' and, when expression, genotype, or marker components are supplied,
#' expression-feature, genotype-column, marker-map, and genotype-marker
#' alignment diagnostics. When `annotation_id` is supplied, it reports
#' expression-feature annotation coverage diagnostics. When `environment_id` is
#' supplied, it also reports environment metadata coverage diagnostics.
#' @examples
#' ped <- data.frame(
#'   id = c("a", "b", "c", "d"),
#'   sire = c(NA, NA, "a", "a"),
#'   dam = c(NA, NA, "b", "c")
#' )
#' dat <- data.frame(
#'   y = c(1, 2, 3),
#'   sex = c("f", "m", "f"),
#'   id = c("a", "c", "d")
#' )
#' bundle <- hs_data(phenotypes = dat, pedigree = ped)
#' bundle
#' summary(bundle)
#' @export
hs_data <- function(
  phenotypes,
  pedigree = NULL,
  genotypes = NULL,
  markers = NULL,
  expression = NULL,
  annotation = NULL,
  annotation_id = NULL,
  environment = NULL,
  environment_id = NULL,
  id = "id"
) {
  if (!is.data.frame(phenotypes)) {
    stop("`phenotypes` must be a data frame.", call. = FALSE)
  }
  if (!is.character(id) || length(id) != 1L || is.na(id) || id == "") {
    stop("`id` must be one non-empty column name.", call. = FALSE)
  }
  if (!id %in% names(phenotypes)) {
    stop("`id` column `", id, "` was not found in `phenotypes`.", call. = FALSE)
  }

  phenotype_ids <- hs_checked_ids(phenotypes[[id]], "`phenotypes`")
  pedigree_ids <- hs_data_pedigree_ids(pedigree, phenotype_ids)
  genotype_ids <- hs_optional_component_ids(genotypes, id, "`genotypes`")
  expression_ids <- hs_optional_component_ids(expression, id, "`expression`")

  marker_spec <- hs_validate_marker_map(markers)
  genotype_marker_spec <- hs_validate_genotype_marker_alignment(
    genotypes,
    id,
    marker_spec
  )
  annotation_spec <- hs_validate_annotation(
    annotation,
    expression,
    id,
    annotation_id
  )
  environment_spec <- hs_validate_environment(
    environment,
    phenotypes,
    environment_id
  )

  structure(
    list(
      phenotypes = phenotypes,
      pedigree = pedigree,
      genotypes = genotypes,
      markers = markers,
      marker_spec = marker_spec,
      genotype_marker_spec = genotype_marker_spec,
      expression = expression,
      annotation = annotation,
      annotation_id = annotation_id,
      annotation_spec = annotation_spec,
      environment = environment,
      environment_id = environment_id,
      environment_spec = environment_spec,
      id = id,
      id_map = list(
        phenotype_ids = unique(phenotype_ids),
        pedigree_ids = pedigree_ids,
        genotype_ids = genotype_ids,
        expression_ids = expression_ids,
        phenotypes_without_pedigree = setdiff(
          unique(phenotype_ids),
          pedigree_ids
        ),
        phenotypes_without_genotypes = setdiff(
          unique(phenotype_ids),
          genotype_ids
        ),
        genotypes_without_phenotypes = setdiff(
          genotype_ids,
          unique(phenotype_ids)
        ),
        phenotypes_without_expression = setdiff(
          unique(phenotype_ids),
          expression_ids
        ),
        expression_without_phenotypes = setdiff(
          expression_ids,
          unique(phenotype_ids)
        )
      )
    ),
    class = "hs_data"
  )
}

#' @export
print.hs_data <- function(x, ...) {
  cat("<hs_data>\n")
  cat("  phenotypes: ", nrow(x$phenotypes), " rows\n", sep = "")
  cat("  phenotype IDs: ", length(x$id_map$phenotype_ids), "\n", sep = "")
  hs_print_component_count("pedigree IDs", x$id_map$pedigree_ids)
  hs_print_component_count("genotype IDs", x$id_map$genotype_ids)
  hs_print_component_count("expression IDs", x$id_map$expression_ids)
  invisible(x)
}

#' Inspect hsquared data-container status
#'
#' `r lifecycle::badge("experimental")`
#'
#' `data_status()` gives a direct user-facing view of the checks stored in an
#' [hs_data()] object. It reports component presence, ID overlap diagnostics,
#' pedigree diagnostics, expression-feature diagnostics, genotype-column
#' diagnostics, and marker-map/genotype-marker alignment diagnostics when those
#' inputs are supplied. When `annotation_id` is supplied, it reports
#' expression-feature annotation coverage diagnostics. When `environment_id` is
#' supplied, it also reports environment metadata coverage diagnostics. It does
#' not fit models, build genomic relationship matrices, add eQTL terms, or add
#' environment-effect terms.
#'
#' @param data An [hs_data()] object.
#'
#' @return An `"hs_data_status"` object.
#' @examples
#' ped <- data.frame(
#'   id = c("a", "b", "c", "d"),
#'   sire = c(NA, NA, "a", "a"),
#'   dam = c(NA, NA, "b", "c")
#' )
#' dat <- data.frame(
#'   y = c(1, 2, 3),
#'   sex = c("f", "m", "f"),
#'   id = c("a", "c", "d")
#' )
#' bundle <- hs_data(phenotypes = dat, pedigree = ped)
#' data_status(bundle)
#' @export
data_status <- function(data) {
  UseMethod("data_status")
}

#' @export
data_status.default <- function(data) {
  stop(
    "`data_status()` currently supports `hs_data()` objects.",
    call. = FALSE
  )
}

#' @export
data_status.hs_data <- function(data) {
  out <- summary(data)
  structure(
    list(
      components = out$components,
      id_overlap = out$id_overlap,
      pedigree_status = out$pedigree_status,
      expression_status = out$expression_status,
      genotype_status = out$genotype_status,
      marker_status = out$marker_status,
      annotation_status = out$annotation_status,
      environment_status = out$environment_status
    ),
    class = "hs_data_status"
  )
}

#' @export
print.hs_data_status <- function(x, ...) {
  cat("<hs_data_status>\n")
  cat("  components: ", paste(x$components, collapse = ", "), "\n", sep = "")
  cat("  ID overlap:\n", sep = "")
  print.data.frame(x$id_overlap, row.names = FALSE)
  if (is.null(x$pedigree_status)) {
    cat("  pedigree status: not available\n", sep = "")
  } else {
    cat("  pedigree status:\n", sep = "")
    print.data.frame(x$pedigree_status, row.names = FALSE)
  }
  if (is.null(x$expression_status)) {
    cat("  expression status: not available\n", sep = "")
  } else {
    cat("  expression status:\n", sep = "")
    print.data.frame(x$expression_status, row.names = FALSE)
  }
  if (is.null(x$genotype_status)) {
    cat("  genotype status: not available\n", sep = "")
  } else {
    cat("  genotype status:\n", sep = "")
    print.data.frame(x$genotype_status, row.names = FALSE)
  }
  if (is.null(x$marker_status)) {
    cat("  marker status: not available\n", sep = "")
  } else {
    cat("  marker status:\n", sep = "")
    print.data.frame(x$marker_status, row.names = FALSE)
  }
  if (is.null(x$annotation_status)) {
    cat("  annotation status: not available\n", sep = "")
  } else {
    cat("  annotation status:\n", sep = "")
    print.data.frame(x$annotation_status, row.names = FALSE)
  }
  if (is.null(x$environment_status)) {
    cat("  environment status: not available\n", sep = "")
  } else {
    cat("  environment status:\n", sep = "")
    print.data.frame(x$environment_status, row.names = FALSE)
  }
  invisible(x)
}

#' @export
summary.hs_data <- function(object, ...) {
  structure(
    list(
      components = names(Filter(Negate(is.null), unclass(object)))[
        names(Filter(Negate(is.null), unclass(object))) %in%
          c(
            "phenotypes",
            "pedigree",
            "genotypes",
            "markers",
            "expression",
            "annotation",
            "environment"
          )
      ],
      id_map = object$id_map,
      id_overlap = hs_data_id_overlap(object$id_map),
      pedigree_status = hs_data_pedigree_status(object),
      expression_status = hs_data_expression_status(object),
      genotype_status = hs_data_genotype_status(object),
      marker_status = hs_data_marker_status(object),
      annotation_status = hs_data_annotation_status(object),
      environment_status = hs_data_environment_status(object)
    ),
    class = "summary_hs_data"
  )
}

#' @export
print.summary_hs_data <- function(x, ...) {
  cat("<summary_hs_data>\n")
  cat("  components: ", paste(x$components, collapse = ", "), "\n", sep = "")
  cat("  phenotype IDs: ", length(x$id_map$phenotype_ids), "\n", sep = "")
  cat("  ID overlap:\n", sep = "")
  print.data.frame(x$id_overlap, row.names = FALSE)
  if (!is.null(x$pedigree_status)) {
    cat("  pedigree status:\n", sep = "")
    print.data.frame(x$pedigree_status, row.names = FALSE)
  }
  if (!is.null(x$expression_status)) {
    cat("  expression status:\n", sep = "")
    print.data.frame(x$expression_status, row.names = FALSE)
  }
  if (!is.null(x$genotype_status)) {
    cat("  genotype status:\n", sep = "")
    print.data.frame(x$genotype_status, row.names = FALSE)
  }
  if (!is.null(x$marker_status)) {
    cat("  marker status:\n", sep = "")
    print.data.frame(x$marker_status, row.names = FALSE)
  }
  if (!is.null(x$annotation_status)) {
    cat("  annotation status:\n", sep = "")
    print.data.frame(x$annotation_status, row.names = FALSE)
  }
  if (!is.null(x$environment_status)) {
    cat("  environment status:\n", sep = "")
    print.data.frame(x$environment_status, row.names = FALSE)
  }
  invisible(x)
}

hs_data_id_overlap <- function(id_map) {
  data.frame(
    metric = c(
      "phenotype_ids",
      "pedigree_ids",
      "genotype_ids",
      "expression_ids",
      "phenotypes_without_pedigree",
      "phenotypes_without_genotypes",
      "genotypes_without_phenotypes",
      "phenotypes_without_expression",
      "expression_without_phenotypes"
    ),
    count = c(
      length(id_map$phenotype_ids),
      length(id_map$pedigree_ids),
      length(id_map$genotype_ids),
      length(id_map$expression_ids),
      length(id_map$phenotypes_without_pedigree),
      length(id_map$phenotypes_without_genotypes),
      length(id_map$genotypes_without_phenotypes),
      length(id_map$phenotypes_without_expression),
      length(id_map$expression_without_phenotypes)
    ),
    stringsAsFactors = FALSE
  )
}

hs_data_pedigree_status <- function(object) {
  if (is.null(object$pedigree)) {
    return(NULL)
  }

  pedigree <- object$pedigree
  cols <- hs_pedigree_columns(pedigree)
  ids <- as.character(pedigree[[cols$id]])
  sire <- hs_normalize_parent(pedigree[[cols$sire]])
  dam <- hs_normalize_parent(pedigree[[cols$dam]])
  known_parent_ids <- unique(c(stats::na.omit(sire), stats::na.omit(dam)))

  duplicate_ids <- unique(ids[duplicated(ids)])
  missing_parents <- setdiff(known_parent_ids, unique(ids))
  self_parent <- (!is.na(sire) & sire == ids) | (!is.na(dam) & dam == ids)
  same_known_parent <- !is.na(sire) & !is.na(dam) & sire == dam
  founders <- is.na(sire) & is.na(dam)
  phenotype_ids <- object$id_map$phenotype_ids
  pedigree_ids <- object$id_map$pedigree_ids

  data.frame(
    metric = c(
      "pedigree_rows",
      "pedigree_ids",
      "phenotype_ids_with_pedigree",
      "pedigree_only_ids",
      "founders",
      "nonfounders",
      "known_sire_links",
      "known_dam_links",
      "missing_known_parent_ids",
      "duplicate_pedigree_ids",
      "self_parent_rows",
      "same_known_parent_rows"
    ),
    count = c(
      nrow(pedigree),
      length(pedigree_ids),
      length(intersect(phenotype_ids, pedigree_ids)),
      length(setdiff(pedigree_ids, phenotype_ids)),
      sum(founders),
      sum(!founders),
      sum(!is.na(sire)),
      sum(!is.na(dam)),
      length(missing_parents),
      length(duplicate_ids),
      sum(self_parent),
      sum(same_known_parent)
    ),
    stringsAsFactors = FALSE
  )
}

hs_data_marker_status <- function(object) {
  marker_spec <- object$marker_spec
  genotype_marker_spec <- object$genotype_marker_spec
  genotype_marker_ids <- hs_data_summary_genotype_marker_ids(object)

  if (is.null(marker_spec) && length(genotype_marker_ids) == 0L) {
    return(NULL)
  }

  marker_count <- if (is.null(marker_spec)) {
    0L
  } else {
    length(marker_spec$marker_ids)
  }
  genotype_count <- length(genotype_marker_ids)
  aligned_count <- if (is.null(genotype_marker_spec)) {
    0L
  } else {
    length(genotype_marker_spec$marker_ids)
  }
  chromosome_count <- if (is.null(marker_spec)) {
    NA_integer_
  } else {
    length(unique(marker_spec$chromosome))
  }
  position_min <- if (is.null(marker_spec)) {
    NA_real_
  } else {
    min(marker_spec$position)
  }
  position_max <- if (is.null(marker_spec)) {
    NA_real_
  } else {
    max(marker_spec$position)
  }
  alignment <- hs_data_marker_alignment_status(
    marker_spec,
    genotype_marker_spec,
    genotype_count
  )

  data.frame(
    metric = c(
      "marker_map_markers",
      "genotype_marker_columns",
      "aligned_marker_columns",
      "chromosomes",
      "position_min",
      "position_max",
      "alignment"
    ),
    value = c(
      as.character(marker_count),
      as.character(genotype_count),
      as.character(aligned_count),
      hs_optional_summary_value(chromosome_count),
      hs_optional_summary_value(position_min),
      hs_optional_summary_value(position_max),
      alignment
    ),
    stringsAsFactors = FALSE
  )
}

hs_data_expression_status <- function(object) {
  if (is.null(object$expression)) {
    return(NULL)
  }

  feature_ids <- hs_expression_feature_ids(
    object$expression,
    object$id,
    require_names = FALSE
  )
  has_feature_name <- !is.na(feature_ids) & feature_ids != ""
  named_features <- feature_ids[has_feature_name]
  duplicate_features <- unique(named_features[duplicated(named_features)])

  data.frame(
    metric = c(
      "expression_rows",
      "expression_ids",
      "expression_features",
      "named_expression_features",
      "unnamed_expression_features",
      "duplicate_expression_features",
      "component_type"
    ),
    value = c(
      as.character(hs_expression_row_count(object$expression)),
      as.character(length(object$id_map$expression_ids)),
      as.character(length(feature_ids)),
      as.character(length(named_features)),
      as.character(sum(!has_feature_name)),
      as.character(length(duplicate_features)),
      hs_expression_component_type(object$expression)
    ),
    stringsAsFactors = FALSE
  )
}

hs_data_genotype_status <- function(object) {
  if (is.null(object$genotypes)) {
    return(NULL)
  }

  marker_ids <- hs_data_summary_genotype_marker_ids(object)
  has_marker_name <- !is.na(marker_ids) & marker_ids != ""
  named_markers <- marker_ids[has_marker_name]
  duplicate_markers <- unique(named_markers[duplicated(named_markers)])

  data.frame(
    metric = c(
      "genotype_rows",
      "genotype_ids",
      "genotype_marker_columns",
      "named_genotype_marker_columns",
      "unnamed_genotype_marker_columns",
      "duplicate_genotype_marker_columns",
      "missing_genotype_values",
      "component_type"
    ),
    value = c(
      as.character(hs_genotype_row_count(object$genotypes)),
      as.character(length(object$id_map$genotype_ids)),
      as.character(length(marker_ids)),
      as.character(length(named_markers)),
      as.character(sum(!has_marker_name)),
      as.character(length(duplicate_markers)),
      as.character(hs_genotype_missing_value_count(
        object$genotypes,
        object$id
      )),
      hs_genotype_component_type(object$genotypes)
    ),
    stringsAsFactors = FALSE
  )
}

hs_data_annotation_status <- function(object) {
  if (is.null(object$annotation)) {
    return(NULL)
  }

  if (is.null(object$annotation_spec)) {
    return(data.frame(
      metric = c(
        "annotation_rows",
        "annotation_key",
        "annotation_features",
        "expression_features",
        "expression_features_with_annotation",
        "annotation_only_features",
        "expression_features_without_annotation",
        "duplicate_annotation_features"
      ),
      value = c(
        as.character(nrow(object$annotation)),
        "not_checked_no_annotation_id",
        rep("not_available", 6L)
      ),
      stringsAsFactors = FALSE
    ))
  }

  spec <- object$annotation_spec
  data.frame(
    metric = c(
      "annotation_rows",
      "annotation_key",
      "annotation_features",
      "expression_features",
      "expression_features_with_annotation",
      "annotation_only_features",
      "expression_features_without_annotation",
      "duplicate_annotation_features"
    ),
    value = c(
      as.character(nrow(object$annotation)),
      spec$key,
      as.character(length(spec$annotation_features)),
      as.character(length(spec$expression_features)),
      as.character(length(intersect(
        spec$expression_features,
        spec$annotation_features
      ))),
      as.character(length(spec$annotation_without_expression)),
      as.character(length(spec$expression_without_annotation)),
      as.character(length(spec$duplicate_annotation_features))
    ),
    stringsAsFactors = FALSE
  )
}

hs_data_environment_status <- function(object) {
  if (is.null(object$environment)) {
    return(NULL)
  }

  if (is.null(object$environment_spec)) {
    return(data.frame(
      metric = c(
        "environment_rows",
        "environment_key",
        "environment_ids",
        "phenotype_environment_ids",
        "phenotype_environment_ids_with_metadata",
        "environment_only_ids",
        "phenotype_environment_ids_without_metadata",
        "duplicate_environment_ids"
      ),
      value = c(
        as.character(nrow(object$environment)),
        "not_checked_no_environment_id",
        rep("not_available", 6L)
      ),
      stringsAsFactors = FALSE
    ))
  }

  spec <- object$environment_spec
  data.frame(
    metric = c(
      "environment_rows",
      "environment_key",
      "environment_ids",
      "phenotype_environment_ids",
      "phenotype_environment_ids_with_metadata",
      "environment_only_ids",
      "phenotype_environment_ids_without_metadata",
      "duplicate_environment_ids"
    ),
    value = c(
      as.character(nrow(object$environment)),
      spec$key,
      as.character(length(spec$environment_ids)),
      as.character(length(spec$phenotype_environment_ids)),
      as.character(length(intersect(
        spec$phenotype_environment_ids,
        spec$environment_ids
      ))),
      as.character(length(spec$environment_without_phenotypes)),
      as.character(length(spec$phenotypes_without_environment)),
      as.character(length(spec$duplicate_environment_ids))
    ),
    stringsAsFactors = FALSE
  )
}

hs_data_summary_genotype_marker_ids <- function(object) {
  if (is.null(object$genotypes)) {
    return(character())
  }
  if (is.matrix(object$genotypes) && is.null(colnames(object$genotypes))) {
    return(rep("", ncol(object$genotypes)))
  }
  if (is.data.frame(object$genotypes)) {
    marker_columns <- names(object$genotypes)
    if (object$id %in% marker_columns) {
      marker_columns <- marker_columns[marker_columns != object$id]
    }
    return(marker_columns)
  }
  hs_genotype_marker_ids(object$genotypes, object$id)
}

hs_data_marker_alignment_status <- function(
  marker_spec,
  genotype_marker_spec,
  genotype_count
) {
  if (!is.null(genotype_marker_spec)) {
    return("checked")
  }
  if (!is.null(marker_spec) && genotype_count == 0L) {
    return("not_checked_no_genotypes")
  }
  if (is.null(marker_spec) && genotype_count > 0L) {
    return("not_checked_no_marker_map")
  }
  "not_applicable"
}

hs_optional_summary_value <- function(x) {
  if (length(x) == 0L || is.na(x)) {
    return("not_available")
  }
  as.character(x)
}

hs_checked_ids <- function(x, label) {
  ids <- as.character(x)
  if (any(is.na(ids) | ids == "" | ids == "0")) {
    stop(label, " IDs cannot be missing, empty, or `0`.", call. = FALSE)
  }
  ids
}

hs_data_pedigree_ids <- function(pedigree, phenotype_ids) {
  if (is.null(pedigree)) {
    return(character())
  }
  if (!is.data.frame(pedigree)) {
    stop("`pedigree` must be a data frame when supplied.", call. = FALSE)
  }
  if (ncol(pedigree) < 3L) {
    stop(
      "`pedigree` must have at least three columns: `id`, `sire`, and `dam`.",
      call. = FALSE
    )
  }

  cols <- hs_pedigree_columns(pedigree)
  pedigree_ids <- hs_checked_ids(pedigree[[cols$id]], "`pedigree`")
  missing <- setdiff(unique(phenotype_ids), pedigree_ids)
  if (length(missing) > 0L) {
    stop(
      "`phenotypes` contain ID",
      if (length(missing) > 1L) "s" else "",
      " not present in `pedigree`: ",
      paste(missing, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  unique(pedigree_ids)
}

hs_optional_component_ids <- function(x, id, label) {
  if (is.null(x)) {
    return(character())
  }
  if (is.data.frame(x)) {
    if (id %in% names(x)) {
      return(unique(hs_checked_ids(x[[id]], label)))
    }
    if (hs_has_explicit_rownames(x)) {
      return(unique(hs_checked_ids(rownames(x), label)))
    }
    stop(
      label,
      " must contain ID column `",
      id,
      "` or non-missing row names.",
      call. = FALSE
    )
  }
  if (is.matrix(x)) {
    if (is.null(rownames(x))) {
      stop(
        label,
        " matrix must have individual IDs as row names.",
        call. = FALSE
      )
    }
    return(unique(hs_checked_ids(rownames(x), label)))
  }

  stop(label, " must be a data frame or matrix when supplied.", call. = FALSE)
}

hs_has_explicit_rownames <- function(x) {
  row_names <- attr(x, "row.names")
  !is.null(row_names) && !is.integer(row_names)
}

hs_validate_optional_data_frame <- function(x, label) {
  if (!is.null(x) && !is.data.frame(x)) {
    stop(label, " must be a data frame when supplied.", call. = FALSE)
  }
  invisible(TRUE)
}

hs_validate_annotation <- function(annotation, expression, id, annotation_id) {
  hs_validate_optional_data_frame(annotation, "`annotation`")

  if (is.null(annotation)) {
    if (!is.null(annotation_id)) {
      stop(
        "`annotation_id` can be supplied only when `annotation` is supplied.",
        call. = FALSE
      )
    }
    return(NULL)
  }
  if (is.null(annotation_id)) {
    return(NULL)
  }
  if (
    !is.character(annotation_id) ||
      length(annotation_id) != 1L ||
      is.na(annotation_id) ||
      annotation_id == ""
  ) {
    stop(
      "`annotation_id` must be one non-empty column name when supplied.",
      call. = FALSE
    )
  }
  if (!annotation_id %in% names(annotation)) {
    stop(
      "`annotation_id` column `",
      annotation_id,
      "` was not found in `annotation`.",
      call. = FALSE
    )
  }

  annotation_features <- hs_checked_component_keys(
    annotation[[annotation_id]],
    "`annotation`",
    annotation_id
  )
  expression_features <- hs_expression_feature_ids(
    expression,
    id,
    require_names = TRUE
  )
  duplicate_annotation_features <- unique(annotation_features[
    duplicated(annotation_features)
  ])

  structure(
    list(
      key = annotation_id,
      annotation_features = unique(annotation_features),
      expression_features = unique(expression_features),
      expression_without_annotation = setdiff(
        unique(expression_features),
        unique(annotation_features)
      ),
      annotation_without_expression = setdiff(
        unique(annotation_features),
        unique(expression_features)
      ),
      duplicate_annotation_features = duplicate_annotation_features
    ),
    class = "hs_annotation_spec"
  )
}

hs_expression_feature_ids <- function(
  expression,
  id,
  require_names = FALSE
) {
  if (is.null(expression)) {
    return(character())
  }
  if (is.matrix(expression)) {
    feature_ids <- colnames(expression)
    if (is.null(feature_ids)) {
      if (isTRUE(require_names)) {
        stop(
          "`expression` matrix must have feature IDs as column names when ",
          "`annotation_id` is supplied.",
          call. = FALSE
        )
      }
      return(rep("", ncol(expression)))
    }
    return(feature_ids)
  }
  if (is.data.frame(expression)) {
    feature_ids <- names(expression)
    if (id %in% feature_ids) {
      feature_ids <- feature_ids[feature_ids != id]
    }
    if (length(feature_ids) == 0L && isTRUE(require_names)) {
      stop(
        "`expression` must contain at least one feature column when ",
        "`annotation_id` is supplied.",
        call. = FALSE
      )
    }
    return(feature_ids)
  }
  character()
}

hs_expression_row_count <- function(expression) {
  if (is.matrix(expression) || is.data.frame(expression)) {
    return(nrow(expression))
  }
  NA_integer_
}

hs_expression_component_type <- function(expression) {
  if (is.matrix(expression)) {
    return("matrix")
  }
  if (is.data.frame(expression)) {
    return("data.frame")
  }
  "unknown"
}

hs_validate_environment <- function(environment, phenotypes, environment_id) {
  hs_validate_optional_data_frame(environment, "`environment`")

  if (is.null(environment)) {
    if (!is.null(environment_id)) {
      stop(
        "`environment_id` can be supplied only when `environment` is supplied.",
        call. = FALSE
      )
    }
    return(NULL)
  }
  if (is.null(environment_id)) {
    return(NULL)
  }
  if (
    !is.character(environment_id) ||
      length(environment_id) != 1L ||
      is.na(environment_id) ||
      environment_id == ""
  ) {
    stop(
      "`environment_id` must be one non-empty column name when supplied.",
      call. = FALSE
    )
  }
  if (!environment_id %in% names(phenotypes)) {
    stop(
      "`environment_id` column `",
      environment_id,
      "` was not found in `phenotypes`.",
      call. = FALSE
    )
  }
  if (!environment_id %in% names(environment)) {
    stop(
      "`environment_id` column `",
      environment_id,
      "` was not found in `environment`.",
      call. = FALSE
    )
  }

  phenotype_environment_ids <- hs_checked_component_keys(
    phenotypes[[environment_id]],
    "`phenotypes`",
    environment_id
  )
  environment_ids <- hs_checked_component_keys(
    environment[[environment_id]],
    "`environment`",
    environment_id
  )
  duplicate_environment_ids <- unique(environment_ids[
    duplicated(environment_ids)
  ])

  structure(
    list(
      key = environment_id,
      phenotype_environment_ids = unique(phenotype_environment_ids),
      environment_ids = unique(environment_ids),
      phenotypes_without_environment = setdiff(
        unique(phenotype_environment_ids),
        unique(environment_ids)
      ),
      environment_without_phenotypes = setdiff(
        unique(environment_ids),
        unique(phenotype_environment_ids)
      ),
      duplicate_environment_ids = duplicate_environment_ids
    ),
    class = "hs_environment_spec"
  )
}

hs_checked_component_keys <- function(x, label, key) {
  values <- as.character(x)
  if (any(is.na(values) | values == "")) {
    stop(
      label,
      " column `",
      key,
      "` cannot contain missing or empty values.",
      call. = FALSE
    )
  }
  values
}

hs_validate_marker_map <- function(markers) {
  if (is.null(markers)) {
    return(NULL)
  }
  hs_validate_optional_data_frame(markers, "`markers`")

  nm <- names(markers)
  lower <- tolower(nm)
  pick <- function(candidates) {
    hit <- which(lower %in% candidates)
    if (length(hit) == 0L) {
      return(NA_integer_)
    }
    hit[[1L]]
  }

  cols <- list(
    marker = pick(c("marker", "marker_id", "snp", "snp_id", "id")),
    chromosome = pick(c("chromosome", "chr", "chrom")),
    position = pick(c("position", "pos", "bp", "base_pair"))
  )

  missing <- names(cols)[is.na(unlist(cols, use.names = FALSE))]
  if (length(missing) > 0L) {
    stop(
      "`markers` must contain marker, chromosome, and position columns. ",
      "Recognized aliases include `marker`, `snp`, or `id`; `chromosome`, ",
      "`chr`, or `chrom`; and `position`, `pos`, `bp`, or `base_pair`. ",
      "Missing: ",
      paste(missing, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  marker_ids <- hs_checked_ids(markers[[cols$marker]], "`markers`")
  if (anyDuplicated(marker_ids)) {
    stop("`markers` contains duplicate marker IDs.", call. = FALSE)
  }

  chromosome <- as.character(markers[[cols$chromosome]])
  if (any(is.na(chromosome) | chromosome == "")) {
    stop(
      "`markers` chromosome column cannot contain missing or empty values.",
      call. = FALSE
    )
  }

  position <- suppressWarnings(as.numeric(as.character(
    markers[[cols$position]]
  )))
  if (
    any(is.na(position)) ||
      any(!is.finite(position)) ||
      any(position < 0)
  ) {
    stop(
      "`markers` position column must contain finite non-negative numeric ",
      "positions.",
      call. = FALSE
    )
  }

  structure(
    list(
      columns = cols,
      marker_ids = marker_ids,
      chromosome = chromosome,
      position = position
    ),
    class = "hs_marker_map_spec"
  )
}

hs_validate_genotype_marker_alignment <- function(genotypes, id, marker_spec) {
  if (is.null(genotypes) || is.null(marker_spec)) {
    return(NULL)
  }

  genotype_markers <- hs_genotype_marker_ids(genotypes, id)
  duplicated <- unique(genotype_markers[duplicated(genotype_markers)])
  if (length(duplicated) > 0L) {
    stop(
      "`genotypes` contains duplicate marker column",
      if (length(duplicated) > 1L) "s" else "",
      ": ",
      paste(duplicated, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  missing_from_map <- setdiff(genotype_markers, marker_spec$marker_ids)
  missing_from_genotypes <- setdiff(marker_spec$marker_ids, genotype_markers)
  if (length(missing_from_map) > 0L || length(missing_from_genotypes) > 0L) {
    details <- c(
      if (length(missing_from_map) > 0L) {
        paste0(
          "missing from `markers`: ",
          paste(missing_from_map, collapse = ", ")
        )
      },
      if (length(missing_from_genotypes) > 0L) {
        paste0(
          "missing from `genotypes`: ",
          paste(missing_from_genotypes, collapse = ", ")
        )
      }
    )
    stop(
      "`genotypes` marker columns must match `markers` marker IDs exactly; ",
      paste(details, collapse = "; "),
      ".",
      call. = FALSE
    )
  }

  structure(
    list(
      marker_ids = genotype_markers,
      marker_map_index = match(genotype_markers, marker_spec$marker_ids)
    ),
    class = "hs_genotype_marker_spec"
  )
}

hs_genotype_marker_ids <- function(genotypes, id) {
  if (is.matrix(genotypes)) {
    marker_ids <- colnames(genotypes)
    if (is.null(marker_ids)) {
      stop(
        "`genotypes` matrix must have marker IDs as column names when ",
        "`markers` is supplied.",
        call. = FALSE
      )
    }
    return(hs_checked_ids(marker_ids, "`genotypes` marker columns"))
  }

  if (is.data.frame(genotypes)) {
    marker_columns <- names(genotypes)
    if (id %in% marker_columns) {
      marker_columns <- marker_columns[marker_columns != id]
    }
    if (length(marker_columns) == 0L) {
      stop(
        "`genotypes` must contain at least one marker column when `markers` ",
        "is supplied.",
        call. = FALSE
      )
    }
    return(hs_checked_ids(marker_columns, "`genotypes` marker columns"))
  }

  character()
}

hs_genotype_row_count <- function(genotypes) {
  if (is.matrix(genotypes) || is.data.frame(genotypes)) {
    return(nrow(genotypes))
  }
  NA_integer_
}

hs_genotype_missing_value_count <- function(genotypes, id) {
  if (is.matrix(genotypes)) {
    return(sum(is.na(genotypes)))
  }
  if (is.data.frame(genotypes)) {
    marker_columns <- seq_along(genotypes)
    if (id %in% names(genotypes)) {
      marker_columns <- marker_columns[names(genotypes) != id]
    }
    payload <- genotypes[marker_columns]
    return(sum(vapply(payload, function(x) sum(is.na(x)), integer(1L))))
  }
  NA_integer_
}

hs_genotype_component_type <- function(genotypes) {
  if (is.matrix(genotypes)) {
    return("matrix")
  }
  if (is.data.frame(genotypes)) {
    return("data.frame")
  }
  "unknown"
}

hs_print_component_count <- function(label, ids) {
  if (length(ids) > 0L) {
    cat("  ", label, ": ", length(ids), "\n", sep = "")
  }
}
