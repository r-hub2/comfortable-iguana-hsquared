hs_build_model_spec <- function(
  formula,
  data,
  family,
  REML,
  allow_families = "gaussian"
) {
  env <- environment(formula)
  if (is.null(env)) {
    env <- parent.frame()
  }
  model_data <- hs_model_data_context(data, env)
  data <- model_data$data
  env <- model_data$env

  hs_validate_model_inputs(formula, data, family, REML, allow_families)

  rhs_terms <- hs_split_additive_rhs(formula[[3L]])

  # A bare `.` (all-other-columns shorthand) would otherwise reach
  # `model.frame()` and abort with the cryptic base-R error "'.' in formula and
  # no 'data' argument". Reject it by name and point to explicit fixed terms.
  hs_check_dot_term(rhs_terms)

  # Planned QG markers the parser now consumes: `permanent()`/`common_env()`/
  # `maternal_genetic()` as the second random effect of an opt-in two-effect
  # model, and `genomic()` as an opt-in primary genomic effect. Every other
  # planned marker still errors as not implemented.
  planned_pos <- which(vapply(
    rhs_terms,
    function(e) {
      hs_is_planned_marker_call(e) &&
        !hs_is_second_effect_call(e) &&
        !hs_is_relinv_primary_call(e) &&
        !hs_is_metafounder_primary_call(e)
    },
    logical(1L)
  ))
  if (length(planned_pos) > 0L) {
    hs_stop_planned_marker(rhs_terms[[planned_pos[[1L]]]])
  }

  # The primary effect is exactly one of `animal()` (pedigree), a
  # supplied-relationship-inverse term (`genomic()` with `Ginv`, or
  # `single_step()` with `Hinv`), or the opt-in supplied-Gamma metafounder
  # relationship term.
  animal_pos <- which(vapply(rhs_terms, hs_is_animal_call, logical(1L)))
  relinv_pos <- which(vapply(
    rhs_terms,
    hs_is_relinv_primary_call,
    logical(1L)
  ))
  metafounder_pos <- which(vapply(
    rhs_terms,
    hs_is_metafounder_primary_call,
    logical(1L)
  ))
  primary_pos <- c(animal_pos, relinv_pos, metafounder_pos)

  if (length(primary_pos) == 0L) {
    hs_abort_unsupported_syntax(
      "`formula` must contain exactly one primary term: ",
      "`animal(1 | id, pedigree = ped)`, `genomic(1 | id, Ginv = Ginv)`, ",
      "`metafounder(1 | id, pedigree = ped, group = mf_group, Gamma = Gamma)`, ",
      "`relmat(1 | id, K = K)`, or `precision(1 | id, Q = Q)`.",
      call. = FALSE
    )
  }
  if (length(primary_pos) > 1L) {
    hs_abort_unsupported_syntax(
      "`formula` can contain only one primary effect ",
      "(`animal()`, `genomic()`, `single_step()`, `metafounder()`, ",
      "`relmat()`, or `precision()`).",
      call. = FALSE
    )
  }

  if (length(animal_pos) == 1L) {
    primary_type <- "animal"
    primary_spec <- hs_parse_animal_call(
      rhs_terms[[animal_pos]],
      data,
      env,
      model_data = model_data
    )
  } else if (length(relinv_pos) == 1L) {
    primary_spec <- hs_parse_relinv_primary_call(
      rhs_terms[[relinv_pos]],
      data,
      env,
      model_data = model_data
    )
    primary_type <- primary_spec$type
  } else {
    primary_type <- "metafounder"
    primary_spec <- hs_parse_metafounder_call(
      rhs_terms[[metafounder_pos]],
      data,
      env,
      model_data = model_data
    )
  }

  second_pos <- which(vapply(rhs_terms, hs_is_second_effect_call, logical(1L)))
  second_spec <- NULL
  if (length(second_pos) > 0L && !identical(primary_type, "animal")) {
    hs_abort_unsupported_syntax(
      "A second random effect (`permanent()`/`common_env()`/",
      "`maternal_genetic()`) requires an `animal()` primary term, not `",
      primary_type,
      "()`.",
      call. = FALSE
    )
  }
  if (length(second_pos) > 1L) {
    hs_abort_unsupported_syntax(
      "`formula` can contain at most one additional random effect ",
      "(`permanent()`, `common_env()`, or `maternal_genetic()`) alongside ",
      "`animal()`.",
      call. = FALSE
    )
  }
  if (length(second_pos) == 1L) {
    second_spec <- hs_parse_second_effect_call(
      rhs_terms[[second_pos]],
      data,
      primary_spec
    )
  }

  # Any leftover bar term (a bare `(... | group)` random effect that is not a
  # recognized named effect). A bare random INTERCEPT `(1 | g)` is accepted as an
  # opt-in i.i.d. random effect (the SAME model as `common_env(1 | g)`); it is
  # routed to the engine's arbitrary-N independent-effect estimator. A random
  # SLOPE `(x | g)` or a correlated/uncorrelated slope term `(x || g)` is still
  # REJECTED: the engine has no slope/correlated estimator, so accepting one
  # would over-claim a capability that is not implemented.
  leftover_pos <- setdiff(seq_along(rhs_terms), c(primary_pos, second_pos))
  # A `(... || group)` correlated/uncorrelated slope term (top-level `||`) is not
  # a `|` bar term, so it must be rejected explicitly before it reaches
  # model.frame() with a cryptic base-R error.
  double_bar_pos <- leftover_pos[vapply(
    rhs_terms[leftover_pos],
    hs_is_double_bar_expr,
    logical(1L)
  )]
  if (length(double_bar_pos) > 0L) {
    hs_stop_unsupported_random_effect(rhs_terms[[double_bar_pos[[1L]]]])
  }
  bar_pos <- leftover_pos[vapply(
    rhs_terms[leftover_pos],
    hs_is_bar_expr,
    logical(1L)
  )]
  iid_effects <- list()
  for (pos in bar_pos) {
    if (hs_is_bare_iid_intercept(rhs_terms[[pos]])) {
      iid_effects[[length(iid_effects) + 1L]] <- hs_parse_bare_iid_call(
        rhs_terms[[pos]],
        data
      )
    } else {
      hs_stop_unsupported_random_effect(rhs_terms[[pos]])
    }
  }

  # A recognized effect/marker call nested inside an interaction or function
  # term (e.g. `x:dominance(1 | id)`) is not a top-level random effect, so it
  # would otherwise leak into the fixed-effect model.frame and abort with a
  # cryptic base-R error. Reject it here with the named term.
  for (pos in leftover_pos) {
    hs_check_nested_effect(rhs_terms[[pos]])
  }

  # Accepted bare-i.i.d. `(1 | g)` terms are random effects, not fixed terms;
  # drop them (alongside the primary and second effects) so they never leak into
  # the fixed-effect model.matrix().
  fixed_terms <- rhs_terms[-c(primary_pos, second_pos, bar_pos)]

  # An `offset()` term is silently dropped from the fixed design (and the bridge
  # payload) by `model.matrix()`. Reject it by name, mirroring the bare-bar and
  # nested-effect guards, rather than quietly ignoring the user's offset.
  hs_check_offset_term(fixed_terms)

  fixed_formula <- formula
  fixed_formula[[3L]] <- hs_rebuild_additive_rhs(fixed_terms)
  environment(fixed_formula) <- env

  model_frame <- stats::model.frame(
    fixed_formula,
    data = data,
    na.action = stats::na.pass,
    drop.unused.levels = FALSE
  )

  response_pos <- attr(stats::terms(fixed_formula), "response")
  fixed_frame <- model_frame
  if (!is.null(response_pos) && response_pos > 0L) {
    fixed_frame <- model_frame[-response_pos]
  }
  if (length(fixed_frame) > 0L && anyNA(fixed_frame)) {
    stop(
      "Missing values in fixed-effect variables are not implemented. ",
      "For multivariate responses, use `NA` only inside the `cbind()` ",
      "response matrix to mark missing trait records.",
      call. = FALSE
    )
  }

  response <- hs_build_response_spec(
    formula[[2L]],
    stats::model.response(
      model_frame
    ),
    family
  )

  fixed_terms_obj <- stats::terms(fixed_formula)

  # `model.matrix()` aborts with a cryptic base-R error ("contrasts can be
  # applied only to factors with 2 or more levels") on a zero-row model.frame or
  # a factor/character fixed effect with fewer than two observed levels. Catch
  # both first and name the offending term, consistent with the rank-deficient
  # wording below.
  hs_validate_fixed_model_frame(fixed_frame)

  fixed_design <- stats::model.matrix(fixed_terms_obj, data = model_frame)
  hs_validate_fixed_design(fixed_design)

  random <- list()
  random[[primary_type]] <- primary_spec
  bridge_target <- if (primary_type %in% c("relmat", "precision")) {
    # Experimental supplied-relationship primary effect: the supplied inverse
    # (Kinv from relmat's K, or the precision Q) fits through the SAME covered
    # relationship-inverse REML path as genomic()'s Ginv.
    "fit_ai_reml(y, X, Z, Kinv; method = :REML)"
  } else if (primary_type %in% c("genomic", "single_step")) {
    if (identical(primary_spec$source, "metafounder_construct")) {
      paste0(
        "fit_metafounder_single_step_reml(y, X, Z, pedigree, ",
        "group_of, Gamma, G, genotyped_rows; method = :REML)"
      )
    } else if (identical(primary_spec$source, "construct")) {
      paste0(
        "fit_single_step_reml(y, X, Z, Ainv, A, G, genotyped_rows; ",
        "method = :REML)"
      )
    } else if (identical(primary_spec$source, "markers")) {
      paste0(
        "fit_ai_reml(y, X, Z, ",
        "genomic_relationship_inverse(genomic_relationship_matrix(markers)); ",
        "method = :REML)"
      )
    } else {
      sprintf(
        "fit_ai_reml(y, X, Z, %s; method = :REML)",
        if (identical(primary_type, "genomic")) "Ginv" else "Hinv"
      )
    }
  } else if (identical(primary_type, "metafounder")) {
    paste0(
      "metafounder_animal_model(y, X, Z, pedigree, group_of, Gamma, ",
      "sigma_a2, sigma_e2)"
    )
  } else if (identical(primary_spec$design, "random_regression")) {
    "fit_random_regression_reml(y, X, Phi, Z, Ainv; ids = ped.ids)"
  } else {
    "fit_animal_model(y, X, Z, Ainv; method = :REML)"
  }
  if (!is.null(second_spec)) {
    random[[second_spec$type]] <- second_spec
    bridge_target <- if (identical(second_spec$type, "permanent")) {
      "fit_repeatability_reml(y, X, Z, Ainv; method = :REML)"
    } else {
      "fit_two_effect_reml(y, X, Z, Ainv, Z2, Ainv2; method = :REML)"
    }
  }
  # Bare `(1 | g)` i.i.d. random effects (opt-in arbitrary-N independent-effect
  # model). Held as a LIST (not keyed by a string type) so that two `(1 | g)`
  # terms on different grouping columns are BOTH retained rather than colliding
  # on a shared `random[[type]]` slot. Scope: an `animal()` pedigree primary plus
  # independent i.i.d. blocks (the covered `fit_multi_effect_reml` estimand); a
  # genomic/single-step/metafounder primary or a random-regression design does
  # not (yet) combine with additional i.i.d. effects in the R bridge.
  if (length(iid_effects) > 0L) {
    if (!identical(primary_type, "animal")) {
      hs_abort_unsupported_syntax(
        "A bare `(1 | group)` i.i.d. random effect currently combines only with ",
        "an `animal(1 | id, pedigree = ped)` primary term, not `",
        primary_type,
        "()`. Genomic, single-step, and metafounder multi-effect models are ",
        "planned, not implemented.",
        call. = FALSE
      )
    }
    if (identical(primary_spec$design, "random_regression")) {
      hs_abort_unsupported_syntax(
        "A bare `(1 | group)` i.i.d. random effect does not (yet) combine with ",
        "an `animal(rr(...) | id)` random-regression design; both are opt-in, ",
        "and the combination is planned, not implemented.",
        call. = FALSE
      )
    }
    random$iid_effects <- iid_effects
    bridge_target <- paste0(
      "fit_multi_effect_reml(y, X, [(Z_i, Ainv_i) ...]; ",
      "ids = [block.ids ...])"
    )
  }
  # The opt-in random-regression design is a single-effect, univariate model. It
  # does not (yet) combine with a second random effect or a multivariate cbind()
  # response; both are planned.
  if (
    identical(primary_spec$design, "random_regression") && !is.null(second_spec)
  ) {
    hs_abort_unsupported_syntax(
      "A `rr(...)` random-regression term in `animal()` is a single-effect, ",
      "opt-in model; combining it with a second random effect ",
      "(`permanent()`/`common_env()`/`maternal_genetic()`) is planned, not ",
      "implemented.",
      call. = FALSE
    )
  }
  if (isTRUE(response$multivariate)) {
    if (
      !identical(primary_type, "animal") ||
        !is.null(second_spec) ||
        length(iid_effects) > 0L
    ) {
      # D3 / hsquared#237: name closest live paths rather than a bare planned
      # fence when the blocker is a permanent-environment second effect.
      if (!is.null(second_spec) && identical(second_spec$type, "permanent")) {
        hs_abort_unsupported_syntax(
          "`cbind(...)` with `permanent()` is not implemented (hsquared#237). ",
          "Closest live paths: (1) univariate repeatability -- drop `cbind` and ",
          "use `control = hs_control(engine = \"julia\", engine_control = list(",
          "target = \"repeatability\"))` with `animal(...) + permanent(1 | id)`; ",
          "(2) multivariate without PE -- keep `cbind(...) ~ fixed + animal(1 | ",
          "id, pedigree = ped)` and drop `permanent()`. Full MV+PE fit is ",
          "deferred; public_covered_count stays 7.",
          call. = FALSE
        )
      }
      hs_abort_unsupported_syntax(
        "The multivariate path currently supports only ",
        "`cbind(...) ~ fixed + animal(1 | id, pedigree = ped)`. ",
        "Multivariate genomic, single-step, second-effect, and multi-effect ",
        "models are planned, not implemented.",
        call. = FALSE
      )
    }
    if (identical(primary_spec$design, "random_regression")) {
      hs_abort_unsupported_syntax(
        "Multivariate random regression (`cbind(...)` response with a ",
        "`rr(...)` term) is planned, not implemented. Use a univariate ",
        "response with the opt-in `target = \"random_regression\"` path.",
        call. = FALSE
      )
    }
    bridge_target <- "fit_multivariate_reml(Y, X, Z, Ainv; method = :REML)"
  }

  list(
    formula = formula,
    fixed_formula = fixed_formula,
    family = list(family = family$family, link = family$link),
    method = if (isTRUE(REML)) "REML" else "ML",
    response = response,
    fixed = list(
      terms = attr(fixed_terms_obj, "term.labels"),
      design = fixed_design,
      contrasts = attr(fixed_design, "contrasts")
    ),
    random = random,
    bridge = list(
      status = "planned",
      engine = "HSquared.jl",
      target = bridge_target
    )
  )
}

hs_build_response_spec <- function(lhs, response, family = NULL) {
  is_cbind <- hs_is_call(hs_unwrap_parentheses(lhs), "cbind")

  # A `cbind(successes, failures)` matrix under `family = binomial()` is the
  # canonical binomial-COUNTS response (the R glm convention), NOT a two-trait
  # multivariate Gaussian. Detect it before the family-blind multivariate branch
  # so it is never silently fitted as two Gaussian traits.
  if (
    is.matrix(response) &&
      ncol(response) == 2L &&
      is_cbind &&
      !is.null(family) &&
      identical(family$family, "binomial") &&
      identical(family$link, "logit")
  ) {
    return(hs_build_binomial_counts_response(lhs, response))
  }

  multivariate <- is.matrix(response)
  if (multivariate && !is_cbind) {
    stop(
      "Multivariate responses must use `cbind(trait1, trait2, ...)` on the ",
      "left-hand side.",
      call. = FALSE
    )
  }

  if (multivariate) {
    hs_validate_cbind_bare_columns(hs_unwrap_parentheses(lhs))
    values <- unname(as.matrix(response))
    if (!is.numeric(values)) {
      stop("Multivariate `cbind()` responses must be numeric.", call. = FALSE)
    }
    if (ncol(values) < 2L) {
      stop(
        "Multivariate `cbind()` responses require at least two trait columns.",
        call. = FALSE
      )
    }
    observed <- !is.na(values)
    if (any(observed & !is.finite(values))) {
      stop(
        "Observed multivariate response values must be finite. Use `NA` for ",
        "missing trait records.",
        call. = FALSE
      )
    }
    if (any(colSums(observed) == 0L)) {
      stop(
        "Each trait in a multivariate `cbind()` response must have at least ",
        "one observed value.",
        call. = FALSE
      )
    }
    trait_names <- colnames(response)
    if (
      is.null(trait_names) || anyNA(trait_names) || any(!nzchar(trait_names))
    ) {
      trait_names <- all.vars(lhs)
    }
    if (length(trait_names) != ncol(values)) {
      trait_names <- paste0("trait", seq_len(ncol(values)))
    }
    hs_validate_multivariate_trait_names(trait_names)
    colnames(values) <- trait_names
    return(list(
      name = hs_deparse(lhs),
      values = values,
      trait_names = trait_names,
      multivariate = TRUE
    ))
  }

  if (!is.numeric(response)) {
    stop(
      "The v0.1 parser supports numeric Gaussian responses only.",
      call. = FALSE
    )
  }
  if (anyNA(response)) {
    stop(
      "Missing values in the response are not implemented for univariate ",
      "models. Use a multivariate `cbind()` response only when missing trait ",
      "records are intended.",
      call. = FALSE
    )
  }
  if (any(!is.finite(response))) {
    stop("Response values must be finite.", call. = FALSE)
  }

  list(
    name = hs_deparse(lhs),
    values = as.numeric(response),
    trait_names = NULL,
    multivariate = FALSE
  )
}

# `cbind(successes, failures) ~ ...` under `family = binomial()`: a binomial-counts
# response. The success counts are the modelled response; the per-record number of
# trials is successes + failures, kept as a per-record integer vector that may
# vary across records. A constant vector is the common-trial special case; an
# all-ones vector is the Bernoulli reduction.
hs_build_binomial_counts_response <- function(lhs, response) {
  hs_validate_cbind_bare_columns(hs_unwrap_parentheses(lhs))
  values <- unname(as.matrix(response))
  if (!is.numeric(values)) {
    stop(
      "A `cbind(successes, failures)` binomial response must be numeric.",
      call. = FALSE
    )
  }
  if (anyNA(values) || any(!is.finite(values))) {
    stop(
      "A `cbind(successes, failures)` binomial response must be finite (no NA).",
      call. = FALSE
    )
  }
  successes <- values[, 1L]
  failures <- values[, 2L]
  if (any(values < 0) || any(values != round(values))) {
    stop(
      "`cbind(successes, failures)` binomial counts must be non-negative ",
      "integers.",
      call. = FALSE
    )
  }
  totals <- successes + failures
  if (any(totals < 1)) {
    stop(
      "Each `cbind(successes, failures)` record must have at least one trial ",
      "(successes + failures >= 1).",
      call. = FALSE
    )
  }
  # Per-record trial counts (may vary across records). A constant vector is the
  # common-trial special case; an all-ones vector reduces to Bernoulli. The R
  # bridge parses and carries the full vector here; live varying-trial fitting is
  # verified separately through the engine round-trip.
  n_trials <- as.integer(totals)
  trait_names <- all.vars(hs_unwrap_parentheses(lhs))
  list(
    name = hs_deparse(lhs),
    values = as.numeric(successes),
    trait_names = if (length(trait_names) >= 1L) trait_names[[1L]] else NULL,
    multivariate = FALSE,
    binomial_counts = TRUE,
    n_trials = as.integer(n_trials)
  )
}

hs_validate_multivariate_trait_names <- function(trait_names) {
  missing_names <- is.na(trait_names) | !nzchar(trait_names)
  duplicate_names <- duplicated(trait_names) |
    duplicated(trait_names, fromLast = TRUE)

  if (!any(missing_names) && !any(duplicate_names)) {
    return(invisible(trait_names))
  }

  detail <- character()
  if (any(missing_names)) {
    detail <- c(detail, "empty or missing names")
  }
  if (any(duplicate_names)) {
    duplicates <- unique(trait_names[duplicate_names & !missing_names])
    detail <- c(
      detail,
      paste0("duplicate names: ", paste(duplicates, collapse = ", "))
    )
  }

  stop(
    "Multivariate `cbind()` responses require unique, non-empty trait names",
    if (length(detail) > 0L) {
      paste0(" (", paste(detail, collapse = "; "), ")")
    },
    ". Rename or wrap response columns before fitting.",
    call. = FALSE
  )
}

# Reject a multivariate `cbind()` response whose arguments are not bare trait
# column symbols (e.g. `cbind(y1, y1 + y2)` or `cbind(log(y1), y2)`). A derived
# column would otherwise be mislabelled from `all.vars(lhs)`, giving a
# confidently wrong trait name. Bare-column `cbind(y1, y2)` is left untouched.
hs_validate_cbind_bare_columns <- function(lhs) {
  args <- as.list(lhs)[-1L]
  bare <- vapply(args, is.symbol, logical(1L))
  if (all(bare)) {
    return(invisible(TRUE))
  }
  offending <- vapply(args[!bare], hs_deparse, character(1L))
  hs_abort_unsupported_syntax(
    "v0.1 multivariate responses require bare trait columns inside ",
    "`cbind()`. Derived or transformed column",
    if (length(offending) > 1L) "s " else " ",
    paste(sprintf("`%s`", offending), collapse = ", "),
    if (length(offending) > 1L) " are " else " is ",
    "not supported. Create the columns in `data` first, then list them as ",
    "`cbind(trait1, trait2, ...)`.",
    call. = FALSE
  )
}

# Reject a zero-row model.frame or a factor/character fixed effect with fewer
# than two observed levels before `model.matrix()` is built. `fixed_frame` is
# the model.frame with the response column removed (the fixed-effect columns).
# Without this guard base R aborts with "contrasts can be applied only to
# factors with 2 or more levels", which names neither the variable nor the
# cause.
hs_validate_fixed_model_frame <- function(fixed_frame) {
  if (length(fixed_frame) == 0L) {
    return(invisible(TRUE))
  }
  if (nrow(fixed_frame) == 0L) {
    stop(
      "The fixed-effect model frame has zero rows. Provide `data` with at ",
      "least one observed record before fitting.",
      call. = FALSE
    )
  }
  for (nm in names(fixed_frame)) {
    column <- fixed_frame[[nm]]
    if (!is.factor(column) && !is.character(column)) {
      next
    }
    n_levels <- length(unique(column[!is.na(column)]))
    if (n_levels < 2L) {
      stop(
        "The fixed-effect term `",
        nm,
        "` has ",
        if (n_levels == 0L) "no observed levels" else "only one observed level",
        ", so it cannot be coded as a contrast. Drop it or supply at least two ",
        "levels before fitting.",
        call. = FALSE
      )
    }
  }
  invisible(TRUE)
}

# Reject an `offset()` term in the fixed RHS. `model.matrix()` would otherwise
# silently drop it from both the design matrix and the bridge payload, so the
# offset the user wrote would have no effect. Mirror the bare-bar and
# nested-effect guards by naming the offending term.
hs_check_offset_term <- function(fixed_terms) {
  offset_pos <- which(vapply(
    fixed_terms,
    function(e) hs_is_call(hs_unwrap_parentheses(e), "offset"),
    logical(1L)
  ))
  if (length(offset_pos) == 0L) {
    return(invisible(NULL))
  }
  hs_abort_unsupported_syntax(
    "Offset terms (`",
    hs_deparse(fixed_terms[[offset_pos[[1L]]]]),
    "`) are planned, not implemented in v0.1. Remove the `offset()` term ",
    "before fitting. Run `formula_status()` for the live list of which terms ",
    "parse and which fit.",
    call. = FALSE
  )
}

# Reject a bare `.` (the all-other-columns shorthand) among the split RHS
# terms. Left alone it would reach `model.frame()` and abort with the cryptic
# base-R error "'.' in formula and no 'data' argument".
hs_check_dot_term <- function(rhs_terms) {
  is_dot <- vapply(
    rhs_terms,
    function(e) is.symbol(e) && identical(as.character(e), "."),
    logical(1L)
  )
  if (!any(is_dot)) {
    return(invisible(NULL))
  }
  hs_abort_unsupported_syntax(
    "The `.` all-columns shorthand is not supported. List the fixed-effect ",
    "terms explicitly, for example ",
    "`y ~ sex + age + animal(1 | id, pedigree = ped)`.",
    call. = FALSE
  )
}

hs_validate_fixed_design <- function(fixed_design) {
  if (ncol(fixed_design) == 0L) {
    return(invisible(TRUE))
  }
  if (qr(fixed_design)$rank < ncol(fixed_design)) {
    stop(
      "The fixed-effect design matrix is rank deficient. Remove redundant ",
      "fixed-effect columns before fitting.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

hs_model_data_context <- function(data, env) {
  if (!inherits(data, "hs_data")) {
    return(list(
      data = data,
      env = env,
      is_hs_data = FALSE,
      components = list()
    ))
  }
  if (!is.data.frame(data$phenotypes)) {
    stop(
      "`data` is an `hs_data` object, but `data$phenotypes` is not a ",
      "data frame.",
      call. = FALSE
    )
  }

  components <- list(
    phenotypes = data$phenotypes,
    pedigree = data$pedigree,
    genotypes = data$genotypes,
    markers = data$markers,
    expression = data$expression,
    annotation = data$annotation,
    environment = data$environment
  )
  components <- components[!vapply(components, is.null, logical(1L))]

  list(
    data = data$phenotypes,
    env = list2env(components, parent = env),
    is_hs_data = TRUE,
    components = components,
    id = data$id
  )
}

hs_validate_model_inputs <- function(
  formula,
  data,
  family,
  REML,
  allow_families = "gaussian"
) {
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    stop("`formula` must be a two-sided formula.", call. = FALSE)
  }
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }
  if (!inherits(family, "family")) {
    stop(
      "`family` must be an R family object such as `gaussian()`.",
      call. = FALSE
    )
  }
  # gaussian fits everywhere (identity link); poisson(log)/binomial(logit) are
  # accepted only on the opt-in non-Gaussian path (`allow_families` widened by
  # `hsquared()` when `target = "nongaussian"`).
  expected_link <- c(gaussian = "identity", poisson = "log", binomial = "logit")
  link_for_family <- expected_link[family$family]
  family_ok <- family$family %in%
    allow_families &&
    !is.na(link_for_family) &&
    identical(family$link, unname(link_for_family))
  if (!family_ok) {
    hs_abort_unsupported_syntax(
      "The requested family `",
      hs_family_label(family),
      "` is not fitted on this path. The default `hsquared()` path fits ",
      "`family = gaussian()` (identity link). Non-Gaussian `poisson(log)` and ",
      "`binomial(logit)` (binary 0/1, or `cbind(successes, failures)` counts) ",
      "fit through the experimental, opt-in ",
      "`hs_control(engine = \"julia\", engine_control = list(target = ",
      "\"nongaussian\"))` path: a latent-scale GLMM (engine row V6-LAPLACE/VA, ",
      "partial) with a Laplace marginal likelihood or variational ELBO; ",
      "the 0.9 conditional three-field h2 contract remains experimental and ",
      "not coverage-calibrated. Use `model_spec()` ",
      "with `family = gaussian()` to inspect the contract without fitting.",
      call. = FALSE
    )
  }
  if (!is.logical(REML) || length(REML) != 1L || is.na(REML)) {
    stop("`REML` must be `TRUE` or `FALSE`.", call. = FALSE)
  }
  # Public hsquared() validate and default fit share one REML = TRUE rule.
  # Internal builders and model_spec() are not that door.
  hs_abort_reml_false_on_public_path(REML)

  invisible(TRUE)
}

hs_family_label <- function(family) {
  family_name <- family$family
  link_name <- family$link
  if (is.null(family_name) || !nzchar(family_name)) {
    family_name <- "<unknown>"
  }
  if (is.null(link_name) || !nzchar(link_name)) {
    link_name <- "<unknown>"
  }
  paste0(family_name, "(", link_name, ")")
}

hs_parse_metafounder_call <- function(call, data, env, model_data) {
  call <- hs_unwrap_parentheses(call)
  args <- as.list(call)[-1L]
  arg_names <- names(args)
  if (is.null(arg_names)) {
    arg_names <- rep("", length(args))
  }

  bar_candidates <- which(arg_names == "" | arg_names == "formula")
  if (length(bar_candidates) != 1L) {
    stop(
      "`metafounder()` must have one random-effect expression, for example ",
      "`metafounder(1 | id, pedigree = ped, group = mf_group, Gamma = Gamma)`.",
      call. = FALSE
    )
  }

  bar <- hs_unwrap_parentheses(args[[bar_candidates]])
  if (!hs_is_call(bar, "|") || length(bar) != 3L) {
    stop(
      "The first `metafounder()` argument must be a random-effect expression ",
      "such as `1 | id`.",
      call. = FALSE
    )
  }

  lhs <- hs_unwrap_parentheses(bar[[2L]])
  group_expr <- hs_unwrap_parentheses(bar[[3L]])
  if (!hs_is_one(lhs)) {
    hs_abort_unsupported_syntax(
      "Only random-intercept syntax `metafounder(1 | id, ...)` is implemented. ",
      "Metafounder slopes are planned, not implemented.",
      call. = FALSE
    )
  }
  if (!is.symbol(group_expr)) {
    stop(
      "The grouping variable in `metafounder()` must be a bare column name.",
      call. = FALSE
    )
  }

  group <- as.character(group_expr)
  if (!group %in% names(data)) {
    stop(
      "`metafounder()` grouping variable `",
      group,
      "` was not found in `data`.",
      call. = FALSE
    )
  }

  named_args <- args[arg_names != ""]
  accepted <- c("pedigree", "group", "Gamma")
  unsupported <- setdiff(names(named_args), accepted)
  if (length(unsupported) > 0L) {
    stop(
      "`metafounder()` argument",
      if (length(unsupported) > 1L) "s " else " ",
      paste(sprintf("`%s`", unsupported), collapse = ", "),
      if (length(unsupported) > 1L) {
        " are planned, not implemented."
      } else {
        " is planned, not implemented."
      },
      call. = FALSE
    )
  }
  if (!all(c("group", "Gamma") %in% names(named_args))) {
    stop(
      "`metafounder()` requires both `group` and `Gamma`: ",
      "`metafounder(1 | id, pedigree = ped, group = mf_group, Gamma = Gamma)`.",
      call. = FALSE
    )
  }

  pedigree_input <- hs_resolve_animal_pedigree(
    named_args,
    data,
    env,
    model_data
  )
  pedigree_spec <- hs_validate_pedigree(
    pedigree_input$data,
    data_ids = data[[group]],
    group = group
  )
  group_of <- hs_eval_metafounder_group(
    named_args$group,
    data,
    env,
    term = "metafounder"
  )
  group_of <- hs_validate_metafounder_group(
    group_of,
    pedigree_spec,
    term = "metafounder"
  )
  gamma <- hs_eval_metafounder_gamma(
    named_args$Gamma,
    data,
    env,
    term = "metafounder"
  )
  gamma <- hs_validate_metafounder_gamma(
    gamma,
    group_of,
    pedigree_spec,
    term = "metafounder"
  )

  list(
    type = "metafounder",
    term = hs_deparse(call),
    design = "intercept",
    group = group,
    values = as.character(data[[group]]),
    relationship = "metafounder",
    covariance = "scalar",
    pedigree_source = pedigree_input$source,
    pedigree = pedigree_spec,
    group_of = group_of,
    Gamma = gamma$Gamma,
    gamma_labels = gamma$labels,
    gamma_source = "supplied"
  )
}

hs_parse_animal_call <- function(call, data, env, model_data) {
  call <- hs_unwrap_parentheses(call)
  args <- as.list(call)[-1L]
  arg_names <- names(args)
  if (is.null(arg_names)) {
    arg_names <- rep("", length(args))
  }

  bar_candidates <- which(arg_names == "" | arg_names == "formula")
  if (length(bar_candidates) != 1L) {
    stop(
      "`animal()` must have one random-effect expression, for example ",
      "`animal(1 | id, pedigree = ped)`.",
      call. = FALSE
    )
  }

  bar <- hs_unwrap_parentheses(args[[bar_candidates]])
  if (!hs_is_call(bar, "|") || length(bar) != 3L) {
    stop(
      "The first `animal()` argument must be a random-effect expression ",
      "such as `1 | id`.",
      call. = FALSE
    )
  }

  lhs <- hs_unwrap_parentheses(bar[[2L]])
  group_expr <- hs_unwrap_parentheses(bar[[3L]])
  # `rr(covariate, order = k)` on the left-hand side is the opt-in, experimental
  # random-regression (reaction-norm) design: a per-record Legendre polynomial of
  # a within-individual covariate. Any other non-intercept left-hand side is the
  # planned random-slope / long-format syntax, still rejected.
  rr_spec <- NULL
  if (!hs_is_one(lhs)) {
    if (hs_is_call(lhs, "rr")) {
      rr_spec <- hs_parse_rr_lhs(lhs, data)
    } else {
      hs_stop_animal_non_intercept()
    }
  }
  if (!is.symbol(group_expr)) {
    stop(
      "The grouping variable in `animal()` must be a bare column name.",
      call. = FALSE
    )
  }

  group <- as.character(group_expr)
  if (!group %in% names(data)) {
    stop(
      "`animal()` grouping variable `",
      group,
      "` was not found in `data`.",
      call. = FALSE
    )
  }

  named_args <- args[arg_names != ""]
  unsupported <- setdiff(names(named_args), "pedigree")
  if (length(unsupported) > 0L) {
    if ("cov" %in% unsupported) {
      hs_stop_animal_covariance_arg()
    }
    stop(
      "`animal()` argument",
      if (length(unsupported) > 1L) "s " else " ",
      paste(sprintf("`%s`", unsupported), collapse = ", "),
      if (length(unsupported) > 1L) {
        " are planned, not implemented in v0.1."
      } else {
        " is planned, not implemented in v0.1."
      },
      call. = FALSE
    )
  }
  pedigree_input <- hs_resolve_animal_pedigree(
    named_args,
    data,
    env,
    model_data
  )
  pedigree_spec <- hs_validate_pedigree(
    pedigree_input$data,
    data_ids = data[[group]],
    group = group
  )

  spec <- list(
    type = "animal",
    term = hs_deparse(call),
    design = "intercept",
    group = group,
    values = as.character(data[[group]]),
    relationship = "pedigree",
    covariance = "scalar",
    pedigree_source = pedigree_input$source,
    pedigree = pedigree_spec
  )
  if (!is.null(rr_spec)) {
    spec$design <- "random_regression"
    spec$covariance <- "random_regression"
    spec$random_regression <- rr_spec
  }
  spec
}

# Parse `rr(covariate, order = k)` from the left-hand side of an `animal()`
# random-effect expression. `covariate` is a bare data column that varies within
# individual (repeated records); `order` is the number of Legendre coefficients
# (default 2 = intercept + slope). The covariate is standardized to [-1, 1] via
# its observed data range in the engine; the (lower, upper) bounds are recorded
# here so extractors can re-standardize a user-supplied `at =` on the original
# scale. The grammar is ratified by the twin on HSquared.jl#61.
hs_parse_rr_lhs <- function(lhs, data) {
  lhs <- hs_unwrap_parentheses(lhs)
  args <- as.list(lhs)[-1L]
  arg_names <- names(args)
  if (is.null(arg_names)) {
    arg_names <- rep("", length(args))
  }

  covariate_candidates <- which(arg_names == "" | arg_names == "covariate")
  if (length(covariate_candidates) != 1L) {
    stop(
      "`rr()` must name exactly one covariate column, for example ",
      "`rr(age, order = 2)`.",
      call. = FALSE
    )
  }
  covariate_expr <- hs_unwrap_parentheses(args[[covariate_candidates]])
  if (!is.symbol(covariate_expr)) {
    stop(
      "The `rr()` covariate must be a bare column name, for example ",
      "`rr(age, order = 2)`.",
      call. = FALSE
    )
  }
  covariate <- as.character(covariate_expr)
  if (!covariate %in% names(data)) {
    stop(
      "`rr()` covariate `",
      covariate,
      "` was not found in `data`.",
      call. = FALSE
    )
  }

  named_args <- args[arg_names != "" & arg_names != "covariate"]
  unsupported <- setdiff(names(named_args), "order")
  if (length(unsupported) > 0L) {
    stop(
      "`rr()` argument",
      if (length(unsupported) > 1L) "s " else " ",
      paste(sprintf("`%s`", unsupported), collapse = ", "),
      if (length(unsupported) > 1L) {
        " are planned, not implemented."
      } else {
        " is planned, not implemented."
      },
      " The opt-in random-regression grammar currently accepts only ",
      "`rr(covariate, order = k)`.",
      call. = FALSE
    )
  }

  order <- if ("order" %in% names(named_args)) named_args$order else 2
  order <- suppressWarnings(as.integer(order))
  if (length(order) != 1L || is.na(order) || order < 1L) {
    stop(
      "`rr(order = ...)` must be a single positive integer (the number of ",
      "Legendre coefficients; 2 = intercept + slope).",
      call. = FALSE
    )
  }

  values <- data[[covariate]]
  if (!is.numeric(values)) {
    stop(
      "`rr()` covariate `",
      covariate,
      "` must be numeric.",
      call. = FALSE
    )
  }
  if (anyNA(values) || any(!is.finite(values))) {
    stop(
      "`rr()` covariate `",
      covariate,
      "` must contain only finite, non-missing values.",
      call. = FALSE
    )
  }
  lower <- min(values)
  upper <- max(values)
  if (!(upper > lower)) {
    stop(
      "`rr()` covariate `",
      covariate,
      "` must vary (its observed range has zero width), so it cannot be ",
      "standardized to [-1, 1].",
      call. = FALSE
    )
  }

  list(
    covariate = covariate,
    values = as.numeric(values),
    order = order,
    lower = as.numeric(lower),
    upper = as.numeric(upper)
  )
}

hs_stop_animal_non_intercept <- function() {
  hs_abort_unsupported_syntax(
    "Only random-intercept syntax `animal(1 | id, pedigree = ped)` and the ",
    "opt-in random-regression syntax `animal(rr(covariate, order = k) | id, ",
    "pedigree = ped)` are implemented inside `animal()`. For the current opt-in ",
    "multivariate animal model, put traits on the left-hand side as ",
    "`cbind(trait1, trait2) ~ ... + animal(1 | id, pedigree = ped)` and use ",
    "`engine_control = list(target = \"multivariate\")`. Random regression uses ",
    "`engine_control = list(target = \"random_regression\")`. Long-format ",
    "`animal(trait | id, cov = ...)` and random-slope syntax are planned, not ",
    "implemented.",
    call. = FALSE
  )
}

hs_stop_animal_covariance_arg <- function() {
  hs_abort_unsupported_syntax(
    "`animal()` argument `cov` is reserved, not parsed. For the current ",
    "multivariate animal model, put traits on the left-hand side as ",
    "`cbind(trait1, trait2) ~ ... + animal(1 | id, pedigree = ped)`; that ",
    "routes to the multivariate fitter on the default path, with no ",
    "`engine`/`target` argument. Long-format `cov = us()`, `cov = diag()`, ",
    "and `cov = lowrank(K = 2)` remain planned grammar. `cov = fa(K = 2)` ",
    "is planned R grammar for factor-analytic G. Julia V4-FA is ",
    "engine-covered (HSquared.jl 60895208); that is not an R-public fit. ",
    "The R formula does not parse or fit it.",
    call. = FALSE
  )
}

hs_resolve_animal_pedigree <- function(named_args, data, env, model_data) {
  if ("pedigree" %in% names(named_args)) {
    return(list(
      data = hs_eval_pedigree(named_args$pedigree, data, env),
      source = "formula"
    ))
  }

  if (
    isTRUE(model_data$is_hs_data) &&
      !is.null(model_data$components$pedigree)
  ) {
    return(list(
      data = model_data$components$pedigree,
      source = "hs_data"
    ))
  }

  stop(
    "`animal()` requires `pedigree = ped` in the v0.1 parser, unless ",
    "`data` is an `hs_data()` object with a pedigree component.",
    call. = FALSE
  )
}

hs_validate_pedigree <- function(pedigree, data_ids, group) {
  if (!is.data.frame(pedigree) || ncol(pedigree) < 3L) {
    stop(
      "`pedigree` must be a data frame with at least three columns: ",
      "`id`, `sire`, and `dam`.",
      call. = FALSE
    )
  }

  cols <- hs_pedigree_columns(pedigree)
  ids <- as.character(pedigree[[cols$id]])
  sire <- hs_normalize_parent(pedigree[[cols$sire]])
  dam <- hs_normalize_parent(pedigree[[cols$dam]])

  if (any(is.na(ids) | ids == "" | ids == "0")) {
    stop(
      "`pedigree` individual IDs cannot be missing, empty, or `0`.",
      call. = FALSE
    )
  }
  if (anyDuplicated(ids)) {
    stop("`pedigree` contains duplicate individual IDs.", call. = FALSE)
  }

  parent_values <- unique(c(stats::na.omit(sire), stats::na.omit(dam)))
  missing_parents <- setdiff(parent_values, ids)
  if (length(missing_parents) > 0L) {
    stop(
      "`pedigree` contains known parent ID",
      if (length(missing_parents) > 1L) "s" else "",
      " not present as individual IDs: ",
      paste(missing_parents, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  self_parent <- (!is.na(sire) & sire == ids) | (!is.na(dam) & dam == ids)
  if (any(self_parent)) {
    stop(
      "`pedigree` cannot list an individual as its own parent.",
      call. = FALSE
    )
  }

  same_known_parent <- !is.na(sire) & !is.na(dam) & sire == dam
  if (any(same_known_parent)) {
    stop(
      "`pedigree` rows with the same known sire and dam are not supported ",
      "in v0.1. Selfing and non-standard inheritance are planned later.",
      call. = FALSE
    )
  }

  observed_ids <- as.character(data_ids)
  if (any(is.na(observed_ids) | observed_ids == "" | observed_ids == "0")) {
    stop(
      "`data` column `",
      group,
      "` cannot contain missing, empty, or `0` IDs.",
      call. = FALSE
    )
  }

  observed_missing <- setdiff(unique(observed_ids), ids)
  if (length(observed_missing) > 0L) {
    stop(
      "`data` column `",
      group,
      "` contains ID",
      if (length(observed_missing) > 1L) "s" else "",
      " not present in `pedigree`: ",
      paste(observed_missing, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  normalized <- hs_topological_pedigree(ids, sire, dam)

  structure(
    list(
      data = normalized$data,
      columns = cols,
      ids = normalized$data$id,
      observed_ids = unique(observed_ids),
      parent_index = normalized$parent_index,
      original_order = normalized$original_order
    ),
    class = "hs_pedigree_spec"
  )
}

hs_pedigree_columns <- function(pedigree) {
  nm <- names(pedigree)
  pick <- function(candidates) {
    hit <- which(tolower(nm) %in% candidates)
    if (length(hit) == 0L) {
      return(NA_integer_)
    }
    hit[[1L]]
  }

  cols <- list(
    id = pick(c("id", "animal")),
    sire = pick(c("sire", "father")),
    dam = pick(c("dam", "mother"))
  )

  if (anyNA(unlist(cols, use.names = FALSE))) {
    cols <- list(id = 1L, sire = 2L, dam = 3L)
  }

  cols
}

hs_eval_pedigree <- function(expr, data, env) {
  tryCatch(
    eval(expr, envir = data, enclos = env),
    error = function(err) {
      stop(
        "Could not evaluate `pedigree = ",
        hs_deparse(expr),
        "`. ",
        "Provide a pedigree data frame in the formula environment.",
        call. = FALSE
      )
    }
  )
}

hs_is_genomic_primary_call <- function(expr) {
  expr <- hs_unwrap_parentheses(expr)
  hs_is_call(expr, "genomic")
}

hs_is_single_step_primary_call <- function(expr) {
  expr <- hs_unwrap_parentheses(expr)
  hs_is_call(expr, "single_step")
}

# Opt-in, EXPERIMENTAL supplied-relationship primary effects that let a user
# bring their own relationship matrix for the animal-model random effect.
# `relmat(1 | id, K = K)` supplies a dense symmetric positive-definite
# relationship/covariance matrix K (the parser marshals the inverse Kinv =
# solve(K)); `precision(1 | id, Q = Q)` supplies the precision (inverse)
# directly. Both are parsed as supplied-relationship-inverse primary effects and
# fit through the SAME covered engine path as genomic()'s Ginv.
hs_is_relmat_primary_call <- function(expr) {
  expr <- hs_unwrap_parentheses(expr)
  hs_is_call(expr, "relmat")
}

hs_is_precision_primary_call <- function(expr) {
  expr <- hs_unwrap_parentheses(expr)
  hs_is_call(expr, "precision")
}

# The opt-in primary effects whose relationship is a user-supplied inverse:
# `genomic()` (a genomic relationship inverse `Ginv`), `single_step()` (a
# single-step relationship inverse `Hinv`), and the experimental
# `relmat()`/`precision()` supplied-relationship terms (a supplied covariance
# inverted to `Kinv`, or a supplied precision). All fit by REML on a
# relationship-inverse-based animal_model_spec.
hs_is_relinv_primary_call <- function(expr) {
  hs_is_genomic_primary_call(expr) ||
    hs_is_single_step_primary_call(expr) ||
    hs_is_relmat_primary_call(expr) ||
    hs_is_precision_primary_call(expr)
}

hs_parse_relinv_primary_call <- function(call, data, env, model_data = NULL) {
  call <- hs_unwrap_parentheses(call)
  # Experimental supplied-relationship terms take their own parse: `relmat()`
  # supplies a covariance K (inverted to Kinv), `precision()` supplies the
  # inverse Q directly. Both return a supplied-inverse primary spec.
  if (hs_is_relmat_primary_call(call) || hs_is_precision_primary_call(call)) {
    return(hs_parse_supplied_relmat_call(call, data, env))
  }
  if (hs_is_genomic_primary_call(call)) {
    term <- "genomic"
    arg_name <- "Ginv"
  } else {
    term <- "single_step"
    arg_name <- "Hinv"
  }
  example <- sprintf("`%s(1 | id, %s = %s)`", term, arg_name, arg_name)

  args <- as.list(call)[-1L]
  arg_names <- names(args)
  if (is.null(arg_names)) {
    arg_names <- rep("", length(args))
  }

  bar_candidates <- which(arg_names == "" | arg_names == "formula")
  if (length(bar_candidates) != 1L) {
    stop(
      "`",
      term,
      "()` must have one random-effect expression, for example ",
      example,
      ".",
      call. = FALSE
    )
  }

  bar <- hs_unwrap_parentheses(args[[bar_candidates]])
  if (!hs_is_call(bar, "|") || length(bar) != 3L) {
    stop(
      "The first `",
      term,
      "()` argument must be a random-effect expression such as `1 | id`.",
      call. = FALSE
    )
  }

  lhs <- hs_unwrap_parentheses(bar[[2L]])
  group_expr <- hs_unwrap_parentheses(bar[[3L]])
  if (!hs_is_one(lhs)) {
    hs_abort_unsupported_syntax(
      "Only random-intercept syntax ",
      example,
      " is implemented. ",
      term,
      " slopes are planned, not implemented.",
      call. = FALSE
    )
  }
  if (!is.symbol(group_expr)) {
    stop(
      "The grouping variable in `",
      term,
      "()` must be a bare column name.",
      call. = FALSE
    )
  }

  group <- as.character(group_expr)
  if (!group %in% names(data)) {
    stop(
      "`",
      term,
      "()` grouping variable `",
      group,
      "` was not found in `data`.",
      call. = FALSE
    )
  }

  # single_step() construction path: `pedigree = ped` + `markers = M` builds the
  # single-step H^-1 engine-side (Aguilar et al. 2010), distinct from the
  # supplied-`Hinv` path below. With an `hs_data()` container that bundles a
  # pedigree and genotypes, `single_step(1 | id)` resolves both from the bundle
  # (the `animal(1 | id)` precedent), so neither argument is required.
  if (identical(term, "single_step")) {
    has_explicit_construct <- any(c("pedigree", "markers") %in% arg_names)
    has_hinv <- any(c("Hinv", "H") %in% arg_names)
    bundle_construct <- !has_hinv &&
      isTRUE(model_data$is_hs_data) &&
      !is.null(model_data$components$pedigree) &&
      !is.null(model_data$components$genotypes)
    if (has_explicit_construct || bundle_construct) {
      return(hs_parse_single_step_construct(
        call,
        args,
        arg_names,
        group,
        data,
        env,
        model_data = model_data
      ))
    }
  }

  # `genomic()` accepts either a supplied `Ginv` or a marker matrix `markers`
  # (the engine builds the genomic relationship from markers); `single_step()`
  # accepts a supplied `Hinv`.
  accepted <- if (identical(term, "genomic")) c("Ginv", "markers") else arg_name
  named_args <- args[arg_names != ""]
  unsupported <- setdiff(names(named_args), accepted)
  if (length(unsupported) > 0L) {
    stop(
      "`",
      term,
      "()` argument",
      if (length(unsupported) > 1L) "s " else " ",
      paste(sprintf("`%s`", unsupported), collapse = ", "),
      " are planned, not implemented.",
      call. = FALSE
    )
  }
  supplied <- intersect(accepted, names(named_args))
  if (length(supplied) == 0L) {
    extra <- if (identical(term, "genomic")) {
      " (or `markers` to build one)"
    } else {
      # single_step: also point at the construction on-ramps (explicit args or an
      # hs_data() bundle that carries BOTH a pedigree and genotypes), so a bare
      # single_step(1 | id) or a partial bundle is not left at a Hinv-only message.
      paste0(
        ", or `pedigree = ped` + `markers = M` (or an `hs_data()` bundle with ",
        "both a pedigree and genotypes) to construct one"
      )
    }
    stop(
      "`",
      term,
      "()` requires a `",
      arg_name,
      "` argument (a relationship inverse with id dimnames)",
      extra,
      ".",
      call. = FALSE
    )
  }
  if (length(supplied) > 1L) {
    stop(
      "`genomic()` takes exactly one of `Ginv` or `markers`.",
      call. = FALSE
    )
  }
  arg_used <- supplied[[1L]]

  observed_ids <- as.character(data[[group]])

  if (identical(arg_used, "markers")) {
    markers <- hs_eval_genomic_ginv(
      named_args$markers,
      data,
      env,
      what = "markers"
    )
    markers <- hs_validate_genomic_markers(markers)
    ids <- rownames(markers)
    unknown <- setdiff(unique(observed_ids), ids)
    if (length(unknown) > 0L) {
      shown <- unknown[seq_len(min(5L, length(unknown)))]
      stop(
        "`genomic()` ids must be present in the `markers` row names. ID(s) not ",
        "in the `markers`: ",
        paste(sprintf("`%s`", shown), collapse = ", "),
        ".",
        call. = FALSE
      )
    }
    return(list(
      type = term,
      term = hs_deparse(call),
      design = "intercept",
      group = group,
      values = observed_ids,
      ids = ids,
      markers = markers,
      source = "markers",
      relationship = term,
      covariance = "scalar"
    ))
  }

  relinv <- hs_eval_genomic_ginv(
    named_args[[arg_used]],
    data,
    env,
    what = arg_used
  )
  relinv <- hs_validate_genomic_ginv(relinv, what = arg_used)
  relinv_ids <- rownames(relinv)
  unknown <- setdiff(unique(observed_ids), relinv_ids)
  if (length(unknown) > 0L) {
    shown <- unknown[seq_len(min(5L, length(unknown)))]
    stop(
      "`",
      term,
      "()` ids must be present in the `",
      arg_used,
      "` dimnames. ID(s) not in the `",
      arg_used,
      "`: ",
      paste(sprintf("`%s`", shown), collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  list(
    type = term,
    term = hs_deparse(call),
    design = "intercept",
    group = group,
    values = observed_ids,
    ids = relinv_ids,
    ginv = relinv,
    source = "supplied",
    relationship = term,
    covariance = "scalar"
  )
}

# Coerce an `hs_data()` genotypes component into the numeric dosage matrix (rows =
# genotyped ids) that single-step construction expects, for the
# `single_step(1 | id)` bundle shorthand. Mirrors the id resolution `hs_data()`
# itself enforces at construction: a matrix carries ids as row names; a data frame
# uses the bundle's `id` column or explicit row names. The error branches are
# defensive (a valid `hs_data()` already satisfies one of these forms).
hs_single_step_bundle_markers <- function(genotypes, id) {
  if (is.matrix(genotypes)) {
    return(genotypes)
  }
  if (is.data.frame(genotypes)) {
    if (!is.null(id) && id %in% names(genotypes)) {
      ids <- as.character(genotypes[[id]])
      mat <- as.matrix(genotypes[setdiff(names(genotypes), id)])
      rownames(mat) <- ids
      return(mat)
    }
    if (hs_has_explicit_rownames(genotypes)) {
      return(as.matrix(genotypes))
    }
    stop(
      "`hs_data()` `genotypes` must carry genotyped ids as a `",
      id %||% "id",
      "` column or row names to use the `single_step(1 | id)` shorthand; ",
      "otherwise pass `markers =` explicitly.",
      call. = FALSE
    )
  }
  stop(
    "`hs_data()` `genotypes` must be a numeric matrix or data frame to use the ",
    "`single_step(1 | id)` shorthand.",
    call. = FALSE
  )
}

# Parse the single-step CONSTRUCTION call `single_step(1 | id, pedigree = ped,
# markers = M, ...)`. Combines the pedigree path (id/sire/dam -> Ainv, A) with a
# genotyped-subset marker matrix (-> G among the genotyped animals). The crux is
# the `genotyped_rows` alignment (docs/design/25 section 3): the genotyped ids must be a
# subset of the pedigree ids (NOT the observed ids -- ungenotyped phenotyped
# animals are the whole point of single-step), and the marker rows are reordered
# to the genotyped animals' sorted pedigree-row positions so G aligns with the
# `genotyped_rows` the engine expects.
hs_parse_single_step_construct <- function(
  call,
  args,
  arg_names,
  group,
  data,
  env,
  model_data = NULL
) {
  named_args <- args[arg_names != ""]
  # The construction path and the supplied-inverse path are mutually exclusive.
  if (any(c("Hinv", "H") %in% names(named_args))) {
    stop(
      "`single_step()` takes EITHER a precomputed `Hinv` ",
      "(`single_step(1 | id, Hinv = Hinv)`) OR `pedigree` + `markers` to ",
      "construct H^-1 -- not both.",
      call. = FALSE
    )
  }
  accepted <- c(
    "pedigree",
    "markers",
    "group",
    "Gamma",
    "tau",
    "omega",
    "blend_weight",
    "ridge"
  )
  unsupported <- setdiff(names(named_args), accepted)
  if (length(unsupported) > 0L) {
    stop(
      "`single_step()` construction argument",
      if (length(unsupported) > 1L) "s " else " ",
      paste(sprintf("`%s`", unsupported), collapse = ", "),
      if (length(unsupported) > 1L) " are" else " is",
      " planned, not implemented.",
      call. = FALSE
    )
  }
  components <- model_data$components
  is_bundle <- isTRUE(model_data$is_hs_data)
  has_bundle_ped <- is_bundle && !is.null(components$pedigree)
  has_bundle_geno <- is_bundle && !is.null(components$genotypes)

  # Pedigree: an explicit `pedigree =` wins; otherwise the `hs_data()` bundle.
  if (!is.null(named_args$pedigree)) {
    ped_df <- hs_eval_pedigree(named_args$pedigree, data, env)
  } else if (has_bundle_ped) {
    ped_df <- components$pedigree
  } else {
    stop(
      "`single_step()` construction needs a pedigree (an explicit ",
      "`pedigree = ped` or an `hs_data()` bundle with one): ",
      "`single_step(1 | id, pedigree = ped, markers = M)`. For a precomputed ",
      "inverse, write `single_step(1 | id, Hinv = Hinv)`.",
      call. = FALSE
    )
  }

  # Markers: an explicit `markers =` wins; otherwise the `hs_data()` genotypes.
  if (!is.null(named_args$markers)) {
    markers <- hs_eval_genomic_ginv(
      named_args$markers,
      data,
      env,
      what = "markers"
    )
  } else if (has_bundle_geno) {
    markers <- hs_single_step_bundle_markers(
      components$genotypes,
      model_data$id
    )
  } else {
    stop(
      "`single_step(1 | id, pedigree = ped, markers = M)` construction requires ",
      "a `markers` matrix alongside `pedigree` (or an `hs_data()` bundle with ",
      "genotypes). Supply a precomputed `Hinv` instead via ",
      "`single_step(1 | id, Hinv = Hinv)`.",
      call. = FALSE
    )
  }
  markers <- hs_validate_genomic_markers(markers)

  observed_ids <- as.character(data[[group]])
  pedigree <- hs_validate_pedigree(ped_df, observed_ids, group)
  ped_ids <- pedigree$ids
  geno_ids <- rownames(markers)
  not_in_ped <- setdiff(geno_ids, ped_ids)
  if (length(not_in_ped) > 0L) {
    shown <- not_in_ped[seq_len(min(5L, length(not_in_ped)))]
    stop(
      "`single_step()` `markers` row ids must all be present in the pedigree. ",
      "Genotyped id(s) not in the pedigree: ",
      paste(sprintf("`%s`", shown), collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  # genotyped_rows: the genotyped animals' pedigree-row positions, sorted
  # ascending; reorder the marker rows to match so G's rows/cols align.
  ped_pos <- match(geno_ids, ped_ids)
  ord <- order(ped_pos)
  genotyped_rows <- ped_pos[ord]
  markers <- markers[ord, , drop = FALSE]

  has_metafounder_group <- "group" %in% names(named_args)
  has_gamma <- "Gamma" %in% names(named_args)
  if (xor(has_metafounder_group, has_gamma)) {
    stop(
      "`single_step()` H^Gamma construction requires both `group` and ",
      "`Gamma` together: `single_step(1 | id, pedigree = ped, markers = M, ",
      "group = mf_group, Gamma = Gamma)`.",
      call. = FALSE
    )
  }

  knob <- function(name, default) {
    if (is.null(named_args[[name]])) {
      return(default)
    }
    value <- as.numeric(eval(named_args[[name]], envir = data, enclos = env))
    if (length(value) != 1L || !is.finite(value)) {
      stop(
        "`single_step()` `",
        name,
        "` must be a single finite number.",
        call. = FALSE
      )
    }
    value
  }

  group_of <- NULL
  Gamma <- NULL
  gamma_labels <- NULL
  relationship <- "single_step"
  source <- "construct"
  if (has_metafounder_group) {
    group_of <- hs_eval_metafounder_group(named_args$group, data, env)
    group_of <- hs_validate_metafounder_group(group_of, pedigree)
    gamma <- hs_eval_metafounder_gamma(named_args$Gamma, data, env)
    gamma <- hs_validate_metafounder_gamma(gamma, group_of, pedigree)
    Gamma <- gamma$Gamma
    gamma_labels <- gamma$labels
    relationship <- "metafounder_single_step"
    source <- "metafounder_construct"
  }

  list(
    type = "single_step",
    source = source,
    term = hs_deparse(call),
    design = "intercept",
    group = group,
    values = observed_ids,
    ids = ped_ids,
    pedigree = pedigree,
    markers = markers,
    genotyped_rows = genotyped_rows,
    group_of = group_of,
    Gamma = Gamma,
    gamma_labels = gamma_labels,
    tau = knob("tau", 1),
    omega = knob("omega", 1),
    blend_weight = knob("blend_weight", 0),
    ridge = knob("ridge", 0),
    relationship = relationship,
    covariance = "scalar"
  )
}

hs_eval_metafounder_group <- function(
  expr,
  data,
  env,
  term = "single_step"
) {
  label <- paste0("`", term, "()`")
  tryCatch(
    eval(expr, envir = data, enclos = env),
    error = function(err) {
      stop(
        "Could not evaluate ",
        label,
        " metafounder `group = ",
        hs_deparse(expr),
        "`. Provide an ID-named vector of metafounder group labels.",
        call. = FALSE
      )
    }
  )
}

hs_eval_metafounder_gamma <- function(
  expr,
  data,
  env,
  term = "single_step"
) {
  label <- paste0("`", term, "()`")
  tryCatch(
    eval(expr, envir = data, enclos = env),
    error = function(err) {
      stop(
        "Could not evaluate ",
        label,
        " `Gamma = ",
        hs_deparse(expr),
        "`. Provide a supplied metafounder relationship matrix.",
        call. = FALSE
      )
    }
  )
}

hs_validate_metafounder_group <- function(
  group_of,
  pedigree,
  term = "single_step"
) {
  label <- paste0("`", term, "()`")
  ped_ids <- pedigree$ids
  if (is.data.frame(group_of)) {
    nms <- names(group_of)
    id_col <- intersect(c("id", "animal"), nms)
    group_col <- intersect(c("group", "metafounder", "mf_group"), nms)
    if (length(id_col) == 0L || length(group_col) == 0L) {
      stop(
        label,
        " metafounder `group` data frames must contain an ",
        "`id` (or `animal`) column and a `group` column.",
        call. = FALSE
      )
    }
    ids <- as.character(group_of[[id_col[[1L]]]])
    values <- group_of[[group_col[[1L]]]]
    names(values) <- ids
    group_of <- values
  }
  if (is.matrix(group_of) || length(dim(group_of)) > 1L) {
    stop(
      label,
      " metafounder `group` must be a vector or data frame, ",
      "not a matrix.",
      call. = FALSE
    )
  }
  values <- as.character(group_of)
  names_in <- names(group_of)
  if (is.null(names_in) || anyNA(names_in) || any(!nzchar(names_in))) {
    stop(
      label,
      " metafounder `group` must be named by pedigree IDs so ",
      "it can be reordered to normalized pedigree order.",
      call. = FALSE
    )
  }
  if (anyDuplicated(names_in)) {
    dup <- unique(names_in[duplicated(names_in)])[1L]
    stop(
      label,
      " metafounder `group` has duplicated id `",
      dup,
      "`.",
      call. = FALSE
    )
  }
  missing_ids <- setdiff(ped_ids, names_in)
  if (length(missing_ids) > 0L) {
    shown <- missing_ids[seq_len(min(5L, length(missing_ids)))]
    stop(
      label,
      " metafounder `group` is missing pedigree id(s): ",
      paste(sprintf("`%s`", shown), collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  extra_ids <- setdiff(names_in, ped_ids)
  if (length(extra_ids) > 0L) {
    shown <- extra_ids[seq_len(min(5L, length(extra_ids)))]
    stop(
      label,
      " metafounder `group` includes id(s) not in the ",
      "normalized pedigree: ",
      paste(sprintf("`%s`", shown), collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  aligned <- values[match(ped_ids, names_in)]
  needs <- is.na(pedigree$data$sire) | is.na(pedigree$data$dam)
  bad_needed <- needs & (is.na(aligned) | !nzchar(aligned))
  if (any(bad_needed)) {
    shown <- ped_ids[bad_needed][seq_len(min(5L, sum(bad_needed)))]
    stop(
      label,
      " metafounder `group` must provide non-missing labels ",
      "for animals with unknown parents. Missing label(s): ",
      paste(sprintf("`%s`", shown), collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  aligned[is.na(aligned)] <- ""
  names(aligned) <- ped_ids
  aligned
}

hs_validate_metafounder_gamma <- function(
  Gamma,
  group_of,
  pedigree,
  term = "single_step"
) {
  label <- paste0("`", term, "()`")
  Gamma <- as.matrix(Gamma)
  if (!is.numeric(Gamma)) {
    stop(label, " `Gamma` must be a numeric matrix.", call. = FALSE)
  }
  if (nrow(Gamma) != ncol(Gamma)) {
    stop(label, " `Gamma` must be square.", call. = FALSE)
  }
  if (!all(is.finite(Gamma))) {
    stop(
      label,
      " `Gamma` must contain only finite values.",
      call. = FALSE
    )
  }
  if (!isTRUE(all.equal(Gamma, t(Gamma), tolerance = 1e-10))) {
    stop(label, " `Gamma` must be symmetric.", call. = FALSE)
  }
  needs <- is.na(pedigree$data$sire) | is.na(pedigree$data$dam)
  labels <- unique(unname(group_of[needs]))
  labels <- labels[nzchar(labels)]
  if (length(labels) == 0L) {
    stop(
      label,
      " metafounder `group` resolved no unknown-parent ",
      "metafounder labels.",
      call. = FALSE
    )
  }
  if (nrow(Gamma) != length(labels)) {
    stop(
      label,
      " `Gamma` must have one row/column per resolved ",
      "metafounder group (",
      length(labels),
      "); got ",
      nrow(Gamma),
      ".",
      call. = FALSE
    )
  }
  rn <- rownames(Gamma)
  cn <- colnames(Gamma)
  if (!is.null(rn) || !is.null(cn)) {
    if (is.null(rn) || is.null(cn) || !identical(rn, cn)) {
      stop(
        label,
        " `Gamma` row names and column names must match.",
        call. = FALSE
      )
    }
    if (anyDuplicated(rn)) {
      stop(label, " `Gamma` group names must be unique.", call. = FALSE)
    }
    missing_labels <- setdiff(labels, rn)
    extra_labels <- setdiff(rn, labels)
    if (length(missing_labels) > 0L || length(extra_labels) > 0L) {
      stop(
        label,
        " `Gamma` dimnames must match the resolved ",
        "metafounder group labels.",
        call. = FALSE
      )
    }
    Gamma <- Gamma[labels, labels, drop = FALSE]
  }
  ev <- eigen(Gamma, symmetric = TRUE, only.values = TRUE)$values
  scale <- max(1, max(abs(Gamma)))
  if (min(ev) < -1e-10 * scale) {
    stop(
      label,
      " `Gamma` must be positive semidefinite.",
      call. = FALSE
    )
  }
  list(Gamma = unname(Gamma), labels = labels)
}

hs_eval_genomic_ginv <- function(expr, data, env, what = "Ginv") {
  tryCatch(
    eval(expr, envir = data, enclos = env),
    error = function(err) {
      stop(
        "Could not evaluate the `",
        what,
        "` argument (`",
        hs_deparse(expr),
        "`). Provide it as a matrix in the formula environment.",
        call. = FALSE
      )
    }
  )
}

hs_validate_genomic_ginv <- function(ginv, what = "Ginv") {
  if (inherits(ginv, "Matrix")) {
    ginv <- as.matrix(ginv)
  }
  if (!is.matrix(ginv) || !is.numeric(ginv)) {
    stop("`", what, "` must be a numeric matrix.", call. = FALSE)
  }
  if (nrow(ginv) < 1L || ncol(ginv) < 1L || nrow(ginv) != ncol(ginv)) {
    stop("`", what, "` must be a nonempty square matrix.", call. = FALSE)
  }
  if (any(!is.finite(ginv))) {
    stop(
      "`",
      what,
      "` must contain only finite, non-missing values.",
      call. = FALSE
    )
  }
  ids <- rownames(ginv)
  col_ids <- colnames(ginv)
  invalid_ids <- function(x) {
    is.null(x) || any(is.na(x)) || any(!nzchar(x)) || anyDuplicated(x) > 0L
  }
  if (invalid_ids(ids) || invalid_ids(col_ids)) {
    stop(
      "`",
      what,
      "` must have unique, non-missing, nonempty row/column names matching the ",
      "genomic ids.",
      call. = FALSE
    )
  }
  if (!identical(ids, col_ids)) {
    stop(
      "`",
      what,
      "` row and column names must match in order.",
      call. = FALSE
    )
  }
  scale <- max(1, max(abs(ginv)))
  if (max(abs(ginv - t(ginv))) > 1e-10 * scale) {
    stop("`", what, "` must be symmetric.", call. = FALSE)
  }
  if (inherits(try(chol(ginv), silent = TRUE), "try-error")) {
    stop("`", what, "` must be positive definite.", call. = FALSE)
  }
  ginv
}

hs_validate_genomic_markers <- function(markers) {
  if (inherits(markers, "Matrix")) {
    markers <- as.matrix(markers)
  }
  if (!is.matrix(markers) || !is.numeric(markers)) {
    stop("`markers` must be a numeric marker matrix.", call. = FALSE)
  }
  if (nrow(markers) < 1L || ncol(markers) < 1L) {
    stop(
      "`markers` must have at least one row and one marker column.",
      call. = FALSE
    )
  }
  if (any(!is.finite(markers))) {
    stop(
      "`markers` must contain only finite, non-missing dosages.",
      call. = FALSE
    )
  }
  if (any(markers < 0 | markers > 2)) {
    hs_abort_out_of_range(
      "`markers` dosages must lie in the closed interval [0, 2]: the ",
      "required convention is the biallelic allele count 0/1/2 (or an ",
      "imputed dosage in [0, 2]). A matrix in the centered -1/0/1 ",
      "convention must be recoded (add 1) before it is passed in."
    )
  }
  ids <- rownames(markers)
  if (
    is.null(ids) ||
      any(is.na(ids)) ||
      any(!nzchar(ids)) ||
      anyDuplicated(ids) > 0L
  ) {
    stop(
      "`markers` must have unique, non-missing, nonempty row names matching the ids ",
      "(one row per genotyped individual).",
      call. = FALSE
    )
  }
  marker_names <- colnames(markers)
  if (
    !is.null(marker_names) &&
      (any(is.na(marker_names)) ||
        any(!nzchar(marker_names)) ||
        anyDuplicated(marker_names) > 0L)
  ) {
    stop(
      "`markers` column names, when supplied, must be unique, non-missing, and nonempty.",
      call. = FALSE
    )
  }
  p <- colMeans(markers) / 2
  k <- 2 * sum(p * (1 - p))
  if (!is.finite(k) || k <= 0) {
    stop(
      "`markers` must contain polymorphic variation: the sample-frequency ",
      "VanRaden denominator `k` must be positive.",
      call. = FALSE
    )
  }
  markers
}

# Parse an experimental supplied-relationship primary term. `relmat(1 | id,
# K = K)` supplies a relationship/covariance matrix K, marshalled as the inverse
# Kinv = solve(K); `precision(1 | id, Q = Q)` (also `Kinv =`) supplies the
# precision (inverse) directly. Mirrors the genomic supplied-`Ginv` parse: the
# supplied inverse is carried in `ginv` with `source = "supplied"`, so the
# existing supplied-relationship-inverse bridge/engine path fits it unchanged.
hs_parse_supplied_relmat_call <- function(call, data, env) {
  call <- hs_unwrap_parentheses(call)
  is_precision <- hs_is_precision_primary_call(call)
  term <- if (is_precision) "precision" else "relmat"
  # `relmat` supplies a covariance `K` (marshalled as `Kinv = solve(K)`) OR a
  # relationship inverse directly via `Kinv`/`Q` (the precision spelling, no
  # solve); `precision` supplies the inverse `Q` (or the `Kinv` spelling)
  # directly. `primary_arg` names the canonical argument for messages; `accepted`
  # lists every spelling the term takes.
  primary_arg <- if (is_precision) "Q" else "K"
  accepted <- if (is_precision) c("Q", "Kinv") else c("K", "Kinv", "Q")
  example <- sprintf("`%s(1 | id, %s = %s)`", term, primary_arg, primary_arg)

  args <- as.list(call)[-1L]
  arg_names <- names(args)
  if (is.null(arg_names)) {
    arg_names <- rep("", length(args))
  }

  bar_candidates <- which(arg_names == "" | arg_names == "formula")
  if (length(bar_candidates) != 1L) {
    stop(
      "`",
      term,
      "()` must have one random-effect expression, for example ",
      example,
      ".",
      call. = FALSE
    )
  }

  bar <- hs_unwrap_parentheses(args[[bar_candidates]])
  if (!hs_is_call(bar, "|") || length(bar) != 3L) {
    stop(
      "The first `",
      term,
      "()` argument must be a random-effect expression such as `1 | id`.",
      call. = FALSE
    )
  }

  lhs <- hs_unwrap_parentheses(bar[[2L]])
  group_expr <- hs_unwrap_parentheses(bar[[3L]])
  if (!hs_is_one(lhs)) {
    stop(
      "Only random-intercept syntax ",
      example,
      " is implemented. ",
      term,
      " slopes are planned, not implemented.",
      call. = FALSE
    )
  }
  if (!is.symbol(group_expr)) {
    stop(
      "The grouping variable in `",
      term,
      "()` must be a bare column name.",
      call. = FALSE
    )
  }

  group <- as.character(group_expr)
  if (!group %in% names(data)) {
    stop(
      "`",
      term,
      "()` grouping variable `",
      group,
      "` was not found in `data`.",
      call. = FALSE
    )
  }

  named_args <- args[arg_names != ""]
  unsupported <- setdiff(names(named_args), accepted)
  if (length(unsupported) > 0L) {
    stop(
      "`",
      term,
      "()` argument",
      if (length(unsupported) > 1L) "s " else " ",
      paste(sprintf("`%s`", unsupported), collapse = ", "),
      if (length(unsupported) > 1L) " are " else " is ",
      "planned, not implemented.",
      call. = FALSE
    )
  }
  supplied <- intersect(accepted, names(named_args))
  if (length(supplied) == 0L) {
    stop(
      "`",
      term,
      "()` requires a `",
      primary_arg,
      "` argument (a ",
      if (is_precision) "precision (inverse) matrix" else "relationship matrix",
      " with id row/column names).",
      call. = FALSE
    )
  }
  if (length(supplied) > 1L) {
    stop(
      "`",
      term,
      "()` takes exactly one of ",
      paste(sprintf("`%s`", accepted), collapse = " or "),
      ".",
      call. = FALSE
    )
  }
  arg_used <- supplied[[1L]]

  supplied_matrix <- hs_eval_genomic_ginv(
    named_args[[arg_used]],
    data,
    env,
    what = arg_used
  )
  supplied_matrix <- hs_validate_relmat_matrix(supplied_matrix, arg_used)

  # A covariance is inverted to its precision (`Kinv = solve(K)`); a directly
  # supplied inverse (`precision()`'s `Q`, or `relmat()`/`precision()`'s `Kinv`,
  # or `relmat()`'s `Q`) is carried verbatim, no solve. Only the covariance
  # spelling `K` needs inverting. Both travel as `ginv` (the supplied
  # relationship inverse) so the shared supplied-inverse bridge path is reused.
  relinv <- if (identical(arg_used, "K")) {
    solve(supplied_matrix)
  } else {
    supplied_matrix
  }
  relinv_ids <- rownames(supplied_matrix)

  observed_ids <- as.character(data[[group]])
  unknown <- setdiff(unique(observed_ids), relinv_ids)
  if (length(unknown) > 0L) {
    shown <- unknown[seq_len(min(5L, length(unknown)))]
    stop(
      "`",
      term,
      "()` ids must be present in the `",
      arg_used,
      "` dimnames. ID(s) not in the `",
      arg_used,
      "`: ",
      paste(sprintf("`%s`", shown), collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  list(
    type = term,
    term = hs_deparse(call),
    design = "intercept",
    group = group,
    values = observed_ids,
    ids = relinv_ids,
    ginv = relinv,
    source = "supplied",
    relationship = term,
    covariance = "scalar"
  )
}

# Validate a supplied `relmat()`/`precision()` matrix: numeric, square, finite,
# symmetric, positive-definite, and ID-keyed via matching row/column names. These
# are data-validation errors (a different family from unsupported-syntax errors).
hs_validate_relmat_matrix <- function(mat, arg_name) {
  if (inherits(mat, "Matrix")) {
    mat <- as.matrix(mat)
  }
  if (!is.matrix(mat) || !is.numeric(mat)) {
    stop("`", arg_name, "` must be a numeric matrix.", call. = FALSE)
  }
  if (nrow(mat) != ncol(mat)) {
    stop("`", arg_name, "` must be a square matrix.", call. = FALSE)
  }
  if (any(!is.finite(mat))) {
    stop("`", arg_name, "` must contain only finite values.", call. = FALSE)
  }
  ids <- rownames(mat)
  if (is.null(ids) || any(is.na(ids)) || anyDuplicated(ids) > 0L) {
    stop(
      "`",
      arg_name,
      "` must have unique, non-missing row/column names matching the ids.",
      call. = FALSE
    )
  }
  if (!identical(ids, colnames(mat))) {
    stop("`", arg_name, "` row and column names must match.", call. = FALSE)
  }
  if (!isSymmetric(unname(mat))) {
    stop("`", arg_name, "` must be symmetric.", call. = FALSE)
  }
  # Positive definiteness via a Cholesky factorization (fails on any
  # non-positive eigenvalue), so the supplied relationship inverse is a valid
  # precision for the REML animal-model spec.
  chol_ok <- tryCatch(
    {
      chol(mat)
      TRUE
    },
    error = function(err) FALSE
  )
  if (!chol_ok) {
    stop("`", arg_name, "` must be positive definite.", call. = FALSE)
  }
  mat
}

hs_normalize_parent <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == "" | x == "0"] <- NA_character_
  x
}

hs_topological_pedigree <- function(ids, sire, dam) {
  sire_index <- match(sire, ids)
  dam_index <- match(dam, ids)
  sire_index[is.na(sire_index)] <- 0L
  dam_index[is.na(dam_index)] <- 0L

  # Iterative DFS post-order (explicit stack) emulating the former recursion so
  # deep linear pedigrees no longer overflow the call stack or get misreported
  # as cycles. The ordering is identical to the recursive version: for each node
  # we fully process the sire subtree, then the dam subtree, then emit the node.
  # `state` codes 0 = unseen (white), 1 = on the stack (grey), 2 = done (black);
  # re-entering a grey node is a real parent-offspring cycle.
  n <- length(ids)
  state <- integer(n)
  order <- integer(n)
  n_emitted <- 0L

  # Each stack frame is a node `index` with a `phase`: 0 = on entry (mark grey,
  # then descend to sire), 1 = sire subtree done (descend to dam), 2 = dam
  # subtree done (mark black, emit). Frames are stored in parallel vectors with
  # a top pointer to avoid per-push reallocation on deep pedigrees.
  cap <- max(n, 1L)
  stack_index <- integer(cap)
  stack_phase <- integer(cap)
  top <- 0L

  push <- function(index, phase) {
    top <<- top + 1L
    if (top > length(stack_index)) {
      new_cap <- length(stack_index) * 2L
      stack_index <<- c(stack_index, integer(new_cap - length(stack_index)))
      stack_phase <<- c(stack_phase, integer(new_cap - length(stack_phase)))
    }
    stack_index[[top]] <<- index
    stack_phase[[top]] <<- phase
    invisible(NULL)
  }

  for (start in seq_len(n)) {
    if (state[[start]] != 0L) {
      next
    }
    push(start, 0L)
    while (top > 0L) {
      index <- stack_index[[top]]
      phase <- stack_phase[[top]]

      if (phase == 0L) {
        # Entry: a grey node here means we re-entered a node still on the stack,
        # i.e. a true parent-offspring cycle. A black node was already emitted.
        if (state[[index]] == 2L) {
          top <- top - 1L
          next
        }
        if (state[[index]] == 1L) {
          stop(
            "`pedigree` contains a parent-offspring cycle involving ID `",
            ids[[index]],
            "`.",
            call. = FALSE
          )
        }
        state[[index]] <- 1L
        stack_phase[[top]] <- 1L
        if (sire_index[[index]] != 0L) {
          push(sire_index[[index]], 0L)
        }
      } else if (phase == 1L) {
        stack_phase[[top]] <- 2L
        if (dam_index[[index]] != 0L) {
          push(dam_index[[index]], 0L)
        }
      } else {
        state[[index]] <- 2L
        n_emitted <- n_emitted + 1L
        order[[n_emitted]] <- index
        top <- top - 1L
      }
    }
  }

  sorted_position <- integer(length(ids))
  sorted_position[order] <- seq_along(order)
  sorted_sire_index <- sire_index[order]
  sorted_dam_index <- dam_index[order]
  known_sire <- sorted_sire_index != 0L
  known_dam <- sorted_dam_index != 0L
  sorted_sire_index[known_sire] <- sorted_position[sorted_sire_index[
    known_sire
  ]]
  sorted_dam_index[known_dam] <- sorted_position[sorted_dam_index[known_dam]]

  list(
    data = data.frame(
      id = ids[order],
      sire = sire[order],
      dam = dam[order],
      stringsAsFactors = FALSE
    ),
    parent_index = list(
      sire = sorted_sire_index,
      dam = sorted_dam_index
    ),
    original_order = order
  )
}

hs_split_additive_rhs <- function(expr) {
  expr <- hs_unwrap_parentheses(expr)
  if (hs_is_call(expr, "+") && length(expr) == 3L) {
    return(c(
      hs_split_additive_rhs(expr[[2L]]),
      hs_split_additive_rhs(expr[[3L]])
    ))
  }
  list(expr)
}

hs_rebuild_additive_rhs <- function(terms) {
  if (length(terms) == 0L) {
    return(1)
  }
  Reduce(function(left, right) call("+", left, right), terms)
}

hs_is_animal_call <- function(expr) {
  expr <- hs_unwrap_parentheses(expr)
  hs_is_call(expr, "animal")
}

hs_is_metafounder_primary_call <- function(expr) {
  expr <- hs_unwrap_parentheses(expr)
  hs_is_call(expr, "metafounder")
}

hs_is_permanent_call <- function(expr) {
  expr <- hs_unwrap_parentheses(expr)
  hs_is_call(expr, "permanent")
}

hs_is_common_env_call <- function(expr) {
  expr <- hs_unwrap_parentheses(expr)
  hs_is_call(expr, "common_env")
}

hs_is_maternal_genetic_call <- function(expr) {
  expr <- hs_unwrap_parentheses(expr)
  hs_is_call(expr, "maternal_genetic")
}

# The planned QG markers the parser now consumes as the (single) second random
# effect of an opt-in two-effect model: `permanent()` (repeatability),
# `common_env()` (IID common-environment), and `maternal_genetic()` (a
# pedigree-related maternal genetic effect).
hs_is_second_effect_call <- function(expr) {
  hs_is_permanent_call(expr) ||
    hs_is_common_env_call(expr) ||
    hs_is_maternal_genetic_call(expr)
}

hs_parse_second_effect_call <- function(call, data, animal_spec) {
  if (hs_is_permanent_call(call)) {
    hs_parse_permanent_call(call, data, animal_spec)
  } else if (hs_is_common_env_call(call)) {
    hs_parse_common_env_call(call, data)
  } else {
    hs_parse_maternal_genetic_call(call, data, animal_spec)
  }
}

# Parse `common_env(1 | group)` as the common-environment effect of the opt-in
# two-effect model: a random intercept on an environmental grouping (e.g. litter
# or cage) carrying an identity relationship (each level an independent IID
# effect). Unlike `permanent()`, the grouping is a separate environmental column,
# not the animal id.
hs_parse_common_env_call <- function(call, data) {
  call <- hs_unwrap_parentheses(call)
  args <- as.list(call)[-1L]
  arg_names <- names(args)
  if (is.null(arg_names)) {
    arg_names <- rep("", length(args))
  }

  bar_candidates <- which(arg_names == "" | arg_names == "formula")
  if (length(bar_candidates) != 1L) {
    stop(
      "`common_env()` must have one random-effect expression, for example ",
      "`common_env(1 | litter)`.",
      call. = FALSE
    )
  }

  bar <- hs_unwrap_parentheses(args[[bar_candidates]])
  if (!hs_is_call(bar, "|") || length(bar) != 3L) {
    stop(
      "The `common_env()` argument must be a random-effect expression such as ",
      "`1 | litter`.",
      call. = FALSE
    )
  }

  lhs <- hs_unwrap_parentheses(bar[[2L]])
  group_expr <- hs_unwrap_parentheses(bar[[3L]])
  if (!hs_is_one(lhs)) {
    hs_abort_unsupported_syntax(
      "Only random-intercept syntax `common_env(1 | group)` is implemented. ",
      "Common-environment slopes are planned, not implemented.",
      call. = FALSE
    )
  }
  if (!is.symbol(group_expr)) {
    stop(
      "The grouping variable in `common_env()` must be a bare column name.",
      call. = FALSE
    )
  }

  named_args <- args[arg_names != ""]
  if (length(named_args) > 0L) {
    stop(
      "`common_env()` takes no extra arguments in the v0.1 two-effect model.",
      call. = FALSE
    )
  }

  group <- as.character(group_expr)
  if (!group %in% names(data)) {
    stop(
      "`common_env()` grouping variable `",
      group,
      "` was not found in `data`.",
      call. = FALSE
    )
  }

  list(
    type = "common_env",
    term = hs_deparse(call),
    design = "intercept",
    group = group,
    values = as.character(data[[group]]),
    levels = unique(as.character(data[[group]])),
    relationship = "identity",
    covariance = "scalar"
  )
}

# Data requirement plus covered-first next path for maternal_genetic()
# parser stops that fire before the default-path abort. Reuses
# hs_maternal_genetic_default_path_message() so order cannot drift:
# covered direct_maternal first, experimental two_effect second.
# Does not auto-route. Does not take dams from the pedigree table.
hs_maternal_genetic_parser_guidance <- function() {
  paste0(
    "The grouping column (usually `dam`) must be a column of `data`: ",
    "mothers of the records, with IDs already in the animal pedigree. ",
    "The parser does not take dams from the pedigree table and does not ",
    "accept `pedigree=` on `maternal_genetic()`.\n\n",
    hs_maternal_genetic_default_path_message()
  )
}

# Parse `maternal_genetic(1 | dam)` as the maternal genetic effect of the opt-in
# two-effect model: a random intercept expressed through the dam, carrying the
# SAME pedigree relationship as the direct animal effect (A2 = pedigree A). The
# grouping column holds dam ids, which must be animals in the `animal()`
# pedigree; the maternal effect is predicted for every pedigree animal.
hs_parse_maternal_genetic_call <- function(call, data, animal_spec) {
  call <- hs_unwrap_parentheses(call)
  args <- as.list(call)[-1L]
  arg_names <- names(args)
  if (is.null(arg_names)) {
    arg_names <- rep("", length(args))
  }

  bar_candidates <- which(arg_names == "" | arg_names == "formula")
  if (length(bar_candidates) != 1L) {
    stop(
      "`maternal_genetic()` must have one random-effect expression, for ",
      "example `maternal_genetic(1 | dam)`.",
      call. = FALSE
    )
  }

  bar <- hs_unwrap_parentheses(args[[bar_candidates]])
  if (!hs_is_call(bar, "|") || length(bar) != 3L) {
    stop(
      "The `maternal_genetic()` argument must be a random-effect expression ",
      "such as `1 | dam`.",
      call. = FALSE
    )
  }

  lhs <- hs_unwrap_parentheses(bar[[2L]])
  group_expr <- hs_unwrap_parentheses(bar[[3L]])
  if (!hs_is_one(lhs)) {
    hs_abort_unsupported_syntax(
      "Only random-intercept syntax `maternal_genetic(1 | dam)` is implemented. ",
      "Maternal slopes are planned, not implemented.",
      call. = FALSE
    )
  }
  if (!is.symbol(group_expr)) {
    stop(
      "The grouping variable in `maternal_genetic()` must be a bare column name.",
      call. = FALSE
    )
  }

  named_args <- args[arg_names != ""]
  if (length(named_args) > 0L) {
    hs_abort_unsupported_syntax(
      "`maternal_genetic()` takes no extra arguments. ",
      "Dam relationships come from the `animal()` pedigree, not a second ",
      "`pedigree=` argument.\n\n",
      hs_maternal_genetic_parser_guidance(),
      call. = FALSE
    )
  }

  group <- as.character(group_expr)
  if (!group %in% names(data)) {
    stop(
      "`maternal_genetic()` grouping variable `",
      group,
      "` was not found in `data`.\n\n",
      hs_maternal_genetic_parser_guidance(),
      call. = FALSE
    )
  }

  dam_values <- as.character(data[[group]])
  ped_ids <- animal_spec$pedigree$ids
  unknown <- setdiff(unique(dam_values), ped_ids)
  if (length(unknown) > 0L) {
    shown <- unknown[seq_len(min(5L, length(unknown)))]
    stop(
      "`maternal_genetic()` dams must be animals in the `animal()` pedigree. ",
      "Dam(s) not in the pedigree: ",
      paste(sprintf("`%s`", shown), collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  list(
    type = "maternal_genetic",
    term = hs_deparse(call),
    design = "intercept",
    group = group,
    values = dam_values,
    levels = ped_ids,
    relationship = "pedigree",
    covariance = "scalar"
  )
}

# Parse `permanent(1 | id)` as the permanent-environment effect of the opt-in
# repeatability model. It must be a random intercept on the SAME grouping
# variable as `animal()`, because the engine shares the animal incidence matrix
# `Z` between the additive-genetic and permanent-environment effects (the PE
# effect carries an identity relationship, A2 = I).
hs_parse_permanent_call <- function(call, data, animal_spec) {
  call <- hs_unwrap_parentheses(call)
  args <- as.list(call)[-1L]
  arg_names <- names(args)
  if (is.null(arg_names)) {
    arg_names <- rep("", length(args))
  }

  bar_candidates <- which(arg_names == "" | arg_names == "formula")
  if (length(bar_candidates) != 1L) {
    stop(
      "`permanent()` must have one random-effect expression, for example ",
      "`permanent(1 | id)`.",
      call. = FALSE
    )
  }

  bar <- hs_unwrap_parentheses(args[[bar_candidates]])
  if (!hs_is_call(bar, "|") || length(bar) != 3L) {
    stop(
      "The `permanent()` argument must be a random-effect expression such as ",
      "`1 | id`.",
      call. = FALSE
    )
  }

  lhs <- hs_unwrap_parentheses(bar[[2L]])
  group_expr <- hs_unwrap_parentheses(bar[[3L]])
  if (!hs_is_one(lhs)) {
    hs_abort_unsupported_syntax(
      "Only random-intercept syntax `permanent(1 | id)` is implemented. ",
      "Permanent-environment slopes are planned, not implemented.",
      call. = FALSE
    )
  }
  if (!is.symbol(group_expr)) {
    stop(
      "The grouping variable in `permanent()` must be a bare column name.",
      call. = FALSE
    )
  }

  named_args <- args[arg_names != ""]
  if (length(named_args) > 0L) {
    stop(
      "`permanent()` takes no extra arguments in the v0.1 repeatability model.",
      call. = FALSE
    )
  }

  group <- as.character(group_expr)
  if (!identical(group, animal_spec$group)) {
    stop(
      "`permanent()` must use the same grouping variable as `animal()` ",
      "(the permanent-environment effect shares the animal incidence). Got ",
      "`permanent(1 | ",
      group,
      ")` with `animal(1 | ",
      animal_spec$group,
      ")`.",
      call. = FALSE
    )
  }

  list(
    type = "permanent",
    term = hs_deparse(call),
    design = "intercept",
    group = group,
    values = as.character(data[[group]]),
    relationship = "identity",
    covariance = "scalar"
  )
}

hs_is_planned_marker_call <- function(expr) {
  expr <- hs_unwrap_parentheses(expr)
  is.call(expr) &&
    as.character(expr[[1L]]) %in% hs_planned_marker_names()
}

hs_planned_marker_names <- function() {
  c(
    hs_planned_genomic_marker_names(),
    hs_planned_qg_effect_marker_names()
  )
}

hs_planned_genomic_marker_names <- function() {
  c("genomic", "single_step", "markers", "marker_scan", "qtl_scan")
}

hs_planned_qg_effect_marker_names <- function() {
  c(
    "permanent",
    "common_env",
    "maternal_genetic",
    "maternal_env",
    "paternal_genetic",
    "paternal_env",
    "cytoplasmic",
    "imprinting",
    "dominance",
    "epistasis",
    "relmat",
    "precision",
    "group",
    "unknown_parent_group",
    "metafounder",
    "inbreeding"
  )
}

# Nearest live path for a reserved marker. common_env() / two_effect paste
# copy is owned elsewhere; this map does not mention that route.
hs_planned_marker_nearest_path <- function(marker) {
  switch(
    marker,
    maternal_env = paste0(
      "The live maternal sibling is `maternal_genetic(1 | dam)` with ",
      "`target = \"direct_maternal\"` (covered correlated 2x2 G / Willham ",
      "triple) or `target = \"two_effect\"` (experimental independent maternal)."
    ),
    paternal_genetic = paste0(
      "The live genetic effect is `animal(1 | id, pedigree = ped)`."
    ),
    paternal_env = paste0(
      "The live genetic effect is `animal(1 | id, pedigree = ped)`."
    ),
    inbreeding = paste0(
      "Pedigree inbreeding F is already used internally when building Ainv; ",
      "`inbreeding()` is a reserved future *effect*, not a missing coefficient."
    ),
    dominance = paste0(
      "The live 'bring your own kernel' path is `relmat(1 | id, K = )` ",
      "(or `precision(1 | id, Q = )`); that does not implement dominance."
    ),
    epistasis = paste0(
      "The live 'bring your own kernel' path is `relmat(1 | id, K = )` ",
      "(or `precision(1 | id, Q = )`); that does not implement epistasis."
    ),
    group = paste0(
      "The live grouping path is bare `(1 | group)` with ",
      "`target = \"multi_effect\"`."
    ),
    unknown_parent_group = paste0(
      "The live grouping path is bare `(1 | group)` with ",
      "`target = \"multi_effect\"`."
    ),
    cytoplasmic = paste0(
      "The live maternal genetic sibling is `maternal_genetic(1 | dam)`."
    ),
    imprinting = paste0(
      "The live maternal genetic sibling is `maternal_genetic(1 | dam)`."
    ),
    markers = paste0(
      "The live marker path is `genomic(1 | id, markers = )` or ",
      "`gwas(fit, markers)`."
    ),
    marker_scan = paste0(
      "The live marker path is `genomic(1 | id, markers = )` or ",
      "`gwas(fit, markers)`."
    ),
    qtl_scan = paste0(
      "The live marker path is `genomic(1 | id, markers = )` or ",
      "`gwas(fit, markers)`."
    ),
    ""
  )
}

hs_stop_planned_marker <- function(expr) {
  expr <- hs_unwrap_parentheses(expr)
  marker <- as.character(expr[[1L]])
  nearest <- hs_planned_marker_nearest_path(marker)
  hs_abort_unsupported_syntax(
    "`",
    marker,
    "()` is planned, not implemented.",
    if (nzchar(nearest)) paste0(" ", nearest) else "",
    " Run `formula_status()` for the live list of which terms parse and which fit.",
    call. = FALSE
  )
}

# Every recognized effect/marker name, whether implemented (`animal()`) or
# planned. Used to catch such a call nested inside an interaction or function
# term (e.g. `x:dominance(1 | id)`), which is not a top-level random effect.
hs_effect_marker_names <- function() {
  c("animal", hs_planned_marker_names())
}

# The first recognized effect/marker name used as a *call* anywhere inside an
# expression, or NA if none. Only call heads count, so a fixed variable that
# merely shares a name (e.g. a column called `dominance` in `sex:dominance`) is
# not flagged.
hs_nested_effect_call_name <- function(expr) {
  if (!is.call(expr)) {
    return(NA_character_)
  }
  head <- expr[[1L]]
  if (is.symbol(head) && as.character(head) %in% hs_effect_marker_names()) {
    return(as.character(head))
  }
  for (part in as.list(expr)[-1L]) {
    hit <- hs_nested_effect_call_name(part)
    if (!is.na(hit)) {
      return(hit)
    }
  }
  NA_character_
}

# Reject a recognized effect/marker call nested inside a fixed term. The term
# itself is not a top-level effect (those were already consumed), so any such
# call found here is nested and would otherwise abort with a cryptic base-R
# error when the fixed-effect model.frame is built. Ordinary fixed interactions
# such as `sex:age` are left untouched.
hs_check_nested_effect <- function(term) {
  marker <- hs_nested_effect_call_name(term)
  if (is.na(marker)) {
    return(invisible(NULL))
  }
  hs_abort_unsupported_syntax(
    "`",
    marker,
    "()` cannot appear inside an interaction or function call (`",
    hs_deparse(term),
    "`). Name it as a top-level random effect, for example ",
    "`y ~ ... + animal(1 | id, pedigree = ped)`. Run `formula_status()` for ",
    "the live list of which terms parse and which fit."
  )
}

hs_is_call <- function(expr, name) {
  is.call(expr) && identical(as.character(expr[[1L]]), name)
}

hs_is_one <- function(expr) {
  is.numeric(expr) && length(expr) == 1L && identical(as.numeric(expr), 1)
}

hs_unwrap_parentheses <- function(expr) {
  while (hs_is_call(expr, "(") && length(expr) == 2L) {
    expr <- expr[[2L]]
  }
  expr
}

hs_deparse <- function(expr) {
  paste(deparse(expr, width.cutoff = 500L), collapse = " ")
}

# A bare grouping/random-effect expression like `(1 | x)` or `x | id` -- i.e. a
# top-level `|` call that is not wrapped in a recognized effect function such as
# `animal()`/`permanent()`. These must be named, not silently absorbed into the
# fixed-effect design.
hs_is_bar_expr <- function(expr) {
  hs_is_call(hs_unwrap_parentheses(expr), "|")
}

# A correlated/uncorrelated slope term `(... || group)` (lme4's `||` syntax). The
# top-level operator is `||`, not `|`, so `hs_is_bar_expr()` does not catch it; it
# is detected separately and rejected, because hsquared has no random-slope /
# correlated-effect estimator (accepting one would over-claim).
hs_is_double_bar_expr <- function(expr) {
  hs_is_call(hs_unwrap_parentheses(expr), "||")
}

# A bare random INTERCEPT `(1 | group)`: a top-level `|` call whose LHS is the
# literal `1` and whose RHS is a bare grouping column name. This is the only bar
# form accepted as an opt-in i.i.d. random effect. A random slope `(x | group)`
# (LHS not `1`) is NOT accepted here.
hs_is_bare_iid_intercept <- function(expr) {
  bar <- hs_unwrap_parentheses(expr)
  if (!hs_is_call(bar, "|") || length(bar) != 3L) {
    return(FALSE)
  }
  lhs <- hs_unwrap_parentheses(bar[[2L]])
  hs_is_one(lhs)
}

# Parse a bare `(1 | group)` term into an i.i.d. random-intercept effect spec,
# mirroring `hs_parse_common_env_call()` (the SAME identity-relationship i.i.d.
# model). The grouping variable must be a bare column present in `data`.
hs_parse_bare_iid_call <- function(term, data) {
  bar <- hs_unwrap_parentheses(term)
  lhs <- hs_unwrap_parentheses(bar[[2L]])
  group_expr <- hs_unwrap_parentheses(bar[[3L]])
  if (!hs_is_one(lhs)) {
    hs_abort_unsupported_syntax(
      "Only random-intercept syntax `(1 | group)` is accepted as an i.i.d. ",
      "random effect. Random slopes are planned, not implemented.",
      call. = FALSE
    )
  }
  if (!is.symbol(group_expr)) {
    stop(
      "The grouping variable in a bare `(1 | group)` random effect must be a ",
      "bare column name.",
      call. = FALSE
    )
  }
  group <- as.character(group_expr)
  if (!group %in% names(data)) {
    stop(
      "The `(1 | ",
      group,
      ")` grouping variable `",
      group,
      "` was not found in `data`.",
      call. = FALSE
    )
  }
  list(
    type = "iid",
    term = hs_deparse(term),
    design = "intercept",
    group = group,
    values = as.character(data[[group]]),
    levels = unique(as.character(data[[group]])),
    relationship = "identity",
    covariance = "scalar"
  )
}

hs_stop_unsupported_random_effect <- function(term) {
  hs_abort_unsupported_syntax(
    "Unsupported random-effect term `",
    hs_deparse(term),
    "`. Only a bare random INTERCEPT `(1 | group)` is accepted as an i.i.d. ",
    "random effect; hsquared has no random-slope `(x | group)` or correlated ",
    "`(x || group)` estimator. Use a bare `(1 | group)` intercept, or name the ",
    "effect: `animal(1 | id, pedigree = ped)` for the additive ",
    "genetic effect, or an opt-in second effect such as `permanent(1 | id)`, ",
    "`common_env(1 | group)`, or `maternal_genetic(1 | dam)`.",
    call. = FALSE
  )
}
