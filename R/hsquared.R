#' Fit a quantitative-genetic model
#'
#' The default path fits a univariate Gaussian animal model by REML:
#' `y ~ sex + age + animal(1 | id, pedigree = ped)`. Fitting needs a local
#' Julia and HSquared.jl. Preview the same call with
#' `control = hs_control(engine = "validate")`.
#'
#' What you may report is listed on
#' [Can I fit and report this?](https://itchyshin.github.io/hsquared/articles/current-limits.html)
#' -- not in [validation_status()]. Other routes live in [formula_status()]
#' and on that limits page.
#'
#' @param formula A model formula. The default path is
#'   `y ~ fixed + animal(1 | id, pedigree = ped)`. `animal(1 | id)` is also
#'   accepted when `data` is an [hs_data()] object with a pedigree.
#' @param data A data frame of model variables, or an [hs_data()] object
#'   whose `phenotypes` component holds them. Formula arguments such as
#'   `pedigree = pedigree` can name components of the bundle.
#' @param family A response family. The default path accepts `gaussian()`.
#'   Opt-in experimental paths accept `poisson()` and `binomial()`; see
#'   [formula_status()].
#' @param REML Logical; REML estimation. The default path supports REML only
#'   (`TRUE`). `REML = FALSE` (ML) is rejected.
#' @param control An object created by [hs_control()].
#' @param ... Reserved for future arguments.
#'
#' @return A `"hsquared_fit"` object from a successful default-path fit.
#'   Missing Julia produces an install-guidance error. When
#'   `engine = "validate"`, the validated model spec is returned invisibly
#'   as a named list (not the classed object from [model_spec()]).
#' @export
hsquared <- function(
  formula,
  data,
  family = stats::gaussian(),
  REML = TRUE,
  control = hs_control(),
  ...
) {
  if (missing(formula)) {
    stop("`formula` is required.", call. = FALSE)
  }
  if (missing(data)) {
    stop("`data` is required.", call. = FALSE)
  }
  data_name <- hs_deparse_user_expr(substitute(data))
  if (!inherits(control, "hs_control")) {
    stop("`control` must be created by `hs_control()`.", call. = FALSE)
  }
  dots <- list(...)
  force(dots)

  # The opt-in non-Gaussian target widens the accepted families before the
  # model spec is validated; poisson(log)/binomial(logit) fit on the latent
  # scale, every other path stays Gaussian-only.
  julia_target <- if (identical(control$engine, "julia")) {
    hs_engine_control_value(control, "target", "fit_animal_model")
  } else {
    NA_character_
  }
  if (identical(julia_target, "nongaussian")) {
    hs_validate_nongaussian_three_field_v09_dots(dots)
  }
  allow_families <- if (identical(julia_target, "nongaussian")) {
    c("gaussian", "poisson", "binomial")
  } else {
    "gaussian"
  }

  spec <- hs_build_model_spec(
    formula = formula,
    data = data,
    family = family,
    REML = REML,
    allow_families = allow_families
  )
  payload <- hs_build_bridge_payload(spec)
  hs_warn_unmodelled_repeated_records(spec)

  if (identical(control$engine, "fit")) {
    # MV-4 (doc 38): a multivariate `cbind(...)` Gaussian response auto-routes to
    # the multivariate REML fitter (dispatched below), no longer requiring the
    # opt-in `target = "multivariate"`. The spec fence (`R/model-spec.R`) already
    # guarantees an animal-only, single-effect, random-intercept multivariate
    # response, so clauses (4)-(7) of the frozen dispatch key hold here.
    if (identical(spec$random$animal$design, "random_regression")) {
      hs_abort_opt_in_next_call(
        "random_regression",
        formula,
        data_name
      )
    }
    if (!isTRUE(REML)) {
      hs_abort_unsupported_syntax(
        "The v0.1 default fit path estimates variance components by REML ",
        "(average-information REML). ML estimation (`REML = FALSE`) is not yet ",
        "implemented on the fit path; use `REML = TRUE` (the default).",
        call. = FALSE
      )
    }
    opt_in_effect <- setdiff(names(spec$random), "animal")
    if (length(opt_in_effect) > 0L) {
      effect_type <- opt_in_effect[[1L]]
      # Boole: maternal_genetic has two opt-in targets. Name the covered
      # direct_maternal sibling first. Do not auto-route. common_env /
      # permanent paste stays on the sibling helper.
      if (identical(effect_type, "maternal_genetic")) {
        hs_abort_maternal_genetic_default_path()
      }
      if (identical(effect_type, "single_step")) {
        hs_abort_single_step_default_path()
      }
      if (identical(effect_type, "genomic")) {
        hs_abort_genomic_default_path(formula, data_name)
      }
      # Pat UX: the natural second-effect formulas print a pasteable next
      # call. Other opt-in primaries keep the existing knob-name abort.
      if (effect_type %in% c("common_env", "permanent")) {
        hs_abort_opt_in_next_call(effect_type, formula, data_name)
      }
      # Bare `(1 | group)` i.i.d. effects live under the `iid_effects` list slot;
      # name them honestly rather than printing the internal slot name.
      model_label <- if (identical(effect_type, "iid_effects")) {
        "bare `(1 | group)` i.i.d. random-effect (multi-effect)"
      } else {
        paste0("`", effect_type, "`")
      }
      hs_abort_unsupported_syntax(
        "The ",
        model_label,
        " model is experimental and opt-in; the default `engine = \"fit\"` ",
        "path fits the single-effect Gaussian animal model only. Use ",
        "`control = hs_control(engine = \"julia\", engine_control = list(",
        "target = \"",
        hs_second_effect_target(effect_type),
        "\"))`.",
        call. = FALSE
      )
    }
    project <- hs_engine_control_value(
      control,
      "julia_project",
      hs_default_julia_project()
    )
    if (!hs_julia_bridge_available(project)) {
      stop(
        "Fitting the v0.1 Gaussian animal model requires the HSquared.jl Julia ",
        "engine (Julia, the `JuliaCall` R package, and a from-source checkout ",
        "of `HSquared.jl`, which is not yet a registered Julia package). To ",
        "fit: (1) `git clone https://github.com/itchyshin/HSquared.jl`; ",
        "(2) point the bridge at that checkout, either by setting the ",
        "`HSQUARED_JULIA_PROJECT` environment variable to the clone path or by ",
        "passing `control = hs_control(engine = \"julia\", engine_control = ",
        "list(julia_project = \"/path/to/HSquared.jl\"))`. To validate the ",
        "model contract without fitting (no Julia needed), use ",
        "`control = hs_control(engine = \"validate\")`.",
        call. = FALSE
      )
    }
    if (isTRUE(spec$response$multivariate)) {
      hs_warn_cbind_experimental_once()
      return(hs_fit_julia_multivariate_payload(
        payload,
        project = project,
        initial = hs_engine_control_value(control, "initial", NULL),
        iterations = hs_engine_control_value(control, "iterations", 2000L),
        genetic_structure = hs_validate_genetic_structure_control(
          control,
          "multivariate"
        )
      ))
    }
    return(hs_fit_julia_ai_reml_payload(
      payload,
      project = project,
      initial = hs_engine_control_value(
        control,
        "initial",
        c(sigma_a2 = 1, sigma_e2 = 1)
      ),
      iterations = hs_engine_control_value(control, "iterations", 100L),
      em_warmup = hs_engine_control_value(control, "em_warmup", 0L)
    ))
  }

  if (identical(control$engine, "julia")) {
    target <- hs_validate_julia_target(hs_engine_control_value(
      control,
      "target",
      "fit_animal_model"
    ))
    # MV-4 (doc 38 section H2): a multivariate `cbind(...)` response under
    # `engine = "julia"` with no explicit target auto-selects the multivariate
    # target, mirroring the default-path auto-route. An explicit
    # non-multivariate target with a multivariate response still errors below.
    if (
      isTRUE(spec$response$multivariate) &&
        identical(target, "fit_animal_model")
    ) {
      target <- "multivariate"
    }
    genetic_structure <- hs_validate_genetic_structure_control(control, target)
    # hsquared#212: error on an `engine_control` key this target does not
    # honour (named, alongside the target) rather than silently discarding
    # it, before any target-specific dispatch below runs.
    hs_engine_control_forwarding(control, target)
    if (
      isTRUE(spec$response$multivariate) && !identical(target, "multivariate")
    ) {
      hs_abort_unsupported_syntax(
        "A `cbind(...)` multivariate response requires the multivariate Julia ",
        "engine path; the `",
        target,
        "` target fits a univariate response.",
        call. = FALSE
      )
    }
    if (
      identical(target, "multivariate") && !isTRUE(spec$response$multivariate)
    ) {
      hs_abort_unsupported_syntax(
        "`target = \"multivariate\"` requires a `cbind(trait1, trait2, ...)` ",
        "response with `animal(1 | id, pedigree = ped)`.",
        call. = FALSE
      )
    }
    rr_design <- identical(spec$random$animal$design, "random_regression")
    if (rr_design && !identical(target, "random_regression")) {
      hs_abort_unsupported_syntax(
        "An `animal(rr(...) | id)` random-regression term requires the opt-in ",
        "`target = \"random_regression\"` Julia engine path. The `",
        target,
        "` target fits a random-intercept response.",
        call. = FALSE
      )
    }
    if (identical(target, "random_regression") && !rr_design) {
      hs_abort_unsupported_syntax(
        "`target = \"random_regression\"` requires an ",
        "`animal(rr(covariate, order = k) | id, pedigree = ped)` term in the ",
        "formula.",
        call. = FALSE
      )
    }
    # ML estimation is not implemented in v0.1. The estimation targets either
    # run the ML optimizer (`fit_animal_model`) or are REML-only
    # (`sparse_reml`/`ai_reml`, which would otherwise silently ignore the ML
    # request). Reject `REML = FALSE` for all of them. `henderson_mme` solves at
    # supplied variances, so its method label is cosmetic and it is exempt;
    # `snp_blup` is exempt ONLY when variances are supplied -- without them it
    # genuinely REML-estimates `sigma_g2`/`sigma_e2` (`fit_snp_blup_reml`), so
    # `REML = FALSE` must not be silently accepted there.
    snp_blup_supplied <- identical(target, "snp_blup") &&
      !is.null(hs_engine_control_value(control, "variance_components", NULL))
    supplied_variance_exempt <- identical(target, "henderson_mme") ||
      identical(target, "metafounder") ||
      snp_blup_supplied
    if (!isTRUE(REML) && !supplied_variance_exempt) {
      hs_abort_unsupported_syntax(
        "ML estimation (`REML = FALSE`) is not implemented; the v0.1 fit path ",
        "estimates variance components by REML. Use `REML = TRUE`.",
        call. = FALSE
      )
    }
    # A non-default random effect only fits through one of its opt-in targets;
    # the single-effect estimators would silently ignore it.
    second_effect <- setdiff(names(spec$random), "animal")
    if (length(second_effect) > 0L) {
      allowed <- hs_effect_targets(second_effect[[1L]])
      if (!target %in% allowed) {
        # Bare `(1 | group)` i.i.d. effects are held under the `iid_effects` list
        # slot rather than a named `<type>()` call; name them honestly.
        effect_label <- if (identical(second_effect[[1L]], "iid_effects")) {
          "bare `(1 | group)` i.i.d. random"
        } else {
          paste0("`", second_effect[[1L]], "(...)`")
        }
        hs_abort_unsupported_syntax(
          "The formula has a ",
          effect_label,
          " term, so it needs `target = \"",
          allowed[[1L]],
          "\"`",
          if (length(allowed) > 1L) {
            paste0(" (or \"", allowed[[2L]], "\")")
          } else {
            ""
          },
          ". The `",
          target,
          "` target fits the single additive-genetic effect only.",
          if (identical(second_effect[[1L]], "maternal_genetic")) {
            hs_maternal_genetic_allowed_targets_note()
          } else if (identical(second_effect[[1L]], "single_step")) {
            hs_single_step_allowed_targets_note()
          } else {
            ""
          },
          "\n\nClosest working call:\n\n",
          hs_format_next_call(
            formula,
            data_name,
            if (identical(second_effect[[1L]], "maternal_genetic")) {
              "direct_maternal"
            } else {
              allowed[[1L]]
            }
          ),
          call. = FALSE
        )
      }
    }
    if (identical(target, "multivariate")) {
      return(hs_fit_julia_multivariate_payload(
        payload,
        project = hs_engine_control_value(
          control,
          "julia_project",
          hs_default_julia_project()
        ),
        initial = hs_engine_control_value(control, "initial", NULL),
        iterations = hs_engine_control_value(
          control,
          "iterations",
          2000L
        ),
        genetic_structure = genetic_structure
      ))
    }
    if (identical(target, "random_regression")) {
      return(hs_fit_julia_random_regression_payload(
        payload,
        project = hs_engine_control_value(
          control,
          "julia_project",
          hs_default_julia_project()
        ),
        iterations = hs_engine_control_value(
          control,
          "iterations",
          2000L
        )
      ))
    }
    if (identical(target, "nongaussian")) {
      if (identical(family$family, "gaussian")) {
        hs_abort_unsupported_syntax(
          "`target = \"nongaussian\"` fits non-Gaussian families ",
          "(`poisson(log)`, `binomial(logit)`); `family = gaussian()` fits ",
          "through the default path or `target = \"ai_reml\"`.",
          call. = FALSE
        )
      }
      return(hs_fit_julia_nongaussian_payload(
        payload,
        project = hs_engine_control_value(
          control,
          "julia_project",
          hs_default_julia_project()
        ),
        family = family,
        marginal = hs_engine_control_value(control, "marginal", "laplace"),
        iterations = hs_engine_control_value(control, "iterations", 200L),
        initial = hs_engine_control_value(control, "initial", NULL),
        restart_check = hs_engine_control_value(control, "restart_check", FALSE)
      ))
    }
    if (identical(target, "henderson_mme")) {
      return(hs_fit_julia_henderson_mme_payload(
        payload,
        project = hs_engine_control_value(
          control,
          "julia_project",
          hs_default_julia_project()
        ),
        variance_components = hs_engine_control_value(
          control,
          "variance_components",
          NULL
        )
      ))
    }

    if (identical(target, "metafounder")) {
      if (is.null(spec$random$metafounder)) {
        hs_abort_unsupported_syntax(
          "`target = \"metafounder\"` requires ",
          "`metafounder(1 | id, pedigree = ped, group = mf_group, ",
          "Gamma = Gamma)` in the formula.",
          call. = FALSE
        )
      }
      return(hs_fit_julia_metafounder_payload(
        payload,
        project = hs_engine_control_value(
          control,
          "julia_project",
          hs_default_julia_project()
        ),
        variance_components = hs_engine_control_value(
          control,
          "variance_components",
          NULL
        )
      ))
    }

    if (identical(target, "sparse_reml")) {
      return(hs_fit_julia_sparse_reml_payload(
        payload,
        project = hs_engine_control_value(
          control,
          "julia_project",
          hs_default_julia_project()
        ),
        initial = hs_engine_control_value(
          control,
          "initial",
          c(sigma_a2 = 1, sigma_e2 = 1)
        ),
        iterations = hs_engine_control_value(
          control,
          "iterations",
          1000L
        )
      ))
    }

    if (identical(target, "ai_reml")) {
      return(hs_fit_julia_ai_reml_payload(
        payload,
        project = hs_engine_control_value(
          control,
          "julia_project",
          hs_default_julia_project()
        ),
        initial = hs_engine_control_value(
          control,
          "initial",
          c(sigma_a2 = 1, sigma_e2 = 1)
        ),
        iterations = hs_engine_control_value(
          control,
          "iterations",
          100L
        ),
        em_warmup = hs_engine_control_value(control, "em_warmup", 0L)
      ))
    }

    if (identical(target, "repeatability")) {
      if (is.null(spec$random$permanent)) {
        hs_abort_unsupported_syntax(
          "`target = \"repeatability\"` requires a `permanent(1 | id)` term ",
          "alongside `animal(1 | id, ...)` in the formula.",
          call. = FALSE
        )
      }
      # `initial`/`iterations` reach the dense estimator but NOT the sparse one
      # (HSquared.jl#343). Say so rather than accepting them and quietly
      # ignoring them -- a start value the user tuned to rescue convergence,
      # silently dropped, is the worst of both behaviours.
      if (
        identical(
          hs_engine_control_value(control, "scale_method", "dense"),
          "auto"
        )
      ) {
        inert <- c("initial", "iterations")[
          c(
            !is.null(hs_engine_control_value(control, "initial", NULL)),
            !is.null(hs_engine_control_value(control, "iterations", NULL))
          )
        ]
        if (length(inert) > 0L) {
          warning(
            "`target = \"repeatability\"` with `scale_method = \"auto\"` ",
            "ignores ",
            paste0("`", inert, "`", collapse = " and "),
            ": the sparse route picks its own start and iteration count ",
            "(HSquared.jl#343). Use `scale_method = \"dense\"` to control them.",
            call. = FALSE
          )
        }
      }
      return(hs_fit_julia_repeatability_payload(
        payload,
        project = hs_engine_control_value(
          control,
          "julia_project",
          hs_default_julia_project()
        ),
        initial = hs_engine_control_value(
          control,
          "initial",
          c(sigma_a2 = 1, sigma_pe2 = 1, sigma_e2 = 1)
        ),
        iterations = hs_engine_control_value(
          control,
          "iterations",
          200L
        ),
        max_dense_cells = hs_engine_control_value(
          control,
          "max_dense_cells",
          1e6
        ),
        # "dense" (default) keeps the covered validation-scale estimator.
        # "auto" routes the SAME model through the engine's sparse K-effect
        # AI-REML, which is the only way a real repeated-measures pedigree
        # clears the dense nobs^2 + nanimals^2 ceiling.
        scale_method = hs_engine_control_value(
          control,
          "scale_method",
          "dense"
        )
      ))
    }

    if (identical(target, "two_effect")) {
      if (
        is.null(spec$random$common_env) &&
          is.null(spec$random$maternal_genetic)
      ) {
        hs_abort_unsupported_syntax(
          "`target = \"two_effect\"` requires a `common_env(1 | group)` term, ",
          "or a `maternal_genetic(1 | dam)` term, alongside ",
          "`animal(1 | id, ...)` in the formula.",
          call. = FALSE
        )
      }
      return(hs_fit_julia_two_effect_payload(
        payload,
        project = hs_engine_control_value(
          control,
          "julia_project",
          hs_default_julia_project()
        ),
        initial = hs_engine_control_value(
          control,
          "initial",
          c(sigma_a2 = 1, sigma_c2 = 1, sigma_e2 = 1)
        ),
        iterations = hs_engine_control_value(
          control,
          "iterations",
          200L
        )
      ))
    }

    if (identical(target, "multi_effect")) {
      # A `permanent(1 | id)` formula is deliberately NOT accepted here: the
      # term-to-target router above already sends it to `target =
      # "repeatability"` with the closest working call spelled out, and that
      # target keeps the PE-specific labels, `repeatability()` and
      # `permanent_effects()`. Accepting it here too would give the same model
      # a second, worse-labelled surface.
      if (
        is.null(spec$random$iid_effects) ||
          length(spec$random$iid_effects) == 0L
      ) {
        hs_abort_unsupported_syntax(
          "`target = \"multi_effect\"` requires an `animal(1 | id, ...)` term ",
          "plus one or more bare `(1 | group)` i.i.d. random effects in the ",
          "formula.",
          call. = FALSE
        )
      }
      return(hs_fit_julia_n_effect_payload(
        payload,
        project = hs_engine_control_value(
          control,
          "julia_project",
          hs_default_julia_project()
        ),
        # `scale_method = "auto"` (opt-in) lets the engine route large multi-effect
        # problems through the EXPERIMENTAL matrix-free Monte-Carlo fit when the dense
        # factorization is infeasible; the default "dense" is the covered
        # validation-scale path. See ?hsquared (engine_control) and the engine's
        # `docs/design/25-completion-ultraplan.md` (V8.6).
        scale_method = hs_engine_control_value(
          control,
          "scale_method",
          "dense"
        ),
        # hsquared#212: forward `initial`/`iterations`; NULL (unsupplied)
        # reproduces the pre-#212-fix default exactly. NOTE: honoured on the
        # `scale_method = "dense"` route only -- `"auto"` does not forward
        # them yet (HSquared.jl#343).
        initial = hs_engine_control_value(control, "initial", NULL),
        iterations = hs_engine_control_value(control, "iterations", NULL)
      ))
    }

    if (identical(target, "direct_maternal")) {
      if (is.null(spec$random$maternal_genetic)) {
        hs_abort_unsupported_syntax(
          "`target = \"direct_maternal\"` requires a `maternal_genetic(1 | dam)` ",
          "term alongside `animal(1 | id, pedigree = ped)` in the formula. ",
          "This is the opt-in correlated direct-maternal model (Phase 4). ",
          "The INDEPENDENT two-effect maternal model uses ",
          "`target = \"two_effect\"` instead.",
          call. = FALSE
        )
      }
      return(hs_fit_julia_direct_maternal_payload(
        payload,
        project = hs_engine_control_value(
          control,
          "julia_project",
          hs_default_julia_project()
        ),
        # hsquared#212: forward `initial`/`iterations`; NULL (unsupplied)
        # reproduces the pre-#212-fix default exactly.
        initial = hs_engine_control_value(control, "initial", NULL),
        iterations = hs_engine_control_value(control, "iterations", NULL)
      ))
    }

    if (identical(target, "snp_blup")) {
      genomic_effect <- spec$random$genomic
      if (
        is.null(genomic_effect) ||
          !identical(genomic_effect$source, "markers")
      ) {
        hs_abort_unsupported_syntax(
          "`target = \"snp_blup\"` requires a `genomic(1 | id, markers = M)` ",
          "term (a raw marker matrix). SNP-BLUP estimates marker effects, so ",
          "a precomputed `Ginv` cannot be used.",
          call. = FALSE
        )
      }
      vc_raw <- hs_engine_control_value(control, "variance_components", NULL)
      project <- hs_engine_control_value(
        control,
        "julia_project",
        hs_default_julia_project()
      )
      # No supplied variances -> estimate them by REML from the markers; supplied
      # variances -> the supplied-variance solve (byte-identical to before).
      if (is.null(vc_raw)) {
        return(hs_fit_julia_snp_blup_reml_payload(payload, project = project))
      }
      vc <- hs_validate_snp_blup_variances(vc_raw)
      return(hs_fit_julia_snp_blup_payload(
        payload,
        project = project,
        variance_components = vc
      ))
    }

    if (identical(target, "single_step_construct")) {
      ss <- spec$random$single_step
      if (!is.null(ss) && identical(ss$source, "metafounder_construct")) {
        hs_abort_unsupported_syntax(
          "`single_step()` with `group` and `Gamma` is the planned ",
          "`target = \"metafounder_single_step\"` contract, not the ordinary ",
          "`target = \"single_step_construct\"` path.",
          call. = FALSE
        )
      }
      if (is.null(ss) || !identical(ss$source, "construct")) {
        hs_abort_unsupported_syntax(
          "`target = \"single_step_construct\"` requires a ",
          "`single_step(1 | id, pedigree = ped, markers = M)` term in the ",
          "formula.",
          call. = FALSE
        )
      }
      return(hs_fit_julia_single_step_construct_payload(
        payload,
        project = hs_engine_control_value(
          control,
          "julia_project",
          hs_default_julia_project()
        ),
        initial = hs_engine_control_value(
          control,
          "initial",
          c(sigma_a2 = 1, sigma_e2 = 1)
        ),
        iterations = hs_engine_control_value(
          control,
          "iterations",
          100L
        )
      ))
    }

    if (identical(target, "metafounder_single_step")) {
      ss <- spec$random$single_step
      if (is.null(ss) || !identical(ss$source, "metafounder_construct")) {
        hs_abort_unsupported_syntax(
          "`target = \"metafounder_single_step\"` requires ",
          "`single_step(1 | id, pedigree = ped, markers = M, group = ",
          "mf_group, Gamma = Gamma)` in the formula.",
          call. = FALSE
        )
      }
      return(hs_fit_julia_metafounder_single_step_payload(
        payload,
        project = hs_engine_control_value(
          control,
          "julia_project",
          hs_default_julia_project()
        ),
        initial = hs_engine_control_value(
          control,
          "initial",
          c(sigma_a2 = 1, sigma_e2 = 1)
        ),
        iterations = hs_engine_control_value(
          control,
          "iterations",
          100L
        )
      ))
    }

    # Experimental supplied-relationship primary effects. They carry a supplied
    # relationship inverse (relmat's Kinv = solve(K), or precision's Q) in the
    # same spec/payload shape as genomic's Ginv, so they fit through the SAME
    # covered supplied-relationship-inverse REML path (fit_ai_reml). Opt-in,
    # experimental, NOT the default; does not change the covered surface.
    if (target %in% c("relmat", "precision")) {
      if (is.null(spec$random[[target]])) {
        arg <- if (identical(target, "relmat")) "K" else "Q"
        stop(
          "`target = \"",
          target,
          "\"` requires a `",
          target,
          "(1 | id, ",
          arg,
          " = ",
          arg,
          ")` term in the formula.",
          call. = FALSE
        )
      }
      return(hs_fit_julia_genomic_payload(
        payload,
        project = hs_engine_control_value(
          control,
          "julia_project",
          hs_default_julia_project()
        ),
        initial = hs_engine_control_value(
          control,
          "initial",
          c(sigma_a2 = 1, sigma_e2 = 1)
        ),
        iterations = hs_engine_control_value(
          control,
          "iterations",
          100L
        )
      ))
    }

    if (target %in% c("genomic", "single_step")) {
      if (is.null(spec$random[[target]])) {
        arg <- if (identical(target, "genomic")) "Ginv" else "Hinv"
        hs_abort_unsupported_syntax(
          "`target = \"",
          target,
          "\"` requires a `",
          target,
          "(1 | id, ",
          arg,
          " = ",
          arg,
          ")` term in the formula.",
          call. = FALSE
        )
      }
      return(hs_fit_julia_genomic_payload(
        payload,
        project = hs_engine_control_value(
          control,
          "julia_project",
          hs_default_julia_project()
        ),
        initial = hs_engine_control_value(
          control,
          "initial",
          c(sigma_a2 = 1, sigma_e2 = 1)
        ),
        iterations = hs_engine_control_value(
          control,
          "iterations",
          100L
        )
      ))
    }

    return(hs_fit_julia_payload(
      payload,
      project = hs_engine_control_value(
        control,
        "julia_project",
        hs_default_julia_project()
      ),
      initial = hs_engine_control_value(
        control,
        "initial",
        c(sigma_a2 = 1, sigma_e2 = 1)
      ),
      max_dense_cells = hs_engine_control_value(
        control,
        "max_dense_cells",
        1e6
      )
    ))
  }

  # `engine == "validate"`: reaching here means `hs_build_model_spec()` and
  # `hs_build_bridge_payload()` both succeeded, so the v0.1 animal-model
  # contract is validated. Confirm with a message and return the validated spec
  # invisibly so it can be assigned and inspected, rather than stopping. Use the
  # default `control` to fit (requires the HSquared.jl Julia engine), or
  # `hs_control(engine = "julia")` for advanced engine control.
  message(
    "Validated the v0.1 animal-model contract; no model was fitted ",
    "(`engine = \"validate\"`). Julia fit target: `HSquared.",
    spec$bridge$target,
    "`. The validated spec is returned invisibly as a named list; assign it ",
    "to inspect the parsed contract (e.g. `spec$bridge$target`)."
  )
  invisible(spec)
}
