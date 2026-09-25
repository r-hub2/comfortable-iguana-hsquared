hs_julia_bridge_state <- new.env(parent = emptyenv())

hs_default_julia_project <- function() {
  env <- Sys.getenv("HSQUARED_JULIA_PROJECT", unset = "")
  if (nzchar(env)) {
    return(env)
  }

  normalizePath(
    file.path(dirname(system.file(package = "hsquared")), "HSquared.jl"),
    winslash = "/",
    mustWork = FALSE
  )
}

hs_julia_bridge_available <- function(project = hs_default_julia_project()) {
  requireNamespace("JuliaCall", quietly = TRUE) &&
    nzchar(Sys.which("julia")) &&
    file.exists(file.path(project, "Project.toml"))
}

# Post-fit engine calls the bridge guards with a Julia `try` (SEs, intervals,
# plot data) record their failure in `hsq_bridge_errors` instead of swallowing
# it (hsquared / HSquared.jl#351). A route with a recording catch opens its
# fit command with `hs_julia_bridge_errors_reset` and returns through
# `hs_julia_surface_bridge_errors()`, which attaches the record as
# `attr(fit, "bridge_errors")` and warns once (see R/conditions.R).
hs_julia_bridge_errors_reset <- "hsq_bridge_errors = Dict{String,String}();"

# The recording `catch` of one guarded post-fit call; `name` is the engine
# function that threw, as the warning reports it:
#   catch err; hsq_bridge_errors["<name>"] = sprint(showerror, err); nothing; end;
hs_julia_catch_record <- function(name) {
  paste0(
    "catch err; hsq_bridge_errors[\"",
    name,
    "\"] = sprint(showerror, err); nothing; end;"
  )
}

# One guarded single-line call: `<slot> = try; <call>; <recording catch>`.
hs_julia_try_slot <- function(slot, call, name) {
  paste0(slot, " = try; ", call, "; ", hs_julia_catch_record(name))
}

hs_julia_read_bridge_errors <- function() {
  keys <- as.character(unlist(JuliaCall::julia_eval(
    "sort!(collect(keys(hsq_bridge_errors)))"
  )))
  if (length(keys) == 0L) {
    return(character())
  }
  values <- as.character(unlist(JuliaCall::julia_eval(
    "[hsq_bridge_errors[k] for k in sort!(collect(keys(hsq_bridge_errors)))]"
  )))
  stats::setNames(values, keys)
}

hs_julia_surface_bridge_errors <- function(fit) {
  hs_attach_bridge_errors(fit, hs_julia_read_bridge_errors())
}

hs_julia_attach_standard_plot_data <- function() {
  JuliaCall::julia_command(paste(
    "if isdefined(HSquared, :variance_components_plot_data);",
    hs_julia_try_slot(
      "hsq_vcpd",
      "HSquared.variance_components_plot_data(hsq_fit)",
      "variance_components_plot_data"
    ),
    "if hsq_vcpd !== nothing;",
    "hsq_result = merge(hsq_result, (",
    "variance_components_plot_data = hsq_vcpd,));",
    "end;",
    "end;",
    "if isdefined(HSquared, :breeding_values_plot_data);",
    hs_julia_try_slot(
      "hsq_bvpd",
      "HSquared.breeding_values_plot_data(hsq_fit)",
      "breeding_values_plot_data"
    ),
    "if hsq_bvpd !== nothing;",
    "hsq_result = merge(hsq_result, (breeding_values_plot_data = hsq_bvpd,));",
    "end;",
    "end;"
  ))
  invisible(TRUE)
}

hs_fit_julia_payload <- function(
  payload,
  project = hs_default_julia_project(),
  initial = c(sigma_a2 = 1, sigma_e2 = 1),
  max_dense_cells = 1e6
) {
  if (!inherits(payload, "hs_bridge_payload")) {
    stop("`payload` must be an internal `hs_bridge_payload`.", call. = FALSE)
  }
  if (!hs_julia_bridge_available(project)) {
    stop(
      "The experimental Julia bridge requires Julia, the `JuliaCall` R ",
      "package, and a local `HSquared.jl` project.",
      call. = FALSE
    )
  }

  hs_julia_setup(project)
  initial <- hs_validate_initial_variances(initial)
  hs_julia_assign_payload(payload, initial)
  JuliaCall::julia_assign("hsq_mdc", hs_validate_max_dense_cells(max_dense_cells))
  hs_julia_fit(
    JuliaCall::julia_command(paste(
      hs_julia_bridge_errors_reset,
      "hsq_ped = HSquared.normalize_pedigree(hsq_id, hsq_sire, hsq_dam);",
      "hsq_Ainv = HSquared.pedigree_inverse(hsq_ped);",
      "hsq_fit = HSquared.fit_animal_model(",
      "hsq_y, hsq_X, hsq_Z, hsq_Ainv;",
      "ids = hsq_ped.ids,",
      "method = Symbol(hsq_method),",
      "initial = (sigma_a2 = hsq_initial_sigma_a2,",
      "sigma_e2 = hsq_initial_sigma_e2),",
      "max_dense_cells = Int(hsq_mdc)",
      ");",
      "hsq_result = HSquared.result_payload(hsq_fit);",
      # Enrich with PEV/reliability only for older engines whose result_payload
      # does not already carry them; current engines emit them via :selinv, and
      # re-merging would clobber that standard field with a redundant :auto solve.
      "if !hasproperty(hsq_result, :prediction_error_variance) &&",
      "isdefined(HSquared, :prediction_error_variance) &&",
      "isdefined(HSquared, :reliability);",
      "hsq_result = merge(hsq_result, (",
      "prediction_error_variance =",
      "HSquared.prediction_error_variance(hsq_fit),",
      "reliability = HSquared.reliability(hsq_fit)));",
      "end;"
    )),
    hint = hs_dense_route_hint
  )
  hs_julia_attach_standard_plot_data()

  raw <- JuliaCall::julia_eval(
    "Dict(String(k) => getfield(hsq_result, k) for k in keys(hsq_result))"
  )
  result <- hs_normalize_julia_result(raw, payload)
  fit <- hs_new_fit(
    spec = list(
      method = payload$method,
      family = list(family = payload$family, link = "identity")
    ),
    payload = payload,
    result = result,
    engine = "HSquared.jl"
  )
  hs_julia_surface_bridge_errors(fit)
}

hs_fit_julia_henderson_mme_payload <- function(
  payload,
  project = hs_default_julia_project(),
  variance_components
) {
  if (!inherits(payload, "hs_bridge_payload")) {
    stop("`payload` must be an internal `hs_bridge_payload`.", call. = FALSE)
  }
  if (!hs_julia_bridge_available(project)) {
    stop(
      "The experimental Julia bridge requires Julia, the `JuliaCall` R ",
      "package, and a local `HSquared.jl` project.",
      call. = FALSE
    )
  }

  variance_components <- hs_validate_supplied_variances(variance_components)
  hs_julia_setup(project)
  hs_julia_assign_payload(payload, variance_components)
  JuliaCall::julia_assign(
    "hsq_supplied_sigma_a2",
    unname(variance_components[["sigma_a2"]])
  )
  JuliaCall::julia_assign(
    "hsq_supplied_sigma_e2",
    unname(variance_components[["sigma_e2"]])
  )
  hs_julia_fit(
    JuliaCall::julia_command(paste(
      "hsq_ped = HSquared.normalize_pedigree(hsq_id, hsq_sire, hsq_dam);",
      "hsq_Ainv = HSquared.pedigree_inverse(hsq_ped);",
      "hsq_spec = HSquared.animal_model_spec(",
      "hsq_y, hsq_X, hsq_Z, hsq_Ainv;",
      "ids = hsq_ped.ids, method = Symbol(hsq_method));",
      "hsq_mme = HSquared.henderson_mme(",
      "hsq_spec, hsq_supplied_sigma_a2, hsq_supplied_sigma_e2);",
      "hsq_mme_bv = HSquared.breeding_values(hsq_mme);",
      "hsq_mme_raw = Dict(",
      "\"fixed_effects\" => HSquared.fixed_effects(hsq_mme),",
      "\"animal_ids\" => hsq_mme_bv.ids,",
      "\"animal_effects\" => hsq_mme_bv.values,",
      "\"fitted\" => HSquared.fitted_values(hsq_mme),",
      # PEV/reliability are now standard on the Henderson MME result
      # (validation-scale): prediction_error_variance/reliability default to
      # method = :auto (selected inverse by storage; :dense remains the
      # explicit validation oracle), so they are attached unconditionally
      # rather than probed.
      "\"prediction_error_variance\" =>",
      "HSquared.prediction_error_variance(hsq_mme),",
      "\"reliability\" => HSquared.reliability(hsq_mme),",
      "\"nobs\" => length(hsq_y)",
      ");"
    )),
    hint = hs_dense_scale_hint
  )

  raw <- JuliaCall::julia_eval("hsq_mme_raw")
  result <- hs_normalize_julia_henderson_mme_result(
    raw,
    payload,
    variance_components
  )
  hs_new_fit(
    spec = list(
      method = payload$method,
      family = list(family = payload$family, link = "identity"),
      target = "henderson_mme"
    ),
    payload = payload,
    result = result,
    engine = "HSquared.jl"
  )
}

hs_fit_julia_metafounder_payload <- function(
  payload,
  project = hs_default_julia_project(),
  variance_components
) {
  if (!inherits(payload, "hs_bridge_payload")) {
    stop("`payload` must be an internal `hs_bridge_payload`.", call. = FALSE)
  }
  if (
    !identical(payload$relationship_source, "metafounder") ||
      is.null(payload$group_of) ||
      is.null(payload$Gamma)
  ) {
    stop(
      "Internal bridge error: the metafounder payload is incomplete ",
      "(needs group_of and Gamma).",
      call. = FALSE
    )
  }
  variance_components <- hs_validate_supplied_variances(
    variance_components,
    target = "metafounder"
  )
  if (!hs_julia_bridge_available(project)) {
    stop(
      "The experimental Julia bridge requires Julia, the `JuliaCall` R ",
      "package, and a local `HSquared.jl` project.",
      call. = FALSE
    )
  }

  hs_julia_setup(project)
  hs_julia_assign_payload(payload, variance_components)
  JuliaCall::julia_assign("hsq_group_of", unname(payload$group_of))
  JuliaCall::julia_assign("hsq_Gamma_vec", as.numeric(payload$Gamma))
  JuliaCall::julia_assign("hsq_Gamma_n", as.integer(nrow(payload$Gamma)))
  JuliaCall::julia_assign(
    "hsq_supplied_sigma_a2",
    unname(variance_components[["sigma_a2"]])
  )
  JuliaCall::julia_assign(
    "hsq_supplied_sigma_e2",
    unname(variance_components[["sigma_e2"]])
  )
  hs_julia_fit(
    JuliaCall::julia_command(paste(
      "hsq_ped = HSquared.normalize_pedigree(hsq_id, hsq_sire, hsq_dam);",
      "collect(String, hsq_ped.ids) == hsq_id ||",
      "error(\"metafounder: engine pedigree order != R order\");",
      "hsq_Gamma = reshape(collect(Float64, hsq_Gamma_vec),",
      "Int(hsq_Gamma_n), Int(hsq_Gamma_n));",
      "hsq_mme = HSquared.metafounder_animal_model(",
      "hsq_y, hsq_X, hsq_Z, hsq_ped, hsq_group_of, hsq_Gamma,",
      "hsq_supplied_sigma_a2, hsq_supplied_sigma_e2;",
      "ids = hsq_ped.ids);",
      "hsq_mme_bv = HSquared.breeding_values(hsq_mme);",
      "hsq_mme_raw = Dict(",
      "\"fixed_effects\" => HSquared.fixed_effects(hsq_mme),",
      "\"animal_ids\" => hsq_mme_bv.ids,",
      "\"animal_effects\" => hsq_mme_bv.values,",
      "\"fitted\" => HSquared.fitted_values(hsq_mme),",
      "\"prediction_error_variance\" =>",
      "HSquared.prediction_error_variance(hsq_mme),",
      "\"reliability\" => HSquared.reliability(hsq_mme),",
      "\"nobs\" => length(hsq_y)",
      ");"
    )),
    hint = hs_dense_scale_hint
  )

  raw <- JuliaCall::julia_eval("hsq_mme_raw")
  result <- hs_normalize_julia_henderson_mme_result(
    raw,
    payload,
    variance_components
  )
  result$variance_components$component[
    result$variance_components$component == "animal"
  ] <- "metafounder"
  result$heritability$term[result$heritability$term == "animal"] <-
    "metafounder"
  names(result$random_effects)[
    names(result$random_effects) == "animal"
  ] <- "metafounder"
  result$diagnostics$target <- "metafounder"
  result$diagnostics$variance_components <- "supplied_metafounder"
  result$diagnostics$gamma_source <- "supplied"
  hs_new_fit(
    spec = list(
      method = payload$method,
      family = list(family = payload$family, link = "identity"),
      target = "metafounder"
    ),
    payload = payload,
    result = result,
    engine = "HSquared.jl"
  )
}

hs_fit_julia_sparse_reml_payload <- function(
  payload,
  project = hs_default_julia_project(),
  initial = c(sigma_a2 = 1, sigma_e2 = 1),
  iterations = 1000L
) {
  if (!inherits(payload, "hs_bridge_payload")) {
    stop("`payload` must be an internal `hs_bridge_payload`.", call. = FALSE)
  }
  if (!hs_julia_bridge_available(project)) {
    stop(
      "The experimental Julia bridge requires Julia, the `JuliaCall` R ",
      "package, and a local `HSquared.jl` project.",
      call. = FALSE
    )
  }

  initial <- hs_validate_initial_variances(initial)
  iterations <- hs_validate_iterations(iterations)
  hs_julia_setup(project)
  hs_julia_assign_payload(payload, initial)
  JuliaCall::julia_assign("hsq_iterations", iterations)
  hs_julia_fit(
    JuliaCall::julia_command(paste(
      hs_julia_bridge_errors_reset,
      "hsq_ped = HSquared.normalize_pedigree(hsq_id, hsq_sire, hsq_dam);",
      "hsq_Ainv = HSquared.pedigree_inverse(hsq_ped);",
      "hsq_spec = HSquared.animal_model_spec(",
      "hsq_y, hsq_X, hsq_Z, hsq_Ainv;",
      "ids = hsq_ped.ids, method = :REML);",
      "hsq_fit = HSquared.fit_sparse_reml(",
      "hsq_spec;",
      "initial = (sigma_a2 = hsq_initial_sigma_a2,",
      "sigma_e2 = hsq_initial_sigma_e2),",
      "iterations = hsq_iterations);",
      "hsq_result = HSquared.result_payload(hsq_fit);",
      # Enrich with PEV/reliability only for older engines whose result_payload
      # does not already carry them; current engines emit them via :selinv, and
      # re-merging would clobber that standard field with a redundant :auto solve.
      "if !hasproperty(hsq_result, :prediction_error_variance) &&",
      "isdefined(HSquared, :prediction_error_variance) &&",
      "isdefined(HSquared, :reliability) &&",
      "applicable(HSquared.prediction_error_variance, hsq_fit) &&",
      "applicable(HSquared.reliability, hsq_fit);",
      "hsq_result = merge(hsq_result, (",
      "prediction_error_variance =",
      "HSquared.prediction_error_variance(hsq_fit),",
      "reliability = HSquared.reliability(hsq_fit)));",
      "end;"
    )),
    hint = hs_dense_scale_hint
  )
  hs_julia_attach_standard_plot_data()

  raw <- JuliaCall::julia_eval(
    "Dict(String(k) => getfield(hsq_result, k) for k in keys(hsq_result))"
  )
  result <- hs_normalize_julia_result(raw, payload)
  result$diagnostics$variance_components <- "estimated_sparse_reml"
  fit <- hs_new_fit(
    spec = list(
      # fit_sparse_reml is a REML-only optimizer; stamp what was computed
      # rather than echoing the requested method.
      method = "REML",
      family = list(family = payload$family, link = "identity"),
      target = "sparse_reml"
    ),
    payload = payload,
    result = result,
    engine = "HSquared.jl"
  )
  hs_julia_surface_bridge_errors(fit)
}

hs_fit_julia_ai_reml_payload <- function(
  payload,
  project = hs_default_julia_project(),
  initial = c(sigma_a2 = 1, sigma_e2 = 1),
  iterations = 100L,
  em_warmup = 0L
) {
  if (!inherits(payload, "hs_bridge_payload")) {
    stop("`payload` must be an internal `hs_bridge_payload`.", call. = FALSE)
  }
  if (!hs_julia_bridge_available(project)) {
    stop(
      "The experimental Julia bridge requires Julia, the `JuliaCall` R ",
      "package, and a local `HSquared.jl` project.",
      call. = FALSE
    )
  }

  initial <- hs_validate_initial_variances(initial)
  iterations <- hs_validate_iterations(iterations)
  em_warmup <- hs_validate_em_warmup(em_warmup)
  hs_julia_setup(project)
  hs_julia_assign_payload(payload, initial)
  JuliaCall::julia_assign("hsq_iterations", iterations)
  JuliaCall::julia_assign("hsq_em_warmup", em_warmup)
  hs_julia_fit(
    JuliaCall::julia_command(paste(
    hs_julia_bridge_errors_reset,
    "hsq_ped = HSquared.normalize_pedigree(hsq_id, hsq_sire, hsq_dam);",
    "hsq_Ainv = HSquared.pedigree_inverse(hsq_ped);",
    "hsq_spec = HSquared.animal_model_spec(",
    "hsq_y, hsq_X, hsq_Z, hsq_Ainv;",
    "ids = hsq_ped.ids, method = :REML);",
    "hsq_fit = HSquared.fit_ai_reml(",
    "hsq_spec;",
    "initial = (sigma_a2 = hsq_initial_sigma_a2,",
    "sigma_e2 = hsq_initial_sigma_e2),",
    # em_warmup (engine V1-AI-REML): opt-in EM-REML warm-start before the AI step;
    # default 0 = byte-identical to the pre-warm-start engine call.
    "iterations = hsq_iterations, em_warmup = hsq_em_warmup);",
    "hsq_result = HSquared.result_payload(hsq_fit);",
    # Enrich with PEV/reliability only for older engines whose result_payload
    # does not already carry them; current engines emit them via :selinv, and
    # re-merging would clobber that standard field with a redundant :auto solve.
    "if !hasproperty(hsq_result, :prediction_error_variance) &&",
    "isdefined(HSquared, :prediction_error_variance) &&",
    "isdefined(HSquared, :reliability) &&",
    "applicable(HSquared.prediction_error_variance, hsq_fit) &&",
    "applicable(HSquared.reliability, hsq_fit);",
    "hsq_result = merge(hsq_result, (",
    "prediction_error_variance =",
    "HSquared.prediction_error_variance(hsq_fit),",
    "reliability = HSquared.reliability(hsq_fit)));",
    "end;",
    # Experimental, opt-in heritability CI (engine row V1-HERIT-CI, partial).
    # Guarded by a try: the engine throws when h2 is on the (0, 1) boundary,
    # which must not abort the fit (the throw is recorded and warned, #351).
    "if isdefined(HSquared, :heritability_interval) &&",
    "applicable(HSquared.heritability_interval, hsq_fit);",
    hs_julia_try_slot(
      "hsq_hi",
      "HSquared.heritability_interval(hsq_fit)",
      "heritability_interval"
    ),
    "if hsq_hi !== nothing;",
    "hsq_result = merge(hsq_result, (heritability_interval = hsq_hi,));",
    "end;",
    "end;",
    # Experimental, opt-in variance-component and heritability standard errors
    # (engine row V1-HERIT-CI, partial). variance_component_covariance() can
    # throw on a singular/ill-conditioned AI matrix, so each call is wrapped in
    # a try so an SE failure never aborts the fit; the throw is recorded and
    # surfaced as a warning rather than a silent absence (#351).
    "if isdefined(HSquared, :variance_component_standard_errors) &&",
    "applicable(HSquared.variance_component_standard_errors, hsq_fit);",
    hs_julia_try_slot(
      "hsq_vcse",
      "HSquared.variance_component_standard_errors(hsq_fit)",
      "variance_component_standard_errors"
    ),
    "if hsq_vcse !== nothing;",
    "hsq_result = merge(hsq_result, (variance_component_se = hsq_vcse,));",
    "end;",
    "end;",
    "if isdefined(HSquared, :heritability_standard_error) &&",
    "applicable(HSquared.heritability_standard_error, hsq_fit);",
    hs_julia_try_slot(
      "hsq_h2se",
      "HSquared.heritability_standard_error(hsq_fit)",
      "heritability_standard_error"
    ),
    "if hsq_h2se !== nothing;",
    "hsq_result = merge(hsq_result, (heritability_se = hsq_h2se,));",
    "end;",
    "end;"
    )),
    hint = hs_dense_scale_hint
  )
  hs_julia_attach_standard_plot_data()

  raw <- JuliaCall::julia_eval(
    "Dict(String(k) => getfield(hsq_result, k) for k in keys(hsq_result))"
  )
  result <- hs_normalize_julia_result(raw, payload)
  result$diagnostics$variance_components <- "estimated_ai_reml"
  fit <- hs_new_fit(
    spec = list(
      # fit_ai_reml is a REML-only (average-information) optimizer; stamp what
      # was computed rather than echoing the requested method.
      method = "REML",
      family = list(family = payload$family, link = "identity"),
      target = "ai_reml"
    ),
    payload = payload,
    result = result,
    engine = "HSquared.jl"
  )
  hs_julia_surface_bridge_errors(fit)
}

# Map an R `family` object to the engine's non-Gaussian family symbol:
# `poisson(log)` -> "poisson"; `binomial(logit)` -> "bernoulli" for a binary 0/1
# response (or an all-ones `cbind` total), or "binomial" when any per-record trial
# count in `n_trials` exceeds 1 (a `cbind(successes, failures)` counts response;
# `n_trials` may be a per-record vector). Other families are planned. The rule is
# vector-safe: `any(n_trials > 1L)`, so a vector whose first element is 1 (e.g.
# c(1, 4, 5)) is still classified Binomial, not silently reduced to Bernoulli.
hs_nongaussian_family_symbol <- function(family, n_trials = NULL) {
  if (identical(family$family, "poisson") && identical(family$link, "log")) {
    return("poisson")
  }
  if (identical(family$family, "binomial") && identical(family$link, "logit")) {
    if (!is.null(n_trials) && any(n_trials > 1L)) {
      return("binomial")
    }
    return("bernoulli")
  }
  stop(
    "The opt-in non-Gaussian target fits `poisson(log)` and `binomial(logit)` ",
    "(binary 0/1, or `cbind(successes, failures)` counts) only; `",
    hs_family_label(family),
    "` is not implemented. Other families are planned.",
    call. = FALSE
  )
}

# Resolve the non-Gaussian marginal-method name to the engine's canonical symbol.
# "laplace" (the Laplace approximation; default) and "variational" (the
# variational/ELBO marginal) are accepted, with the DRM-style short spellings
# "la"/"va" as aliases (the engine itself accepts :laplace/:LA and
# :variational/:VA). Both objectives are engine-validated (row V6-LAPLACE/VA).
hs_validate_marginal_method <- function(marginal) {
  if (is.null(marginal)) {
    return("laplace")
  }
  canon <- switch(
    tolower(as.character(marginal)),
    laplace = "laplace",
    la = "laplace",
    variational = "variational",
    va = "variational",
    NULL
  )
  if (is.null(canon)) {
    stop(
      "`engine_control$marginal` must be \"laplace\" (Laplace approximation) or ",
      "\"variational\" (variational/ELBO; aliases \"la\"/\"va\"); got `",
      as.character(marginal),
      "`.",
      call. = FALSE
    )
  }
  canon
}

# Opt-in, experimental non-Gaussian animal model.  The versioned 0.9 transport
# reports only the ratified conditional three-field contract: Poisson has latent
# and count-scale observation h2; logit Bernoulli and common-trial Binomial have
# latent, liability, and observation h2, while heterogeneous trial vectors retain
# the explicit non-scalar observation-scale sentinel.  The Laplace
# objective is a marginal likelihood approximation; the variational objective
# is an ELBO, not a REML or AI-REML claim.
hs_fit_julia_nongaussian_payload <- function(
  payload,
  project = hs_default_julia_project(),
  family = stats::binomial(),
  marginal = "laplace",
  iterations = 200L,
  initial = NULL,
  restart_check = FALSE
) {
  if (!inherits(payload, "hs_bridge_payload")) {
    stop("`payload` must be an internal `hs_bridge_payload`.", call. = FALSE)
  }
  if (is.null(payload$pedigree)) {
    stop(
      "Internal bridge error: the non-Gaussian target requires a pedigree ",
      "animal-model payload.",
      call. = FALSE
    )
  }
  admission <- hs_validate_nongaussian_three_field_v09_admission(
    payload = payload,
    family = family,
    marginal = marginal
  )
  if (!hs_julia_bridge_available(project)) {
    stop(
      "The experimental Julia bridge requires Julia, the `JuliaCall` R ",
      "package, and a local `HSquared.jl` project.",
      call. = FALSE
    )
  }

  n_trials <- admission$n_trials
  family_symbol <- admission$family
  marginal <- admission$method
  iterations <- hs_validate_iterations(iterations)
  initial <- hs_validate_nongaussian_initial(initial)
  restart_check <- hs_validate_restart_check(restart_check)
  hs_julia_setup(project)
  JuliaCall::julia_assign("hsq_y", payload$y)
  JuliaCall::julia_assign("hsq_X", payload$X)
  hs_julia_assign_sparse_csc("hsq_Z", payload$Z)
  JuliaCall::julia_assign("hsq_id", payload$pedigree$id)
  JuliaCall::julia_assign(
    "hsq_sire",
    hs_parent_for_julia(payload$pedigree$sire)
  )
  JuliaCall::julia_assign("hsq_dam", hs_parent_for_julia(payload$pedigree$dam))
  JuliaCall::julia_assign("hsq_family", family_symbol)
  JuliaCall::julia_assign("hsq_marginal", marginal)
  JuliaCall::julia_assign("hsq_iterations", iterations)
  if (identical(family_symbol, "binomial")) {
    n_trials_int <- as.integer(n_trials)
    if (length(unique(n_trials_int)) == 1L) {
      JuliaCall::julia_assign("hsq_n_trials", n_trials_int[[1L]])
    } else {
      JuliaCall::julia_assign("hsq_n_trials", n_trials_int)
    }
  }
  hs_julia_fit(
    JuliaCall::julia_command(hs_nongaussian_three_field_julia_command(
      family_symbol = family_symbol,
      marginal = marginal,
      n_trials = n_trials,
      initial = initial,
      restart_check = restart_check
    )),
    hint = hs_dense_scale_hint
  )

  raw <- JuliaCall::julia_eval("hsq_ng_raw")
  result <- hs_normalize_nongaussian_three_field_v09(raw, payload)
  # The engine echoes the canonical marginal objective it actually ran.
  method_label <- if (identical(result$marginal_method, "variational")) {
    "Variational ELBO"
  } else {
    "Laplace marginal likelihood"
  }
  hs_new_fit(
    spec = list(
      method = method_label,
      family = list(family = family$family, link = family$link),
      target = "nongaussian"
    ),
    payload = payload,
    result = result,
    engine = "HSquared.jl"
  )
}

hs_normalize_nongaussian_result <- function(raw, payload) {
  family <- as.character(raw$family)
  method <- hs_validate_marginal_method(raw$method)
  n_trials <- raw$n_trials
  if (!is.null(n_trials)) {
    n_trials <- as.integer(n_trials)
  }
  fixed_effects <- as.numeric(raw$beta)
  fixed_names <- payload$metadata$fixed_colnames
  if (length(fixed_effects) == length(fixed_names)) {
    names(fixed_effects) <- fixed_names
  }
  animal_bv <- data.frame(
    id = as.character(raw$breeding_ids),
    value = as.numeric(raw$breeding_values),
    stringsAsFactors = FALSE
  )
  converged <- isTRUE(raw$converged)
  result <- list(
    # Legacy envelope only: it reports latent-scale additive variance, but no
    # h2. The separately versioned three-field route carries the narrow A3
    # contract after explicit admission validation.
    variance_components = data.frame(
      component = "animal",
      estimate = as.numeric(raw$sigma_a2),
      stringsAsFactors = FALSE
    ),
    breeding_values = animal_bv,
    random_effects = list(animal = animal_bv),
    fixed_effects = fixed_effects,
    nobs = length(payload$y),
    converged = converged,
    family = family,
    marginal_method = method,
    diagnostics = list(
      target = "nongaussian",
      variance_components = if (identical(method, "variational")) {
        "estimated_variational_reml"
      } else {
        "estimated_laplace_reml"
      },
      engine_family = family,
      marginal_method = method,
      latent_scale = TRUE,
      # The Laplace marginal reports the Laplace-approximate marginal loglik; the
      # variational marginal reports the ELBO (a LOWER BOUND on log p(y)), so
      # logLik/AIC are NOT comparable across the two marginals.
      loglik_kind = if (identical(method, "variational")) {
        "elbo (variational lower bound)"
      } else {
        "laplace marginal loglik"
      },
      heritability_note = paste(
        "The legacy non-Gaussian envelope reports no heritability. The",
        "separately versioned three-field route is available only for its",
        "explicitly admitted experimental contract."
      )
    )
  )
  if (!is.null(n_trials)) {
    result$n_trials <- n_trials
  }
  if (converged) {
    result$loglik <- as.numeric(raw$loglik)
    # The objective value: the Laplace-approximate marginal log-likelihood for
    # `marginal = "laplace"`, or the ELBO (a lower bound) for `"variational"` --
    # see diagnostics$loglik_kind; the two are not comparable across marginals.
    result$loglik_kind <- result$diagnostics$loglik_kind
    # df = fixed effects + the single additive-genetic variance component.
    result$df <- as.integer(ncol(payload$X) + 1L)
  }
  result
}

# Private 0.9 wire-contract helpers.  These are intentionally distinct from
# `hs_normalize_nongaussian_result()` above: the old experimental payload is a
# compatibility surface and may not acquire the three-field semantics by
# accident.  A later, explicitly versioned Julia call will select this
# normalizer by its `nongaussian_three_field_v09` schema tag.
hs_ng09_abort <- function(message) {
  stop("Invalid `nongaussian_three_field_v09` envelope: ", message, call. = FALSE)
}

hs_ng09_required <- function(raw, name) {
  if (!is.list(raw) || !name %in% names(raw)) {
    hs_ng09_abort(paste0("missing required `", name, "` member."))
  }
  raw[[name]]
}

hs_ng09_scalar_number <- function(value, name, nonnegative = FALSE) {
  ok <- is.numeric(value) && length(value) == 1L && is.finite(value)
  if (isTRUE(nonnegative)) {
    ok <- ok && value >= 0
  }
  if (!ok) {
    hs_ng09_abort(paste0("`", name, "` must be a finite numeric scalar",
      if (isTRUE(nonnegative)) " >= 0." else "."))
  }
  as.numeric(value)
}

hs_ng09_exact_number <- function(actual, expected, name) {
  actual <- hs_ng09_scalar_number(actual, name)
  if (!identical(unname(actual), unname(as.numeric(expected)))) {
    hs_ng09_abort(paste0("`", name, "` does not equal its ratified identity."))
  }
  actual
}

hs_ng09_components <- function(raw) {
  components <- hs_ng09_required(raw, "components")
  expected_names <- c("V_A", "V_RE", "V_O")
  if (
    is.list(components) && length(components) == 2L &&
      setequal(names(components), c("names", "values"))
  ) {
    components <- stats::setNames(
      components[["values"]], components[["names"]]
    )
  }
  if (
    !is.numeric(components) ||
      !identical(names(components), expected_names) ||
      length(components) != length(expected_names) ||
      any(!is.finite(components)) || any(components < 0) ||
      sum(components) <= 0
  ) {
    hs_ng09_abort(
      "`components` must be named exactly V_A, V_RE, V_O with finite non-negative values and positive total."
    )
  }
  if (components[["V_RE"]] != 0 || components[["V_O"]] != 0) {
    hs_ng09_abort(
      "`V_RE` and `V_O` are structural zero in the current 0.9 one-additive-effect model."
    )
  }
  as.numeric(stats::setNames(components, expected_names))
}

hs_ng09_intercept <- function(raw) {
  fixed_effects <- hs_ng09_required(raw, "fixed_effects")
  if (
    is.list(fixed_effects) && length(fixed_effects) == 2L &&
      setequal(names(fixed_effects), c("names", "values"))
  ) {
    fixed_effects <- stats::setNames(
      fixed_effects[["values"]],
      fixed_effects[["names"]]
    )
  }
  if (
    !is.numeric(fixed_effects) || length(fixed_effects) != 1L ||
      !identical(names(fixed_effects), "(Intercept)") ||
      !is.finite(fixed_effects)
  ) {
    hs_ng09_abort(
      "`fixed_effects` must contain exactly one finite named `(Intercept)`; 0.9 does not average predictors."
    )
  }
  as.numeric(fixed_effects[[1L]])
}

hs_ng09_loglik_kind <- function(method) {
  if (identical(method, "variational")) {
    return("elbo (variational lower bound)")
  }
  "laplace marginal loglik"
}

hs_ng09_breeding_values <- function(raw) {
  ids <- hs_ng09_required(raw, "breeding_ids")
  values <- hs_ng09_required(raw, "breeding_values")
  if (!is.character(ids) || !is.numeric(values) || length(ids) < 1L ||
    length(ids) != length(values) || anyNA(ids) || any(!nzchar(ids)) ||
    any(!is.finite(values))) {
    hs_ng09_abort(
      "`breeding_ids` and `breeding_values` must be non-empty, aligned, finite vectors."
    )
  }
  data.frame(
    id = as.character(ids), value = as.numeric(values),
    stringsAsFactors = FALSE
  )
}

hs_ng09_converged <- function(raw) {
  converged <- hs_ng09_required(raw, "converged")
  if (!is.logical(converged) || length(converged) != 1L || is.na(converged)) {
    hs_ng09_abort("`converged` must be one non-missing logical value.")
  }
  if (!isTRUE(converged)) {
    hs_ng09_abort(
      "`converged` is FALSE; no three-field fit is returned from an unconverged optimization."
    )
  }
  TRUE
}

# hsquared#222: the R bridge must consume `boundary`, not only `converged`
# (HSquared.jl#342 / HSquared.jl#327). The Julia payload builder
# (`nongaussian_three_field_payload`) already refuses a `boundary = true` fit
# with an `ArgumentError` before this envelope is even built -- that refusal
# is translated into a classed `hsquared_julia_error` by `hs_julia_fit()` at
# the call site, and is the intended contract for THIS route. This check is a
# defense-in-depth backstop: if a future route ever constructs the
# `nongaussian_three_field_v09` envelope without going through that refusing
# builder, a `boundary = TRUE` wire value must not silently pass through as a
# usable point estimate.
hs_ng09_boundary <- function(raw) {
  boundary <- hs_ng09_required(raw, "boundary")
  if (!is.logical(boundary) || length(boundary) != 1L || is.na(boundary)) {
    hs_ng09_abort("`boundary` must be one non-missing logical value.")
  }
  if (isTRUE(boundary)) {
    hs_abort_boundary_refused(
      "the non-Gaussian fit is at its search boundary (boundary = TRUE): ",
      "the estimate sits on the rail of the engine's log-scale search bracket ",
      "(log(initial$sigma_a2) +/- 6), or the two starts disagreed under ",
      "`restart_check = TRUE`; either way it is a function of the search, not ",
      "the data. Retry with a different `initial` in `engine_control` -- ",
      "`initial = list(sigma_a2 = <a value near the expected scale>)` recentres ",
      "the bracket and can move the optimum into its interior. ",
      "`restart_check = TRUE` only makes this check stricter; it cannot clear a ",
      "boundary that is already flagged (HSquared.jl#327)."
    )
  }
  FALSE
}

hs_ng09_heritability_table <- function(result) {
  fields <- "h2_latent"
  labels <- result$h2_latent_label
  estimates <- result$h2_latent
  reasons <- NA_character_
  if ("h2_liability" %in% names(result)) {
    fields <- c(fields, "h2_liability")
    labels <- c(labels, result$h2_liability_label)
    estimates <- c(estimates, result$h2_liability)
    reasons <- c(reasons, NA_character_)
  }
  if ("h2_observation" %in% names(result)) {
    fields <- c(fields, "h2_observation")
    labels <- c(labels, result$h2_observation_label)
    estimates <- c(estimates, result$h2_observation)
    reasons <- c(reasons, result$h2_observation_undefined_reason %||% NA_character_)
  }
  data.frame(
    field = fields,
    label = labels,
    estimate = as.numeric(estimates),
    undefined_reason = reasons,
    stringsAsFactors = FALSE
  )
}

# Assemble, but do not execute, the versioned Julia transport.  Keeping this
# string builder pure makes the every-key-present contract testable without a
# Julia process.  The returned `Dict` deliberately contains explicit `nothing`
# values, which JuliaCall converts to R NULL; absence is therefore detectable
# as a schema error by the normalizer.
hs_nongaussian_three_field_julia_command <- function(
  family_symbol,
  marginal,
  n_trials = NULL,
  initial = NULL,
  restart_check = FALSE
) {
  if (!family_symbol %in% c("poisson", "bernoulli", "binomial")) {
    stop("Invalid v0.9 non-Gaussian engine family.", call. = FALSE)
  }
  marginal <- hs_validate_marginal_method(marginal)
  n_trials_kw <- ""
  if (identical(family_symbol, "binomial")) {
    n_trials <- as.integer(n_trials)
    if (length(n_trials) == 1L) {
      n_trials_kw <- "n_trials = Int(hsq_n_trials), "
    } else {
      n_trials_kw <- "n_trials = Vector{Int}(hsq_n_trials), "
    }
  }
  # hsquared#225: an unsupplied `initial`/`restart_check` omits the keyword
  # entirely so the pre-fix command is reproduced byte for byte -- Julia's own
  # defaults (sigma_a2 = 1.0, restart_check = false) apply unchanged.
  initial_kw <- if (is.null(initial)) {
    ""
  } else {
    paste0(
      "initial = (sigma_a2 = ",
      format(initial, digits = 15, scientific = FALSE, trim = TRUE),
      ",), "
    )
  }
  restart_kw <- if (isTRUE(restart_check)) ", restart_check = true" else ""
  paste0(
    "hsq_ped = HSquared.normalize_pedigree(hsq_id, hsq_sire, hsq_dam); ",
    "hsq_Ainv = HSquared.pedigree_inverse(hsq_ped); ",
    "hsq_fit = HSquared.fit_laplace_reml(",
    "hsq_y, hsq_X, hsq_Z, hsq_Ainv; ",
    "family = Symbol(hsq_family), marginal = Symbol(hsq_marginal), ",
    n_trials_kw,
    initial_kw,
    "ids = hsq_ped.ids, iterations = hsq_iterations",
    restart_kw,
    "); ",
    "hsq_result = HSquared.nongaussian_three_field_payload(",
    "hsq_fit; predictor_variance = 0.0, response_length = length(hsq_y)); ",
    "hsq_ng_raw = Dict(",
    "\"schema\" => hsq_result.schema, ",
    "\"family\" => hsq_result.family, ",
    "\"method\" => hsq_result.method, ",
    "\"loglik\" => hsq_result.loglik, ",
    "\"components\" => Dict(\"names\" => [\"V_A\", \"V_RE\", \"V_O\"], ",
    "\"values\" => [hsq_result.components.V_A, hsq_result.components.V_RE, hsq_result.components.V_O]), ",
    "\"fixed_effects\" => Dict(\"names\" => hsq_result.fixed_effects.names, ",
    "\"values\" => hsq_result.fixed_effects.values), ",
    "\"breeding_ids\" => string.(collect(hsq_fit.ids)), ",
    "\"breeding_values\" => collect(Float64, hsq_fit.breeding_values), ",
    "\"h2_latent\" => hsq_result.h2_latent, ",
    "\"h2_liability\" => hsq_result.h2_liability, ",
    "\"h2_observation\" => hsq_result.h2_observation, ",
    "\"h2_observation_undefined_reason\" => hsq_result.h2_observation_undefined_reason, ",
    "\"n_trials\" => hsq_result.n_trials, ",
    "\"converged\" => hsq_fit.converged, ",
    "\"boundary\" => hsq_fit.boundary, ",
    "\"restart_estimate\" => hsq_fit.restart_estimate);"
  )
}

# Normalize a complete v0.9 three-field envelope.  Exact (zero-tolerance)
# identities are checked at the language boundary, so a transport or formula
# mutation cannot become a different scientific estimand in the R result.
hs_normalize_nongaussian_three_field_v09 <- function(raw, payload) {
  schema <- hs_ng09_required(raw, "schema")
  if (!identical(schema, "nongaussian_three_field_v09")) {
    hs_ng09_abort("`schema` must equal `nongaussian_three_field_v09`.")
  }
  if (!is.list(payload) || is.null(payload$y)) {
    hs_ng09_abort("a payload with response `y` is required for normalization.")
  }

  family <- hs_ng09_required(raw, "family")
  if (!is.character(family) || length(family) != 1L ||
    !family %in% c("poisson", "bernoulli", "binomial")) {
    hs_ng09_abort("`family` must be one of poisson, bernoulli, or binomial.")
  }
  method <- hs_validate_marginal_method(hs_ng09_required(raw, "method"))
  loglik <- hs_ng09_scalar_number(hs_ng09_required(raw, "loglik"), "loglik")
  converged <- hs_ng09_converged(raw)
  boundary <- hs_ng09_boundary(raw)
  components <- hs_ng09_components(raw)
  names(components) <- c("V_A", "V_RE", "V_O")
  mu <- hs_ng09_intercept(raw)
  animal_bv <- hs_ng09_breeding_values(raw)
  v_eta_random <- sum(components)
  h2_latent <- hs_ng09_exact_number(
    hs_ng09_required(raw, "h2_latent"),
    components[["V_A"]] / v_eta_random,
    "h2_latent"
  )
  h2_liability <- hs_ng09_required(raw, "h2_liability")
  h2_observation <- hs_ng09_required(raw, "h2_observation")
  undefined_reason <- hs_ng09_required(
    raw,
    "h2_observation_undefined_reason"
  )
  n_trials <- hs_ng09_required(raw, "n_trials")

  result <- list(
    family = family,
    marginal_method = method,
    loglik = loglik,
    loglik_kind = hs_ng09_loglik_kind(method),
    variance_components = data.frame(
      component = names(components), estimate = as.numeric(components),
      stringsAsFactors = FALSE
    ),
    fixed_effects = stats::setNames(mu, "(Intercept)"),
    nobs = length(payload$y),
    converged = converged,
    boundary = boundary,
    breeding_values = animal_bv,
    random_effects = list(animal = animal_bv),
    h2_latent = h2_latent,
    h2_latent_label = "latent-scale h2 (conditional)"
  )

  if (identical(family, "poisson")) {
    if (!is.null(h2_liability) || !is.null(undefined_reason) || !is.null(n_trials)) {
      hs_ng09_abort(
        "Poisson requires present `nothing` for liability, observation reason, and n_trials."
      )
    }
    expected_observation <- components[["V_A"]] / (
      expm1(v_eta_random) + exp(-(mu + v_eta_random / 2))
    )
    result$h2_observation <- hs_ng09_exact_number(
      h2_observation,
      expected_observation,
      "h2_observation"
    )
    result$h2_observation_label <- "count-scale observation h2 (conditional)"
    result$heritability <- hs_ng09_heritability_table(result)
    return(result)
  }

  expected_liability <- components[["V_A"]] /
    (v_eta_random + pi^2 / 3)
  result$h2_liability <- hs_ng09_exact_number(
    h2_liability,
    expected_liability,
    "h2_liability"
  )
  result$h2_liability_label <- "liability-scale h2 (conditional)"

  if (identical(family, "bernoulli")) {
    if (!is.null(n_trials)) {
      hs_ng09_abort("Bernoulli requires present `nothing` for n_trials.")
    }
    varying_trials <- FALSE
  } else {
    n_trials <- hs_ng09_validate_trials(n_trials, length(payload$y))
    varying_trials <- length(n_trials) > 1L && length(unique(n_trials)) > 1L
  }

  if (isTRUE(varying_trials)) {
    if (!is.numeric(h2_observation) || length(h2_observation) != 1L ||
      !is.nan(h2_observation)) {
      hs_ng09_abort(
        "varying-trial logit `h2_observation` must be literal NaN."
      )
    }
    if (!identical(undefined_reason, "varying_trials_no_scalar_estimand")) {
      hs_ng09_abort(
        paste0(
          "varying-trial `h2_observation_undefined_reason` must equal ",
          "`varying_trials_no_scalar_estimand`."
        )
      )
    }
    result$h2_observation <- NaN
    result$h2_observation_label <- "observation-scale h2 (no scalar estimand)"
    result$h2_observation_undefined_reason <- undefined_reason
  } else {
    if (!is.null(undefined_reason)) {
      hs_ng09_abort(
        "defined logit observation h2 requires present `nothing` for its reason."
      )
    }
    h2_observation <- hs_ng09_scalar_number(
      h2_observation,
      "h2_observation",
      nonnegative = TRUE
    )
    if (h2_observation > 1) {
      hs_ng09_abort("logit `h2_observation` must lie in [0, 1].")
    }
    result$h2_observation <- h2_observation
    result$h2_observation_label <- "observation-scale h2 (conditional)"
  }
  result$heritability <- hs_ng09_heritability_table(result)

  if (!identical(family, "bernoulli")) {
    result$n_trials <- n_trials
  }
  result
}

hs_ng09_validate_trials <- function(n_trials, nobs) {
  if (!is.numeric(n_trials) || length(n_trials) < 1L ||
    any(!is.finite(n_trials)) || any(n_trials != round(n_trials)) ||
    any(n_trials < 1)) {
    hs_ng09_abort("Binomial `n_trials` must be positive integer scalar or vector.")
  }
  n_trials <- as.integer(n_trials)
  if (length(n_trials) == 1L) {
    if (n_trials <= 1L) {
      hs_ng09_abort("Binomial scalar `n_trials` must be greater than one.")
    }
    return(n_trials)
  }
  if (length(n_trials) != nobs) {
    hs_ng09_abort("Binomial `n_trials` vector must have response length.")
  }
  if (all(n_trials == 1L)) {
    hs_ng09_abort("All-one trials are Bernoulli and must not carry n_trials.")
  }
  n_trials
}

# Validate the narrow 0.9 admission boundary before any Julia marshalling.
# This helper has no side effects and intentionally does not select an engine
# target; dispatch stays unavailable until the paired Julia entrypoint exists.
hs_validate_nongaussian_three_field_v09_admission <- function(
  payload,
  family,
  marginal = "laplace",
  predictor_variance = 0,
  weights = NULL,
  dots = list()
) {
  if (!is.list(payload) || is.null(payload$y) || is.null(payload$X)) {
    stop("The v0.9 non-Gaussian admission helper requires payload y and X.", call. = FALSE)
  }
  y <- payload$y
  X <- payload$X
  fixed_names <- payload$metadata$fixed_colnames
  if (!is.matrix(X) || nrow(X) != length(y) || ncol(X) != 1L ||
    !identical(fixed_names, "(Intercept)")) {
    stop(
      "The v0.9 non-Gaussian contract permits exactly one intercept and no fixed predictors.",
      call. = FALSE
    )
  }
  if (!is.numeric(predictor_variance) || length(predictor_variance) != 1L ||
    !is.finite(predictor_variance) || predictor_variance != 0) {
    stop("`predictor_variance` must be exactly zero in the v0.9 contract.", call. = FALSE)
  }
  if (!is.null(weights)) {
    stop("The v0.9 non-Gaussian contract does not accept weights.", call. = FALSE)
  }
  if (!is.list(dots) || length(dots) > 0L) {
    stop("The v0.9 non-Gaussian contract does not accept ... controls.", call. = FALSE)
  }
  if (!inherits(family, "family")) {
    stop("`family` must be an R family object.", call. = FALSE)
  }
  method <- hs_validate_marginal_method(marginal)
  if (!is.numeric(y) || anyNA(y) || any(!is.finite(y))) {
    stop("The v0.9 non-Gaussian response must be numeric and finite.", call. = FALSE)
  }

  if (identical(family$family, "poisson") && identical(family$link, "log")) {
    if (any(y < 0) || any(y != round(y)) || !is.null(payload$n_trials)) {
      stop(
        "Poisson(log) requires non-negative integer counts and no n_trials.",
        call. = FALSE
      )
    }
    return(list(family = "poisson", method = method, n_trials = NULL))
  }
  if (!identical(family$family, "binomial") || !identical(family$link, "logit")) {
    stop(
      "The v0.9 non-Gaussian contract admits only poisson(log) and binomial(logit).",
      call. = FALSE
    )
  }

  n_trials <- payload$n_trials
  if (is.null(n_trials)) {
    if (any(!y %in% c(0, 1))) {
      stop("A one-column binomial response must be binary 0/1, not proportions.", call. = FALSE)
    }
    return(list(family = "bernoulli", method = method, n_trials = NULL))
  }
  if (
    is.numeric(n_trials) && length(n_trials) == length(y) &&
      all(is.finite(n_trials)) && all(n_trials == round(n_trials)) &&
      all(n_trials == 1)
  ) {
    if (any(!y %in% c(0, 1))) {
      stop("All-one binomial trials require binary 0/1 successes.", call. = FALSE)
    }
    return(list(family = "bernoulli", method = method, n_trials = NULL))
  }
  n_trials <- hs_ng09_validate_trials(n_trials, length(y))
  if (any(y < 0) || any(y != round(y)) || any(y > rep(n_trials, length.out = length(y)))) {
    stop(
      "Binomial successes must be non-negative integers no greater than n_trials.",
      call. = FALSE
    )
  }
  if (length(n_trials) == 1L) {
    return(list(family = "binomial", method = method, n_trials = n_trials))
  }
  if (length(unique(n_trials)) == 1L) {
    return(list(
      family = "binomial", method = method,
      n_trials = n_trials[[1L]]
    ))
  }
  list(family = "binomial", method = method, n_trials = n_trials)
}

hs_validate_nongaussian_three_field_v09_dots <- function(dots) {
  if (!is.list(dots) || length(dots) == 0L) {
    return(invisible(TRUE))
  }
  if (!is.null(names(dots)) && "weights" %in% names(dots)) {
    stop("The v0.9 non-Gaussian contract does not accept weights.", call. = FALSE)
  }
  stop("The v0.9 non-Gaussian contract does not accept ... controls.", call. = FALSE)
}

# Opt-in, experimental repeatability (permanent-environment) estimator. Surfaces
# the Julia-owned `HSquared.fit_repeatability_reml()` REML-only optimizer through
# the bridge. The permanent-environment effect shares the animal incidence `Z`
# (the engine carries an identity relationship for it), so the existing payload
# is sufficient. Variance components sigma^2_a and sigma^2_pe are only identifiable with
# repeated records per individual.
hs_fit_julia_repeatability_payload <- function(
  payload,
  project = hs_default_julia_project(),
  initial = c(sigma_a2 = 1, sigma_pe2 = 1, sigma_e2 = 1),
  iterations = 200L,
  max_dense_cells = 1e6,
  scale_method = c("dense", "auto")
) {
  scale_method <- match.arg(scale_method)
  if (!inherits(payload, "hs_bridge_payload")) {
    stop("`payload` must be an internal `hs_bridge_payload`.", call. = FALSE)
  }
  if (!hs_julia_bridge_available(project)) {
    stop(
      "The experimental Julia bridge requires Julia, the `JuliaCall` R ",
      "package, and a local `HSquared.jl` project.",
      call. = FALSE
    )
  }

  initial <- hs_validate_repeatability_initial(initial)
  iterations <- hs_validate_iterations(iterations)
  hs_julia_setup(project)
  hs_julia_assign_payload(payload, initial)
  JuliaCall::julia_assign(
    "hsq_initial_sigma_pe2",
    unname(initial[["sigma_pe2"]])
  )
  JuliaCall::julia_assign("hsq_iterations", iterations)
  JuliaCall::julia_assign("hsq_mdc", hs_validate_max_dense_cells(max_dense_cells))
  if (identical(scale_method, "dense")) {
    # Covered validation-scale route, unchanged: the dense three-component
    # optimizer, gated by the engine on nobs^2 + nanimals^2 <= max_dense_cells.
    hs_julia_fit(
      JuliaCall::julia_command(paste(
        hs_julia_bridge_errors_reset,
        "hsq_ped = HSquared.normalize_pedigree(hsq_id, hsq_sire, hsq_dam);",
        "hsq_Ainv = HSquared.pedigree_inverse(hsq_ped);",
        "hsq_fit = HSquared.fit_repeatability_reml(",
        "hsq_y, hsq_X, hsq_Z, hsq_Ainv;",
        "initial = (sigma_a2 = hsq_initial_sigma_a2,",
        "sigma_pe2 = hsq_initial_sigma_pe2,",
        "sigma_e2 = hsq_initial_sigma_e2),",
        "iterations = hsq_iterations, ids = hsq_ped.ids,",
        "max_dense_cells = Int(hsq_mdc));"
      )),
      hint = hs_repeatability_dense_route_hint
    )
  } else {
    # `scale_method = "auto"`: the SAME animal + permanent-environment model,
    # expressed as the K = 2 independent-block problem the engine already
    # solves, so `fit_multi_effect(:auto)` can take the SPARSE-exact AI-REML
    # route (`sparse_multi_effect_aireml`) and escape the dense ceiling. Block 1
    # is the animal effect carrying A^-1; block 2 is the permanent-environment
    # effect on the SAME incidence `Z` with an identity relationship.
    #
    # `initial`/`iterations` are NOT forwarded here: `fit_multi_effect` does not
    # accept them on this route (HSquared.jl#343). The engine picks its own
    # start, so a supplied `initial` is silently inert -- the R dispatch warns
    # rather than letting the user believe it was honoured.
    hs_julia_fit(
      JuliaCall::julia_command(paste(
        hs_julia_bridge_errors_reset,
        "hsq_ped = HSquared.normalize_pedigree(hsq_id, hsq_sire, hsq_dam);",
        "hsq_Ainv = HSquared.pedigree_inverse(hsq_ped);",
        "hsq_Ipe = spdiagm(0 => ones(size(hsq_Ainv, 1)));",
        "hsq_eff = [(hsq_Z, hsq_Ainv), (hsq_Z, hsq_Ipe)];",
        "hsq_fit = HSquared.fit_multi_effect(",
        "hsq_y, hsq_X, hsq_eff;",
        "method = :auto, verbose = false);"
      ))
    )

    # K-effect asymptotic standard errors (HSquared.jl#352). The engine refuses
    # at a flat/boundary optimum rather than returning NaN, so this is guarded
    # the same way as every other post-fit quantity: the error is RECORDED and
    # surfaced as one warning by hs_julia_surface_bridge_errors(), never
    # swallowed into a silent NA (HSquared.jl#351).
    # ONE information matrix, not three. Each of the three post-fit quantities
    # this route wants -- the variance-component SEs, the ratio (h2 / PE
    # proportion) SEs, and the summed-ratio (repeatability) interval -- is
    # derived from the same asymptotic covariance, and that covariance is a
    # finite-difference Hessian of the REML log-likelihood: the single most
    # expensive post-fit quantity in the engine. Calling the three functions
    # separately built it THREE times. Measured on the great tit animal +
    # permanent-environment fit (11,856 records, 10,937 pedigree): 4.37 s for
    # the three separate calls against 0.17 s through multi_effect_uncertainty,
    # with identical returned values (HSquared.jl#370).
    #
    # The engine function is NEW, so this falls back to the three separate calls
    # when it is absent -- hsquared must keep working against an engine checkout
    # that predates it. The fallback is keyed on the function being UNDEFINED,
    # not on the combined call returning nothing: when it exists and fails, all
    # three quantities are unavailable for the SAME reason (a flat or boundary
    # optimum makes the covariance itself unavailable), so re-calling the three
    # would re-pay the cost only to record the same failure twice more.
    JuliaCall::julia_command(paste(
      "hsq_unc_avail = isdefined(HSquared, :multi_effect_uncertainty);",
      "hsq_unc = if hsq_unc_avail;",
      "try; HSquared.multi_effect_uncertainty(",
      "hsq_y, hsq_X, hsq_eff, hsq_fit.variance_components.sigmas,",
      "hsq_fit.variance_components.sigma_e2; which = 1:2);",
      hs_julia_catch_record("post_fit_uncertainty"),
      "else; nothing; end;",
      "hsq_has_unc = hsq_unc !== nothing;",
      "hsq_vcse = if hsq_has_unc; hsq_unc.variance_component_se;",
      "elseif hsq_unc_avail; nothing;",
      "elseif isdefined(",
      "HSquared, :multi_effect_variance_component_standard_errors);",
      "try; HSquared.multi_effect_variance_component_standard_errors(",
      "hsq_y, hsq_X, hsq_eff, hsq_fit.variance_components.sigmas,",
      "hsq_fit.variance_components.sigma_e2);",
      hs_julia_catch_record("variance_component_standard_errors"),
      "else; nothing; end;",
      "hsq_has_vcse = hsq_vcse !== nothing;",
      "hsq_rse = if hsq_has_unc; hsq_unc.ratio_se;",
      "elseif hsq_unc_avail; nothing;",
      "elseif isdefined(HSquared, :multi_effect_ratio_standard_errors);",
      "try; HSquared.multi_effect_ratio_standard_errors(",
      "hsq_y, hsq_X, hsq_eff, hsq_fit.variance_components.sigmas,",
      "hsq_fit.variance_components.sigma_e2);",
      hs_julia_catch_record("heritability_standard_error"),
      "else; nothing; end;",
      "hsq_has_rse = hsq_rse !== nothing;"
    ))
  }

  # Experimental, opt-in repeatability-coefficient CI (engine row V3-REPEAT-REML,
  # partial). repeatability_interval() takes the raw matrices (not a fit) and
  # refits internally; it throws on a non-positive-definite REML information
  # (flat/boundary optimum) or a boundary t, so the try guard keeps an interval
  # failure from aborting the fit. hsq_has_ri gates the eval so a Julia `nothing`
  # never crosses the bridge.
  # NOT run on the sparse route: repeatability_interval() refits INTERNALLY with
  # the dense estimator, so calling it for a fit that only succeeded because it
  # escaped the dense ceiling would re-impose exactly the ceiling we just
  # escaped. The sparse route therefore carries no repeatability interval.
  if (identical(scale_method, "dense")) {
    JuliaCall::julia_command(paste(
      "hsq_ri = if isdefined(HSquared, :repeatability_interval) &&",
      "applicable(HSquared.repeatability_interval, hsq_y, hsq_X, hsq_Z, hsq_Ainv);",
      "try; HSquared.repeatability_interval(",
      "hsq_y, hsq_X, hsq_Z, hsq_Ainv;",
      "initial = (sigma_a2 = hsq_initial_sigma_a2,",
      "sigma_pe2 = hsq_initial_sigma_pe2,",
      "sigma_e2 = hsq_initial_sigma_e2),",
      "iterations = hsq_iterations, ids = hsq_ped.ids);",
      hs_julia_catch_record("repeatability_interval"),
      "else; nothing; end;",
      "hsq_has_ri = hsq_ri !== nothing;"
    ))
  } else {
    # The sparse route CANNOT call repeatability_interval(): that refits
    # densely, so it would re-impose the ceiling this route exists to escape.
    # multi_effect_sum_ratio_interval() instead forms the same logit-scale
    # delta interval for t = (sigma_a2 + sigma_pe2) / sigma_P from the
    # components already fitted, differentiating the SPARSE loglik. It returns
    # NaN endpoints with boundary = TRUE rather than throwing on a rail, so the
    # `boundary` check below decides whether an interval exists at all.
    #
    # When the combined multi_effect_uncertainty() call above ran, this interval
    # is already in it (computed from the SAME covariance, over which = 1:2), so
    # it is read rather than recomputed. Same value, same fields, one less
    # finite-difference Hessian.
    JuliaCall::julia_command(paste(
      "hsq_ri = if hsq_has_unc; hsq_unc.sum_ratio_interval;",
      "elseif hsq_unc_avail; nothing;",
      "elseif isdefined(HSquared, :multi_effect_sum_ratio_interval);",
      "try; HSquared.multi_effect_sum_ratio_interval(",
      "hsq_y, hsq_X, hsq_eff, hsq_fit.variance_components.sigmas,",
      "hsq_fit.variance_components.sigma_e2; which = 1:2);",
      hs_julia_catch_record("repeatability_interval"),
      "else; nothing; end;",
      "hsq_has_ri = hsq_ri !== nothing && !hsq_ri.boundary;"
    ))
  }

  raw <- if (identical(scale_method, "dense")) {
    JuliaCall::julia_eval(paste(
      "Dict(",
      "\"sigma_a2\" => hsq_fit.variance_components.sigma_a2,",
      "\"sigma_pe2\" => hsq_fit.variance_components.sigma_pe2,",
      "\"sigma_e2\" => hsq_fit.variance_components.sigma_e2,",
      "\"repeatability\" => hsq_fit.repeatability,",
      "\"heritability\" => hsq_fit.heritability,",
      "\"beta\" => collect(Float64, hsq_fit.beta),",
      "\"animal_ids\" => string.(collect(hsq_fit.animal_effects.ids)),",
      "\"animal_values\" => collect(Float64, hsq_fit.animal_effects.values),",
      "\"pe_ids\" => string.(collect(hsq_fit.permanent_effects.ids)),",
      "\"pe_values\" => collect(Float64, hsq_fit.permanent_effects.values),",
      "\"loglik\" => hsq_fit.loglik,",
      "\"loglik_convention\" => hasproperty(hsq_fit, :loglik_convention) ? String(hsq_fit.loglik_convention) : \"unknown\",",
      "\"loglik_full_constant_offset\" => hasproperty(hsq_fit, :loglik_full_constant_offset) ? Float64(hsq_fit.loglik_full_constant_offset) : NaN,",
      "\"loglik_comparable_across_routes\" => hasproperty(hsq_fit, :loglik_comparable_across_routes) ? hsq_fit.loglik_comparable_across_routes : false,",
      "\"converged\" => hsq_fit.converged)"
    ))
  } else {
    # `fit_multi_effect` returns a DIFFERENT shape: variance_components.sigmas
    # is a length-K vector in block order (1 = animal, 2 = permanent) and the
    # per-block effects carry INTEGER row indices into the normalized pedigree,
    # not id strings -- hence the `hsq_ped.ids[...]` lookups. h2 and the
    # repeatability coefficient are not returned by this estimator, so they are
    # formed here from the same components the dense route reports, using the
    # identical definitions (h2 = Va/(Va+Vpe+Ve), R = (Va+Vpe)/(Va+Vpe+Ve)).
    JuliaCall::julia_eval(paste(
      "let s = hsq_fit.variance_components.sigmas,",
      "e = hsq_fit.variance_components.sigma_e2,",
      "tot = hsq_fit.variance_components.sigmas[1] +",
      "hsq_fit.variance_components.sigmas[2] +",
      "hsq_fit.variance_components.sigma_e2;",
      "Dict(",
      "\"sigma_a2\" => s[1],",
      "\"sigma_pe2\" => s[2],",
      "\"sigma_e2\" => e,",
      "\"repeatability\" => (s[1] + s[2]) / tot,",
      "\"heritability\" => s[1] / tot,",
      "\"beta\" => collect(Float64, hsq_fit.beta),",
      "\"animal_ids\" => string.(hsq_ped.ids[collect(hsq_fit.effects[1].ids)]),",
      "\"animal_values\" => collect(Float64, hsq_fit.effects[1].values),",
      "\"pe_ids\" => string.(hsq_ped.ids[collect(hsq_fit.effects[2].ids)]),",
      "\"pe_values\" => collect(Float64, hsq_fit.effects[2].values),",
      "\"loglik\" => hsq_fit.loglik,",
      "\"loglik_convention\" => hasproperty(hsq_fit, :loglik_convention) ? String(hsq_fit.loglik_convention) : \"unknown\",",
      "\"loglik_full_constant_offset\" => hasproperty(hsq_fit, :loglik_full_constant_offset) ? Float64(hsq_fit.loglik_full_constant_offset) : NaN,",
      "\"loglik_comparable_across_routes\" => hasproperty(hsq_fit, :loglik_comparable_across_routes) ? hsq_fit.loglik_comparable_across_routes : false,",
      "\"converged\" => hsq_fit.converged) end"
    ))
  }

  result <- hs_normalize_repeatability_result(raw, payload, scale_method)
  if (!identical(scale_method, "dense")) {
    # Three components, in the SAME order the variance_components table uses,
    # so a reader can bind the two side by side.
    if (isTRUE(JuliaCall::julia_eval("hsq_has_vcse"))) {
      vcse <- JuliaCall::julia_eval(
        "Dict(\"sigmas\" => collect(Float64, hsq_vcse.sigmas),
         \"sigma_e2\" => hsq_vcse.sigma_e2)"
      )
      result$variance_component_se <- data.frame(
        component = c("animal", "permanent", "residual"),
        se = c(as.numeric(vcse$sigmas), as.numeric(vcse$sigma_e2)),
        stringsAsFactors = FALSE
      )
    }
    if (isTRUE(JuliaCall::julia_eval("hsq_has_rse"))) {
      rse <- as.numeric(JuliaCall::julia_eval("collect(Float64, hsq_rse)"))
      # Block 1 is the animal block, so its ratio SE is the h2 SE. Block 2's
      # ratio is the permanent-environment proportion, which is NOT a
      # heritability and is reported separately rather than as an h2 row.
      # A SINGLE NUMERIC, matching the documented contract of
      # `heritability_standard_error()` and what the animal-model route stores
      # (`hs_normalize_*`: `as.numeric(raw$heritability_se)`). Returning a data
      # frame here would make the extractor's return type depend on the target.
      result$heritability_se <- rse[[1L]]
      # The permanent block's ratio is the PE proportion of phenotypic
      # variance, not a heritability, so it gets its own name rather than a
      # second row under h2.
      result$permanent_proportion_se <- rse[[2L]]
    }
  }
  if (isTRUE(JuliaCall::julia_eval("hsq_has_ri"))) {
    # The two routes' interval objects differ: the dense one names the point
    # `repeatability` and carries its own `level`; the sparse one names it
    # `estimate` and takes `level` as an argument (default 0.95). Both are
    # normalized to the SAME R-facing one-row frame, so the extractor does not
    # have to know which estimator produced the fit.
    raw_ri <- if (identical(scale_method, "dense")) {
      JuliaCall::julia_eval(paste(
        "Dict(",
        "\"repeatability\" => hsq_ri.repeatability,",
        "\"lower\" => hsq_ri.lower,",
        "\"upper\" => hsq_ri.upper,",
        "\"level\" => hsq_ri.level,",
        "\"se\" => hsq_ri.se)"
      ))
    } else {
      JuliaCall::julia_eval(paste(
        "Dict(",
        "\"repeatability\" => hsq_ri.estimate,",
        "\"lower\" => hsq_ri.lower,",
        "\"upper\" => hsq_ri.upper,",
        "\"level\" => 0.95,",
        "\"se\" => hsq_ri.se)"
      ))
    }
    result$repeatability_interval <- hs_normalize_repeatability_interval(raw_ri)
  }
  fit <- hs_new_fit(
    spec = list(
      method = "REML",
      family = list(family = payload$family, link = "identity"),
      target = "repeatability"
    ),
    payload = payload,
    result = result,
    engine = "HSquared.jl"
  )
  hs_julia_surface_bridge_errors(fit)
}

hs_normalize_repeatability_result <- function(raw, payload, scale_method = "dense") {
  fixed_effects <- as.numeric(raw$beta)
  fixed_names <- payload$metadata$fixed_colnames
  if (length(fixed_effects) == length(fixed_names)) {
    names(fixed_effects) <- fixed_names
  }
  animal_bv <- data.frame(
    id = as.character(raw$animal_ids),
    value = as.numeric(raw$animal_values),
    stringsAsFactors = FALSE
  )
  list(
    variance_components = data.frame(
      component = c("animal", "permanent", "residual"),
      estimate = c(
        as.numeric(raw$sigma_a2),
        as.numeric(raw$sigma_pe2),
        as.numeric(raw$sigma_e2)
      ),
      stringsAsFactors = FALSE
    ),
    heritability = data.frame(
      term = "animal",
      estimate = as.numeric(raw$heritability)
    ),
    repeatability = data.frame(
      term = "individual",
      estimate = as.numeric(raw$repeatability)
    ),
    breeding_values = animal_bv,
    permanent_effects = data.frame(
      id = as.character(raw$pe_ids),
      value = as.numeric(raw$pe_values),
      stringsAsFactors = FALSE
    ),
    random_effects = list(
      animal = animal_bv,
      permanent = data.frame(
        id = as.character(raw$pe_ids),
        value = as.numeric(raw$pe_values),
        stringsAsFactors = FALSE
      )
    ),
    fixed_effects = fixed_effects,
    loglik = as.numeric(raw$loglik),
    nobs = length(payload$y),
    converged = isTRUE(raw$converged),
    # Provenance names the estimator that ACTUALLY ran. The sparse route is not
    # `fit_repeatability_reml`, so it must not claim that estimator's label --
    # a reader checking `variance_components_source` is checking which code
    # produced the numbers, not which model was requested.
    diagnostics = list(
      variance_components = if (identical(scale_method, "dense")) {
        "estimated_repeatability_reml"
      } else {
        "estimated_repeatability_sparse_multi_effect_aireml"
      },
      scale_method = scale_method,
      # HSquared.jl #365: dense omit-2π vs sparse full-constant. Present when
      # the linked engine exposes the fields; otherwise "unknown"/NA/FALSE.
      loglik_convention = if (!is.null(raw$loglik_convention)) {
        as.character(raw$loglik_convention)
      } else {
        "unknown"
      },
      loglik_full_constant_offset = if (!is.null(raw$loglik_full_constant_offset)) {
        as.numeric(raw$loglik_full_constant_offset)
      } else {
        NA_real_
      },
      loglik_comparable_across_routes = isTRUE(raw$loglik_comparable_across_routes)
    )
  )
}

hs_validate_repeatability_initial <- function(initial) {
  if (
    !is.numeric(initial) ||
      !setequal(names(initial), c("sigma_a2", "sigma_pe2", "sigma_e2"))
  ) {
    stop(
      "`initial` for the repeatability target must be a named numeric vector ",
      "with `sigma_a2`, `sigma_pe2`, and `sigma_e2`.",
      call. = FALSE
    )
  }
  if (any(!is.finite(initial)) || any(initial <= 0)) {
    stop(
      "`initial` variance components must be finite and positive.",
      call. = FALSE
    )
  }
  initial[c("sigma_a2", "sigma_pe2", "sigma_e2")]
}

# Opt-in, experimental two-effect (common-environment) estimator. Surfaces the
# Julia-owned `HSquared.fit_two_effect_reml()` REML-only optimizer: effect 1 is
# the additive-genetic animal effect (Z, pedigree Ainv); effect 2 is the
# common-environment effect (Z2 from the environmental grouping, identity
# relationship). Returns three variance components (animal, common_env,
# residual).
hs_fit_julia_two_effect_payload <- function(
  payload,
  project = hs_default_julia_project(),
  initial = c(sigma_a2 = 1, sigma_c2 = 1, sigma_e2 = 1),
  iterations = 200L
) {
  if (!inherits(payload, "hs_bridge_payload")) {
    stop("`payload` must be an internal `hs_bridge_payload`.", call. = FALSE)
  }
  if (is.null(payload$Z2) || is.null(payload$effect2)) {
    stop(
      "Internal bridge error: the two-effect payload is missing its second ",
      "design matrix.",
      call. = FALSE
    )
  }
  if (!inherits(payload$Z2, "dgCMatrix")) {
    stop(
      "Internal bridge error: the two-effect Z2 must be a sparse dgCMatrix.",
      call. = FALSE
    )
  }
  if (!hs_julia_bridge_available(project)) {
    stop(
      "The experimental Julia bridge requires Julia, the `JuliaCall` R ",
      "package, and a local `HSquared.jl` project.",
      call. = FALSE
    )
  }

  initial <- hs_validate_two_effect_initial(initial)
  iterations <- hs_validate_iterations(iterations)
  hs_julia_setup(project)
  hs_julia_assign_payload(payload, initial)
  hs_julia_assign_sparse_csc("hsq_Z2", payload$Z2)
  JuliaCall::julia_assign("hsq_initial_sigma_c2", unname(initial[["sigma_c2"]]))
  JuliaCall::julia_assign("hsq_iterations", iterations)

  if (identical(payload$effect2$relationship, "pedigree")) {
    # Maternal genetic effect: effect 2 shares the pedigree relationship; its
    # columns are the pedigree animals (dams expressed through Z2).
    ainv2_cmd <- "hsq_Ainv2 = hsq_Ainv;"
    ids2_cmd <- "ids2 = hsq_ped.ids"
  } else {
    # IID effect (e.g. common environment): identity relationship over levels.
    JuliaCall::julia_assign("hsq_env_levels", payload$effect2$levels)
    ainv2_cmd <- paste(
      "hsq_n2 = size(hsq_Z2, 2);",
      "hsq_Ainv2 = sparse(collect(1:hsq_n2), collect(1:hsq_n2),",
      "ones(Float64, hsq_n2), hsq_n2, hsq_n2);"
    )
    ids2_cmd <- "ids2 = hsq_env_levels"
  }

  hs_julia_fit(
    JuliaCall::julia_command(paste(
      hs_julia_bridge_errors_reset,
      "hsq_ped = HSquared.normalize_pedigree(hsq_id, hsq_sire, hsq_dam);",
      "hsq_Ainv = HSquared.pedigree_inverse(hsq_ped);",
      ainv2_cmd,
      "hsq_fit = HSquared.fit_two_effect_reml(",
      "hsq_y, hsq_X, hsq_Z, hsq_Ainv, hsq_Z2, hsq_Ainv2;",
      "initial = (sigma1 = hsq_initial_sigma_a2,",
      "sigma2 = hsq_initial_sigma_c2,",
      "sigma_e2 = hsq_initial_sigma_e2),",
      "iterations = hsq_iterations, ids1 = hsq_ped.ids,",
      ids2_cmd,
      ");"
    )),
    hint = hs_dense_scale_hint
  )

  raw <- JuliaCall::julia_eval(paste(
    "Dict(",
    "\"sigma_a2\" => hsq_fit.variance_components.sigma1,",
    "\"sigma_c2\" => hsq_fit.variance_components.sigma2,",
    "\"sigma_e2\" => hsq_fit.variance_components.sigma_e2,",
    "\"heritability\" => hsq_fit.ratio1,",
    "\"c2\" => hsq_fit.ratio2,",
    "\"beta\" => collect(Float64, hsq_fit.beta),",
    "\"animal_ids\" => string.(collect(hsq_fit.effect1.ids)),",
    "\"animal_values\" => collect(Float64, hsq_fit.effect1.values),",
    "\"env_ids\" => string.(collect(hsq_fit.effect2.ids)),",
    "\"env_values\" => collect(Float64, hsq_fit.effect2.values),",
    "\"loglik\" => hsq_fit.loglik,",
    "\"converged\" => hsq_fit.converged)"
  ))

  # Experimental, opt-in ratio interval (engine `two_effect_ratio_interval`):
  # ratio1 = h2, ratio2 = c2/m2, on the SAME inputs used for the fit. It refits
  # internally and forms the observed information as the finite-difference
  # Hessian of the two-effect REML loglik; on a flat/boundary optimum the sub-
  # information can be non-positive-definite, so the try guard keeps an interval
  # failure from aborting the fit. hsq_has_ci gates the eval so a Julia `nothing`
  # never crosses the bridge. Asymptotic delta-method, NOT coverage-calibrated.
  JuliaCall::julia_command(paste(
    "hsq_ci = if isdefined(HSquared, :two_effect_ratio_interval) &&",
    "applicable(HSquared.two_effect_ratio_interval,",
    "hsq_y, hsq_X, hsq_Z, hsq_Ainv, hsq_Z2, hsq_Ainv2);",
    "try; HSquared.two_effect_ratio_interval(",
    "hsq_y, hsq_X, hsq_Z, hsq_Ainv, hsq_Z2, hsq_Ainv2;",
    "initial = (sigma1 = hsq_initial_sigma_a2,",
    "sigma2 = hsq_initial_sigma_c2,",
    "sigma_e2 = hsq_initial_sigma_e2),",
    "iterations = hsq_iterations, ids1 = hsq_ped.ids,",
    ids2_cmd,
    ");",
    hs_julia_catch_record("two_effect_ratio_interval"),
    "else; nothing; end;",
    "hsq_has_ci = hsq_ci !== nothing;"
  ))

  result <- hs_normalize_two_effect_result(raw, payload)
  if (isTRUE(JuliaCall::julia_eval("hsq_has_ci"))) {
    raw_ci <- JuliaCall::julia_eval(paste(
      "Dict(",
      "\"level\" => hsq_ci.level,",
      "\"r1_estimate\" => hsq_ci.ratio1.estimate,",
      "\"r1_lower\" => hsq_ci.ratio1.lower,",
      "\"r1_upper\" => hsq_ci.ratio1.upper,",
      "\"r1_se\" => hsq_ci.ratio1.se,",
      "\"r1_lower_clamped\" => hsq_ci.ratio1.lower_clamped,",
      "\"r1_upper_clamped\" => hsq_ci.ratio1.upper_clamped,",
      "\"r1_boundary\" => hsq_ci.ratio1.boundary,",
      "\"r2_estimate\" => hsq_ci.ratio2.estimate,",
      "\"r2_lower\" => hsq_ci.ratio2.lower,",
      "\"r2_upper\" => hsq_ci.ratio2.upper,",
      "\"r2_se\" => hsq_ci.ratio2.se,",
      "\"r2_lower_clamped\" => hsq_ci.ratio2.lower_clamped,",
      "\"r2_upper_clamped\" => hsq_ci.ratio2.upper_clamped,",
      "\"r2_boundary\" => hsq_ci.ratio2.boundary)"
    ))
    result <- hs_attach_two_effect_intervals(result, raw_ci, payload)
  }
  fit <- hs_new_fit(
    spec = list(
      method = "REML",
      family = list(family = payload$family, link = "identity"),
      target = "two_effect"
    ),
    payload = payload,
    result = result,
    engine = "HSquared.jl"
  )
  hs_julia_surface_bridge_errors(fit)
}

# ---------------------------------------------------------------------------
# Opt-in correlated direct--maternal model (Phase 4)
# ---------------------------------------------------------------------------
# Fits the 2x2 G_dm (direct-maternal) model through the engine's payload-v2
# correlated-block path. The payload carries a SINGLE block of type="correlated"
# with Z = Zd (record->animal direct incidence) and partner_incidence = Zm
# (record->dam maternal incidence); the engine builds Ainv from the pedigree and
# calls fit_direct_maternal_reml(y, X, Zd, Zm, Ainv). Returns four variance
# components: sigma_ad (direct), sigma_am (maternal/partner), sigma_dm
# (covariance), sigma_e2 (residual), plus the genetic correlation r_am.
# EXPERIMENTAL: asymptotic estimator, NOT coverage-calibrated; a negative r_am
# is real and biologically expected in many livestock traits (antagonistic
# direct-maternal covariance, Willham 1963, 1972). STRICTLY opt-in via
# target = "direct_maternal".
hs_fit_julia_direct_maternal_payload <- function(
  payload,
  project = hs_default_julia_project(),
  initial = NULL,
  iterations = NULL
) {
  if (!inherits(payload, "hs_bridge_payload")) {
    stop("`payload` must be an internal `hs_bridge_payload`.", call. = FALSE)
  }
  # The generic payload from hs_build_bridge_payload will have TWO blocks for a
  # maternal_genetic formula: block1 = animal (pedigree), block2 = maternal
  # (pedigree). We reassemble them as a SINGLE correlated block here, mirroring
  # the engine's payload-v2 section 2 correlated-block schema:
  #   type="correlated", Z=Zd (block1), partner_incidence=Zm (block2), pedigree.
  blocks <- payload$random_effects
  if (
    is.null(blocks) ||
      length(blocks) < 2L ||
      !identical(blocks[[1L]]$type, "pedigree") ||
      !identical(blocks[[2L]]$type, "pedigree")
  ) {
    stop(
      "Internal bridge error: the direct-maternal payload expects two pedigree ",
      "blocks (animal + maternal). Ensure the formula includes ",
      "`animal(1 | id, pedigree = ped) + maternal_genetic(1 | dam)` and ",
      "`target = \"direct_maternal\"`.",
      call. = FALSE
    )
  }
  if (!hs_julia_bridge_available(project)) {
    stop(
      "The experimental Julia bridge requires Julia, the `JuliaCall` R ",
      "package, and a local `HSquared.jl` project.",
      call. = FALSE
    )
  }

  # hsquared#212: `initial`/`iterations` default to NULL (Julia `nothing`) so
  # an unsupplied control reproduces the exact pre-#212-fix call byte for
  # byte.
  if (!is.null(initial)) {
    initial <- hs_validate_direct_maternal_initial(initial)
  }
  if (!is.null(iterations)) {
    iterations <- hs_validate_iterations(iterations)
  }

  hs_julia_setup(project)
  JuliaCall::julia_assign("hsq_y", payload$y)
  JuliaCall::julia_assign("hsq_X", payload$X)
  JuliaCall::julia_assign("hsq_method", payload$method)

  # block1: animal pedigree block -- Z = Zd (record->animal)
  b_animal <- blocks[[1L]]
  # block2: maternal pedigree block -- Z = Zm (record->dam)
  b_maternal <- blocks[[2L]]

  # Zd: record->animal (direct effect incidence)
  hs_julia_assign_sparse_csc("hsq_Zd", b_animal$Z)
  # Zm: record->dam (maternal effect incidence)
  hs_julia_assign_sparse_csc("hsq_Zm", b_maternal$Z)
  JuliaCall::julia_assign("hsq_blkids", as.character(b_animal$ids))

  # Pedigree columns (from the animal block; maternal block shares the same
  # pedigree since dams are pedigree animals)
  JuliaCall::julia_assign("hsq_blkped_id", as.character(b_animal$pedigree$id))
  JuliaCall::julia_assign(
    "hsq_blkped_sire",
    hs_parent_for_julia(b_animal$pedigree$sire)
  )
  JuliaCall::julia_assign(
    "hsq_blkped_dam",
    hs_parent_for_julia(b_animal$pedigree$dam)
  )

  # Build the correlated block Dict and the full payload Dict on the Julia side,
  # matching the frozen payload-v2 correlated-block schema:
  #   type="correlated", Z=Zd, partner_incidence=Zm, partner_name="maternal",
  #   relmat_status="build_in_julia", pedigree=(id,sire,dam), ids
  JuliaCall::julia_command(paste(
    "hsq_blk = Dict{String,Any}(",
    "\"name\" => \"animal\",",
    "\"type\" => \"correlated\",",
    "\"Z\" => hsq_Zd,",
    "\"partner_incidence\" => hsq_Zm,",
    "\"partner_name\" => \"maternal\",",
    "\"relmat_status\" => \"build_in_julia\",",
    "\"pedigree\" => Dict{String,Any}(",
    "\"id\" => hsq_blkped_id,",
    "\"sire\" => hsq_blkped_sire,",
    "\"dam\" => hsq_blkped_dam),",
    "\"ids\" => hsq_blkids);"
  ))
  JuliaCall::julia_command(paste(
    "hsq_payload_dm = Dict{String,Any}(",
    "\"payload_version\" => 2,",
    "\"y\" => hsq_y, \"X\" => hsq_X,",
    "\"method\" => hsq_method,",
    "\"random_effects\" => Any[hsq_blk]);"
  ))
  JuliaCall::julia_command(
    "hsq_parsed_dm = HSquared.parse_payload_v2(hsq_payload_dm);"
  )
  # hsquared#212: forward `initial`/`iterations` to `fit_payload_v2` only when
  # supplied, so the default call (neither supplied) is the exact pre-fix
  # command string, byte for byte.
  dm_fit_kwargs <- character(0)
  if (!is.null(initial)) {
    JuliaCall::julia_assign("hsq_initial_G_dm", initial$G_dm)
    JuliaCall::julia_assign("hsq_initial_sigma_e2_dm", unname(initial$sigma_e2))
    dm_fit_kwargs <- c(
      dm_fit_kwargs,
      "initial = (G_dm = hsq_initial_G_dm, sigma_e2 = hsq_initial_sigma_e2_dm)"
    )
  }
  if (!is.null(iterations)) {
    JuliaCall::julia_assign("hsq_iterations_dm", iterations)
    dm_fit_kwargs <- c(dm_fit_kwargs, "iterations = hsq_iterations_dm")
  }
  dm_kwargs_str <- if (length(dm_fit_kwargs) > 0L) {
    paste0("; ", paste(dm_fit_kwargs, collapse = ", "))
  } else {
    ""
  }
  hs_julia_fit(
    JuliaCall::julia_command(sprintf(
      "hsq_fit_dm = HSquared.fit_payload_v2(hsq_payload_dm%s);",
      dm_kwargs_str
    )),
    hint = hs_dense_scale_hint
  )
  JuliaCall::julia_command(
    "hsq_res_dm = HSquared.result_payload_v2(hsq_fit_dm, hsq_parsed_dm);"
  )

  # Pull the correlated block variance fields and genetic correlation.
  raw <- JuliaCall::julia_eval(paste(
    "let vb = hsq_res_dm.variance_components.blocks[1];",
    "Dict(",
    "\"direct_variance\"  => Float64(vb.direct_variance),",
    "\"partner_variance\" => Float64(vb.partner_variance),",
    "\"covariance\"       => Float64(vb.covariance),",
    "\"correlation\"      => Float64(hsq_fit_dm.genetic_correlation),",
    "\"residual\"         => Float64(hsq_res_dm.variance_components.residual),",
    "\"loglik\"           => Float64(hsq_res_dm.loglik),",
    "\"converged\"        => hsq_res_dm.converged)",
    "end"
  ))
  # BLUPs: direct (animal) and maternal (dam) effects, pulled per-block.
  direct_ids <- JuliaCall::julia_eval(
    "string.(collect(hsq_res_dm.random_effects[1].ids))"
  )
  direct_vals <- JuliaCall::julia_eval(
    "collect(Float64, hsq_res_dm.random_effects[1].values)"
  )
  maternal_ids <- JuliaCall::julia_eval(
    "string.(collect(hsq_res_dm.random_effects[2].ids))"
  )
  maternal_vals <- JuliaCall::julia_eval(
    "collect(Float64, hsq_res_dm.random_effects[2].values)"
  )
  beta <- JuliaCall::julia_eval("collect(Float64, hsq_fit_dm.beta)")

  result <- hs_normalize_direct_maternal_result(
    raw,
    direct_ids,
    direct_vals,
    maternal_ids,
    maternal_vals,
    beta,
    payload
  )

  hs_new_fit(
    spec = list(
      method = "REML",
      family = list(family = payload$family, link = "identity"),
      target = "direct_maternal"
    ),
    payload = payload,
    result = result,
    engine = "HSquared.jl"
  )
}

hs_normalize_direct_maternal_result <- function(
  raw,
  direct_ids,
  direct_vals,
  maternal_ids,
  maternal_vals,
  beta,
  payload
) {
  fixed_names <- payload$metadata$fixed_colnames
  fe <- as.numeric(beta)
  if (length(fe) == length(fixed_names)) {
    names(fe) <- fixed_names
  }

  sigma_ad <- as.numeric(raw$direct_variance)
  sigma_am <- as.numeric(raw$partner_variance)
  sigma_dm <- as.numeric(raw$covariance)
  sigma_e2 <- as.numeric(raw$residual)
  r_am <- as.numeric(raw$correlation)
  converged <- isTRUE(raw$converged)

  sigma_P <- sigma_ad + sigma_am + sigma_dm + sigma_e2
  # Direct narrow-sense heritability: h2_d = sigma_ad / sigma_P
  # sigma_P = sigma_ad + sigma_am + sigma_dm + sigma_e2 = Var(y_i) for a
  # non-inbred base (2*A[i,dam] = 2*(1/2) = 1, so the covariance contributes
  # coefficient 1 to phenotypic variance; Willham 1963, 1972).  sigma_dm is
  # included because it is a legitimate part of Var(y) -- dropping it would
  # misstate the denominator.  h2 is denominator-dependent under maternal
  # effects; see the conditioning_caveat and total_heritability().
  h2_direct <- if (sigma_P > 0) sigma_ad / sigma_P else NA_real_

  direct_bv <- data.frame(
    id = as.character(direct_ids),
    value = as.numeric(direct_vals),
    stringsAsFactors = FALSE
  )
  maternal_bv <- data.frame(
    id = as.character(maternal_ids),
    value = as.numeric(maternal_vals),
    stringsAsFactors = FALSE
  )

  list(
    variance_components = data.frame(
      component = c("direct", "maternal", "covariance", "residual"),
      estimate = c(sigma_ad, sigma_am, sigma_dm, sigma_e2),
      stringsAsFactors = FALSE
    ),
    # A labelled heritability: direct h2 only. heritability() on a
    # direct_maternal fit returns this with an interpretation fence
    # (NOT a bare scalar).
    heritability = data.frame(
      term = "direct",
      estimate = h2_direct,
      stringsAsFactors = FALSE
    ),
    # The genetic correlation r_am between direct and maternal effects.
    # A NEGATIVE value is real and expected in many livestock traits
    # (antagonistic direct-maternal covariance; Willham 1963, 1972).
    # Never conflate with a multivariate trait-to-trait genetic correlation.
    genetic_correlation = data.frame(
      term_1 = "direct",
      term_2 = "maternal",
      estimate = r_am,
      stringsAsFactors = FALSE
    ),
    # Explicitly labelled direct and partner variance fields for accessors.
    direct_variance = sigma_ad,
    partner_variance = sigma_am,
    covariance = sigma_dm,
    breeding_values = direct_bv,
    random_effects = list(
      animal = direct_bv,
      maternal = maternal_bv
    ),
    maternal_effects = maternal_bv,
    fixed_effects = fe,
    loglik = if (converged) as.numeric(raw$loglik) else NA_real_,
    nobs = length(payload$y),
    converged = converged,
    diagnostics = list(
      target = "direct_maternal",
      variance_components = "estimated_direct_maternal_reml",
      optimizer_status = if (converged) "converged" else "not_converged",
      conditioning_caveat = paste(
        "Experimental direct-maternal 2x2 G_dm REML estimator (Phase 4).",
        "Asymptotic delta-method intervals are NOT coverage-calibrated.",
        "A negative genetic correlation (r_am) is real and expected in many",
        "livestock traits (antagonistic direct-maternal covariance;",
        "Willham 1963, 1972). Identifiability requires multiple offspring per dam",
        "with sires recorded; shallow pedigrees produce boundary G_dm.",
        "sigma^2_P = sigma_ad + sigma_am + sigma_dm + sigma_e2 (Willham 1972);",
        "h2 is denominator-dependent under maternal effects. ASReml/BLUPF90/",
        "WOMBAT report raw components and leave sigma^2_P to the user; sommer/",
        "MCMCglmm h2 depends on the user's chosen denominator - compare",
        "(co)variance components, not h2 values, across software.",
        "Use validate = TRUE to inspect the contract before fitting."
      )
    )
  )
}

# Fit an arbitrary-N independent-random-effect model (animal + >= 2 i.i.d.
# blocks) through the engine's payload-v2 entry points. Rather than re-deriving
# the `(Z_i, Ainv_i)` effects vector in R, this assembles the block-structured
# `random_effects` payload on the Julia side and calls the exported
# `HSquared.fit_payload_v2(payload)` (which parses the block list and routes
# K >= 3 independent blocks to `fit_multi_effect_reml`) and
# `HSquared.result_payload_v2(fit, parsed)` (the block-structured result).
# After the fit, an EXPERIMENTAL per-component ratio interval is attached from
# the engine's `multi_effect_ratio_interval` (guarded try, mirroring the two-
# effect `hs_attach_two_effect_intervals` path): the ANIMAL block's ratio
# populates `heritability_interval` (so `heritability_interval()` resolves for a
# K >= 3 fit) and each block's ratio populates a per-component field. Asymptotic
# logit delta-method, NOT coverage-calibrated (engine row `V3-NEFFECT-REML`).
hs_fit_julia_n_effect_payload <- function(
  payload,
  project = hs_default_julia_project(),
  scale_method = c("dense", "auto"),
  initial = NULL,
  iterations = NULL
) {
  scale_method <- match.arg(scale_method)
  if (!inherits(payload, "hs_bridge_payload")) {
    stop("`payload` must be an internal `hs_bridge_payload`.", call. = FALSE)
  }
  blocks <- payload$random_effects
  # K >= 2, not K >= 3. The old three-block floor was an R-side accident, not a
  # design: every engine entry point (`fit_multi_effect_reml`,
  # `fit_sparse_multi_effect_aireml`, `fit_multi_effect`) asserts only K >= 1.
  # The floor made `animal + ONE i.i.d. effect` unreachable, which is the
  # standard repeated-measures animal model, so users saw "one random effect,
  # or at least three" and V_A silently absorbed the second component
  # (HSquared.jl#352).
  if (is.null(blocks) || length(blocks) < 2L) {
    stop(
      "Internal bridge error: the multi-effect payload needs at least two ",
      "random-effect blocks (animal + one or more i.i.d. effects).",
      call. = FALSE
    )
  }
  if (!hs_julia_bridge_available(project)) {
    stop(
      "The experimental Julia bridge requires Julia, the `JuliaCall` R ",
      "package, and a local `HSquared.jl` project.",
      call. = FALSE
    )
  }

  # hsquared#212: `initial`/`iterations` default to NULL (Julia `nothing`) so
  # an unsupplied control reproduces the exact pre-#212-fix call byte for
  # byte. When supplied, `initial` must be a length K+1 vector (one value per
  # block, in formula order, plus the residual).
  if (!is.null(initial)) {
    initial <- hs_validate_multi_effect_initial(initial, length(blocks) + 1L)
  }
  if (!is.null(iterations)) {
    iterations <- hs_validate_iterations(iterations)
  }

  hs_julia_setup(project)
  JuliaCall::julia_assign("hsq_y", payload$y)
  JuliaCall::julia_assign("hsq_X", payload$X)
  JuliaCall::julia_assign("hsq_method", payload$method)

  # Assign each block's engine inputs and build a Julia Dict per block, matching
  # the payload-v2 section 2 field table read by `parse_payload_v2`:
  #   pedigree block: type="pedigree", relmat_status="build_in_julia",
  #                   pedigree=(id,sire,dam), ids
  #   iid block:      type="iid",      relmat_status="identity", ids
  block_syms <- character(length(blocks))
  for (i in seq_along(blocks)) {
    b <- blocks[[i]]
    z_name <- paste0("hsq_blkZ_", i)
    hs_julia_assign_sparse_csc(z_name, b$Z)
    ids_name <- paste0("hsq_blkids_", i)
    JuliaCall::julia_assign(ids_name, as.character(b$ids))
    blk_sym <- paste0("hsq_blk_", i)
    block_syms[[i]] <- blk_sym
    if (identical(b$type, "pedigree")) {
      JuliaCall::julia_assign(
        paste0("hsq_blkped_id_", i),
        as.character(b$pedigree$id)
      )
      JuliaCall::julia_assign(
        paste0("hsq_blkped_sire_", i),
        hs_parent_for_julia(b$pedigree$sire)
      )
      JuliaCall::julia_assign(
        paste0("hsq_blkped_dam_", i),
        hs_parent_for_julia(b$pedigree$dam)
      )
      JuliaCall::julia_command(sprintf(
        paste0(
          "%s = Dict{String,Any}(\"name\" => \"%s\", \"type\" => \"pedigree\", ",
          "\"Z\" => %s, \"relmat_status\" => \"build_in_julia\", ",
          "\"pedigree\" => Dict{String,Any}(\"id\" => hsq_blkped_id_%d, ",
          "\"sire\" => hsq_blkped_sire_%d, \"dam\" => hsq_blkped_dam_%d), ",
          "\"ids\" => %s);"
        ),
        blk_sym,
        b$name,
        z_name,
        i,
        i,
        i,
        ids_name
      ))
    } else {
      # iid block: identity relationship, no pedigree.
      JuliaCall::julia_command(sprintf(
        paste0(
          "%s = Dict{String,Any}(\"name\" => \"%s\", \"type\" => \"iid\", ",
          "\"Z\" => %s, \"relmat_status\" => \"identity\", \"ids\" => %s);"
        ),
        blk_sym,
        b$name,
        z_name,
        ids_name
      ))
    }
  }

  JuliaCall::julia_command(sprintf(
    "hsq_blocks = Any[%s];",
    paste(block_syms, collapse = ", ")
  ))
  JuliaCall::julia_command(paste(
    "hsq_payload = Dict{String,Any}(\"payload_version\" => 2, ",
    "\"y\" => hsq_y, \"X\" => hsq_X, \"method\" => hsq_method, ",
    "\"random_effects\" => hsq_blocks);"
  ))

  # fit_payload_v2 parses the block list and dispatches K >= 3 independent
  # blocks. scale_method = "dense" (DEFAULT) uses the dense fit_multi_effect_reml
  # (the covered validation-scale path). scale_method = "auto" routes through the
  # engine's fit_multi_effect(:auto): the SPARSE-EXACT AI-REML (reduces exactly to
  # the dense optimum) at validation scale, and the EXPERIMENTAL matrix-free
  # Monte-Carlo fit for large problems the dense factorization cannot reach.
  # result_payload_v2 builds the block-structured result from the same parse.
  JuliaCall::julia_command(
    "hsq_parsed = HSquared.parse_payload_v2(hsq_payload);"
  )
  # hsquared#212: forward `initial`/`iterations` to `fit_payload_v2` (and,
  # below, to the `multi_effect_ratio_interval` refit on the SAME inputs)
  # only when supplied, so the default call (neither supplied) is the exact
  # pre-fix command string, byte for byte.
  ne_extra_kwargs <- character(0)
  if (!is.null(initial)) {
    JuliaCall::julia_assign("hsq_initial_ne", initial)
    ne_extra_kwargs <- c(ne_extra_kwargs, "initial = hsq_initial_ne")
  }
  if (!is.null(iterations)) {
    JuliaCall::julia_assign("hsq_iterations_ne", iterations)
    ne_extra_kwargs <- c(ne_extra_kwargs, "iterations = hsq_iterations_ne")
  }
  ne_fit_kwargs <- c(
    sprintf("scale_method = :%s", scale_method),
    ne_extra_kwargs
  )
  hs_julia_fit(
    JuliaCall::julia_command(sprintf(
      "hsq_fit = HSquared.fit_payload_v2(hsq_payload; %s);",
      paste(ne_fit_kwargs, collapse = ", ")
    )),
    hint = hs_dense_scale_hint
  )
  JuliaCall::julia_command(paste(
    hs_julia_bridge_errors_reset,
    "hsq_result = HSquared.result_payload_v2(hsq_fit, hsq_parsed);"
  ))

  raw <- JuliaCall::julia_eval(paste(
    "Dict(",
    "\"residual\" => hsq_result.variance_components.residual,",
    "\"block_names\" => [b.name for b in hsq_result.variance_components.blocks],",
    "\"block_types\" => [b.type for b in hsq_result.variance_components.blocks],",
    "\"block_variances\" => [Float64(b.variance) for b in ",
    "hsq_result.variance_components.blocks],",
    "\"re_names\" => [r.name for r in hsq_result.random_effects],",
    "\"beta\" => collect(Float64, hsq_fit.beta),",
    "\"loglik\" => hsq_result.loglik,",
    "\"converged\" => hsq_result.converged)"
  ))
  # Per-block BLUP ids/values (ragged), pulled one block at a time so JuliaCall
  # returns clean per-block vectors rather than a jagged nested structure.
  n_blocks <- length(raw$block_names)
  re_values <- vector("list", n_blocks)
  re_ids <- vector("list", n_blocks)
  for (i in seq_len(n_blocks)) {
    re_ids[[i]] <- JuliaCall::julia_eval(sprintf(
      "string.(collect(hsq_result.random_effects[%d].ids))",
      i
    ))
    re_values[[i]] <- JuliaCall::julia_eval(sprintf(
      "collect(Float64, hsq_result.random_effects[%d].values)",
      i
    ))
  }

  result <- hs_normalize_n_effect_result(raw, re_ids, re_values, payload)

  # Experimental, opt-in per-component ratio interval (engine
  # `multi_effect_ratio_interval`) on the SAME inputs used for the fit. The
  # already-parsed `hsq_parsed` carries the resolved (Z_i, Ainv_i) effects and
  # per-block ids; we rebuild the identical `effects` vector the multi-effect
  # dispatch uses (bridge_payload_v2.jl `_dispatch_fit` :multi_effect) so the
  # interval refit sees exactly the fitted model. The interval refits internally
  # and forms the observed information as the finite-difference Hessian of the
  # K-effect REML loglik; on a flat/boundary optimum the sub-information can be
  # non-positive-definite, so the try guard keeps an interval failure from
  # aborting the fit. hsq_has_nci gates the eval so a Julia `nothing` never
  # crosses the bridge. Asymptotic delta-method, NOT coverage-calibrated.
  # hsquared#212: forward the SAME `initial`/`iterations` used for the main
  # fit (mirroring `two_effect_ratio_interval`'s pattern above) so the
  # interval refit is not silently ignoring them independently of the
  # main-fit defect fixed above.
  nci_extra_kwargs <- if (length(ne_extra_kwargs) > 0L) {
    paste0(", ", paste(ne_extra_kwargs, collapse = ", "))
  } else {
    ""
  }
  JuliaCall::julia_command(paste(
    "hsq_nci = if isdefined(HSquared, :multi_effect_ratio_interval);",
    "try;",
    "hsq_neff = [(Matrix{Float64}(b.Z), Matrix{Float64}(b.relmat_inverse))",
    "for b in hsq_parsed.blocks];",
    "hsq_nids = [b.ids for b in hsq_parsed.blocks];",
    "HSquared.multi_effect_ratio_interval(",
    sprintf(
      "hsq_parsed.y, hsq_parsed.X, hsq_neff; ids = hsq_nids%s);",
      nci_extra_kwargs
    ),
    hs_julia_catch_record("multi_effect_ratio_interval"),
    "else; nothing; end;",
    "hsq_has_nci = hsq_nci !== nothing;"
  ))
  if (isTRUE(JuliaCall::julia_eval("hsq_has_nci"))) {
    raw_nci <- JuliaCall::julia_eval(paste(
      "Dict(",
      "\"level\" => hsq_nci.level,",
      "\"converged\" => hsq_nci.converged,",
      "\"estimate\" => [Float64(r.estimate) for r in hsq_nci.ratios],",
      "\"lower\" => [Float64(r.lower) for r in hsq_nci.ratios],",
      "\"upper\" => [Float64(r.upper) for r in hsq_nci.ratios],",
      "\"se\" => [Float64(r.se) for r in hsq_nci.ratios],",
      "\"lower_clamped\" => [r.lower_clamped for r in hsq_nci.ratios],",
      "\"upper_clamped\" => [r.upper_clamped for r in hsq_nci.ratios],",
      "\"boundary\" => [r.boundary for r in hsq_nci.ratios])"
    ))
    result <- hs_attach_n_effect_intervals(result, raw_nci, raw)
  }

  fit <- hs_new_fit(
    spec = list(
      method = "REML",
      family = list(family = payload$family, link = "identity"),
      target = "multi_effect"
    ),
    payload = payload,
    result = result,
    engine = "HSquared.jl"
  )
  hs_julia_surface_bridge_errors(fit)
}

hs_normalize_two_effect_result <- function(raw, payload) {
  fixed_effects <- as.numeric(raw$beta)
  fixed_names <- payload$metadata$fixed_colnames
  if (length(fixed_effects) == length(fixed_names)) {
    names(fixed_effects) <- fixed_names
  }
  animal_bv <- data.frame(
    id = as.character(raw$animal_ids),
    value = as.numeric(raw$animal_values),
    stringsAsFactors = FALSE
  )
  type2 <- payload$effect2$type
  second_eff <- data.frame(
    id = as.character(raw$env_ids),
    value = as.numeric(raw$env_values),
    stringsAsFactors = FALSE
  )
  result <- list(
    variance_components = data.frame(
      component = c("animal", type2, "residual"),
      estimate = c(
        as.numeric(raw$sigma_a2),
        as.numeric(raw$sigma_c2),
        as.numeric(raw$sigma_e2)
      ),
      stringsAsFactors = FALSE
    ),
    heritability = data.frame(
      term = "animal",
      estimate = as.numeric(raw$heritability)
    ),
    breeding_values = animal_bv,
    random_effects = list(animal = animal_bv),
    fixed_effects = fixed_effects,
    loglik = as.numeric(raw$loglik),
    nobs = length(payload$y),
    converged = isTRUE(raw$converged),
    diagnostics = list(variance_components = "estimated_two_effect_reml")
  )
  result$random_effects[[type2]] <- second_eff
  if (identical(type2, "common_env")) {
    result$common_env_effects <- second_eff
    result$common_env_proportion <- data.frame(
      term = "common_env",
      estimate = as.numeric(raw$c2)
    )
  } else {
    result$maternal_effects <- second_eff
    result$maternal_proportion <- data.frame(
      term = "maternal_genetic",
      estimate = as.numeric(raw$c2)
    )
  }
  result
}

# Normalize the block-structured `result_payload_v2` multi-effect result into the
# flat `hsquared_fit` result shape read by the S3 extractors. `raw` carries the
# per-block variances (aligned vectors `block_names`/`block_types`/
# `block_variances`), the residual, `beta`, `loglik`, and `converged`; `re_ids` /
# `re_values` are the ragged per-block BLUP ids/values.
#
# Falconer fence: `heritability` is the narrow-sense h2 of the ANIMAL (pedigree)
# block only -- its additive-genetic variance over the TOTAL phenotypic variance
# (the sum of ALL block variances plus the residual). Every other block's
# variance ratio (`block_variance / total`) is a variance-explained proportion,
# NOT a heritability. This normalizes the POINT ESTIMATES; the caller attaches
# the experimental ratio interval separately via `hs_attach_n_effect_intervals`.
hs_normalize_n_effect_result <- function(raw, re_ids, re_values, payload) {
  fixed_effects <- as.numeric(raw$beta)
  fixed_names <- payload$metadata$fixed_colnames
  if (length(fixed_effects) == length(fixed_names)) {
    names(fixed_effects) <- fixed_names
  }

  block_names <- as.character(raw$block_names)
  block_variances <- as.numeric(raw$block_variances)
  residual <- as.numeric(raw$residual)
  total <- sum(block_variances) + residual

  variance_components <- data.frame(
    component = c(block_names, "residual"),
    estimate = c(block_variances, residual),
    stringsAsFactors = FALSE
  )

  # Per-block random effects (BLUPs), keyed by block name.
  random_effects <- list()
  for (i in seq_along(block_names)) {
    random_effects[[block_names[[i]]]] <- data.frame(
      id = as.character(re_ids[[i]]),
      value = as.numeric(re_values[[i]]),
      stringsAsFactors = FALSE
    )
  }

  animal_idx <- match("animal", block_names)
  if (is.na(animal_idx)) {
    stop(
      "Internal bridge error: the multi-effect result has no `animal` block to ",
      "form a heritability from.",
      call. = FALSE
    )
  }
  heritability <- data.frame(
    term = "animal",
    estimate = block_variances[[animal_idx]] / total,
    stringsAsFactors = FALSE
  )

  # Per-block variance-explained proportions (NOT heritabilities except for the
  # animal block). Surfaced under a distinct field so extractors never confuse a
  # proportion with h2.
  variance_ratios <- data.frame(
    component = block_names,
    estimate = block_variances / total,
    stringsAsFactors = FALSE
  )

  animal_bv <- random_effects[["animal"]]
  list(
    variance_components = variance_components,
    variance_ratios = variance_ratios,
    heritability = heritability,
    breeding_values = animal_bv,
    random_effects = random_effects,
    fixed_effects = fixed_effects,
    loglik = as.numeric(raw$loglik),
    nobs = length(payload$y),
    converged = isTRUE(raw$converged),
    diagnostics = list(variance_components = "estimated_multi_effect_reml")
  )
}

# Attach the multi-effect (K >= 3) per-component ratio intervals to the result.
# `raw_nci` carries the engine `multi_effect_ratio_interval` output as aligned
# per-block vectors (`estimate`/`lower`/`upper`/`se`/`*_clamped`/`boundary`), in
# the SAME block order as `raw$block_names` (both come from the one
# `parse_payload_v2` parse). The ANIMAL block's ratio populates
# `heritability_interval` (one-row shape identical to the univariate delta CI, so
# `heritability_interval()` resolves identically on a K >= 3 fit); EVERY block's
# ratio (animal included) is also surfaced under `variance_ratio_intervals`, a
# named list keyed by block name, so each variance-explained proportion has its
# matching interval. Falconer fence: only the animal entry is a heritability; the
# others are variance-explained-proportion intervals, NOT heritabilities.
# Boundary-flagged (sigma -> 0) components arrive with NaN bounds from the engine
# and are normalized to NA (not a spurious CI).
hs_attach_n_effect_intervals <- function(result, raw_nci, raw) {
  level <- as.numeric(raw_nci$level)
  block_names <- as.character(raw$block_names)
  k <- length(block_names)

  one_row <- function(i, with_method) {
    hs_normalize_two_effect_ratio_interval(
      estimate = raw_nci$estimate[[i]],
      lower = raw_nci$lower[[i]],
      upper = raw_nci$upper[[i]],
      se = raw_nci$se[[i]],
      lower_clamped = raw_nci$lower_clamped[[i]],
      upper_clamped = raw_nci$upper_clamped[[i]],
      boundary = raw_nci$boundary[[i]],
      level = level,
      with_method = with_method
    )
  }

  ratio_intervals <- list()
  for (i in seq_len(k)) {
    ratio_intervals[[block_names[[i]]]] <- one_row(i, with_method = FALSE)
  }
  result$variance_ratio_intervals <- ratio_intervals

  animal_idx <- match("animal", block_names)
  if (!is.na(animal_idx)) {
    result$heritability_interval <- one_row(animal_idx, with_method = TRUE)
  }
  result
}

# Attach the two-effect ratio intervals to the result. ratio1 (h2) populates
# `heritability_interval` (same one-row shape as the univariate delta CI, so
# heritability_interval() resolves identically on a two-effect fit); ratio2
# (c2/m2) populates the second-effect interval field, keyed by the second-effect
# type (common_env -> common_env_proportion_interval, otherwise ->
# maternal_proportion_interval). Boundary-flagged (sigma -> 0) components arrive
# with NaN bounds from the engine and are normalized to NA (not a spurious CI).
hs_attach_two_effect_intervals <- function(result, raw_ci, payload) {
  level <- as.numeric(raw_ci$level)
  result$heritability_interval <- hs_normalize_two_effect_ratio_interval(
    estimate = raw_ci$r1_estimate,
    lower = raw_ci$r1_lower,
    upper = raw_ci$r1_upper,
    se = raw_ci$r1_se,
    lower_clamped = raw_ci$r1_lower_clamped,
    upper_clamped = raw_ci$r1_upper_clamped,
    boundary = raw_ci$r1_boundary,
    level = level,
    with_method = TRUE
  )
  ratio2 <- hs_normalize_two_effect_ratio_interval(
    estimate = raw_ci$r2_estimate,
    lower = raw_ci$r2_lower,
    upper = raw_ci$r2_upper,
    se = raw_ci$r2_se,
    lower_clamped = raw_ci$r2_lower_clamped,
    upper_clamped = raw_ci$r2_upper_clamped,
    boundary = raw_ci$r2_boundary,
    level = level,
    with_method = FALSE
  )
  if (identical(payload$effect2$type, "common_env")) {
    result$common_env_proportion_interval <- ratio2
  } else {
    result$maternal_proportion_interval <- ratio2
  }
  result
}

# Normalize one ratio's interval from `two_effect_ratio_interval` into a one-row
# data frame. The engine reports NaN bounds when a component is on the boundary
# (sigma -> 0); those become NA here. `with_method = TRUE` appends
# `method = "delta"` so ratio1 (h2) matches the univariate heritability_interval
# shape; ratio2 (c2/m2) omits it (delta is the only method for the ratio).
hs_normalize_two_effect_ratio_interval <- function(
  estimate,
  lower,
  upper,
  se,
  lower_clamped,
  upper_clamped,
  boundary,
  level,
  with_method
) {
  na_if_nan <- function(x) {
    x <- as.numeric(x)
    if (length(x) != 1L || is.nan(x)) NA_real_ else x
  }
  out <- data.frame(
    estimate = na_if_nan(estimate),
    lower = na_if_nan(lower),
    upper = na_if_nan(upper),
    level = as.numeric(level),
    se = na_if_nan(se),
    lower_clamped = isTRUE(as.logical(lower_clamped)),
    upper_clamped = isTRUE(as.logical(upper_clamped)),
    boundary = isTRUE(as.logical(boundary)),
    stringsAsFactors = FALSE
  )
  if (isTRUE(with_method)) {
    out$method <- "delta"
  }
  out
}

hs_fit_julia_multivariate_payload <- function(
  payload,
  project = hs_default_julia_project(),
  initial = NULL,
  iterations = 2000L,
  genetic_structure = "unstructured"
) {
  if (!inherits(payload, "hs_bridge_payload")) {
    stop("`payload` must be an internal `hs_bridge_payload`.", call. = FALSE)
  }
  if (is.null(payload$Y) || !is.matrix(payload$Y)) {
    stop(
      "Internal bridge error: the multivariate payload is missing its `Y` ",
      "response matrix.",
      call. = FALSE
    )
  }
  if (is.null(payload$pedigree)) {
    stop(
      "Internal bridge error: the multivariate target requires a pedigree ",
      "animal-model payload.",
      call. = FALSE
    )
  }
  if (!hs_julia_bridge_available(project)) {
    stop(
      "The experimental Julia bridge requires Julia, the `JuliaCall` R ",
      "package, and a local `HSquared.jl` project.",
      call. = FALSE
    )
  }

  ntraits <- ncol(payload$Y)
  initial <- hs_validate_multivariate_initial(initial, ntraits)
  iterations <- hs_validate_iterations(iterations)
  traits <- payload$metadata$trait_names %||% colnames(payload$Y)
  if (is.null(traits)) {
    traits <- paste0("trait", seq_len(ntraits))
  }

  hs_julia_setup(project)
  # JuliaCall 0.17.6 / R 4.6 delivers R NA_real_ as Julia Missing, not NaN.
  # Convert here so the engine receives the documented NaN missing-trait
  # sentinel (it also accepts Missing, but isnan-based round-trips through
  # julia_eval segfault inside Rcpp precious-preserve on Missing arrays).
  JuliaCall::julia_assign("hsq_Y", hs_y_matrix_for_julia(payload$Y))
  JuliaCall::julia_assign("hsq_X", payload$X)
  hs_julia_assign_sparse_csc("hsq_Z", payload$Z)
  JuliaCall::julia_assign("hsq_id", payload$pedigree$id)
  JuliaCall::julia_assign(
    "hsq_sire",
    hs_parent_for_julia(payload$pedigree$sire)
  )
  JuliaCall::julia_assign("hsq_dam", hs_parent_for_julia(payload$pedigree$dam))
  JuliaCall::julia_assign("hsq_traits", as.character(traits))
  JuliaCall::julia_assign("hsq_initial_G0", initial$G0)
  JuliaCall::julia_assign("hsq_initial_R0", initial$R0)
  JuliaCall::julia_assign("hsq_iterations", iterations)
  JuliaCall::julia_assign(
    "hsq_genetic_structure",
    as.character(genetic_structure)
  )
  hs_julia_fit(
    JuliaCall::julia_command(paste(
    hs_julia_bridge_errors_reset,
    "hsq_ped = HSquared.normalize_pedigree(hsq_id, hsq_sire, hsq_dam);",
    "hsq_Ainv = HSquared.pedigree_inverse(hsq_ped);",
    "hsq_fit = HSquared.fit_multivariate_reml(",
    "hsq_Y, hsq_X, hsq_Z, hsq_Ainv;",
    "initial = (G0 = hsq_initial_G0, R0 = hsq_initial_R0),",
    "iterations = hsq_iterations, ids = hsq_ped.ids, traits = hsq_traits,",
    "genetic_structure = Symbol(hsq_genetic_structure));",
    "hsq_mv_raw = Dict(",
    "\"genetic_covariance\" => Matrix{Float64}(hsq_fit.genetic_covariance),",
    "\"residual_covariance\" => Matrix{Float64}(hsq_fit.residual_covariance),",
    "\"genetic_correlation\" => Matrix{Float64}(hsq_fit.genetic_correlation),",
    "\"residual_correlation\" => Matrix{Float64}(hsq_fit.residual_correlation),",
    "\"heritability\" => collect(Float64, hsq_fit.heritability),",
    "\"beta\" => Matrix{Float64}(hsq_fit.beta),",
    "\"breeding_ids\" => string.(collect(hsq_fit.breeding_values.ids)),",
    "\"breeding_traits\" => string.(collect(hsq_fit.breeding_values.traits)),",
    "\"breeding_values\" => Matrix{Float64}(hsq_fit.breeding_values.values),",
    "\"loglik\" => hsq_fit.loglik,",
    "\"converged\" => hsq_fit.converged,",
    "\"iterations\" => hsq_fit.iterations,",
    "\"traits\" => string.(collect(hsq_fit.traits)),",
    "\"genetic_structure\" => string(hsq_fit.genetic_structure)",
    ");",
    # Number of genetic covariance parameters (contract field for the
    # structure LRT). Read from the engine payload when present; the R
    # normalizer falls back to deriving it from genetic_structure + n_traits.
    "if hasproperty(hsq_fit, :n_genetic_params);",
    "hsq_mv_raw[\"n_genetic_params\"] = hsq_fit.n_genetic_params;",
    "end;",
    # Experimental covariance standard errors (engine row V4-MV-REML, partial;
    # :unstructured only -- the engine throws for structured / factor-analytic
    # fits, and the observed information can be non-positive-definite at a
    # flat/boundary optimum, hence the try guard).
    "if isdefined(HSquared, :multivariate_covariance_standard_errors) &&",
    "hsq_fit.genetic_structure == :unstructured;",
    "hsq_mvse = try; HSquared.multivariate_covariance_standard_errors(",
    "hsq_fit, hsq_Y, hsq_X, hsq_Z, hsq_Ainv);",
    hs_julia_catch_record("multivariate_covariance_standard_errors"),
    "if hsq_mvse !== nothing;",
    "hsq_mv_raw[\"se_genetic_covariance\"] = Matrix{Float64}(hsq_mvse.genetic_covariance);",
    "hsq_mv_raw[\"se_residual_covariance\"] = Matrix{Float64}(hsq_mvse.residual_covariance);",
    "hsq_mv_raw[\"se_genetic_correlation\"] = Matrix{Float64}(hsq_mvse.genetic_correlation);",
    "hsq_mv_raw[\"se_residual_correlation\"] = Matrix{Float64}(hsq_mvse.residual_correlation);",
    "hsq_mv_raw[\"se_heritability\"] = collect(Float64, hsq_mvse.heritability);",
    "end;",
    "end;"
    )),
    hint = hs_dense_scale_hint
  )
  hs_julia_attach_multivariate_plot_data()

  raw <- JuliaCall::julia_eval("hsq_mv_raw")
  result <- hs_normalize_multivariate_result(raw, payload)
  fit <- hs_new_fit(
    spec = list(
      method = "REML",
      family = list(family = payload$family, link = "identity"),
      target = "multivariate"
    ),
    payload = payload,
    result = result,
    engine = "HSquared.jl"
  )
  hs_julia_surface_bridge_errors(fit)
}

hs_normalize_multivariate_result <- function(raw, payload) {
  traits <- as.character(raw$traits %||% payload$metadata$trait_names)
  if (length(traits) == 0L) {
    traits <- paste0("trait", seq_len(ncol(payload$Y)))
  }
  ntraits <- length(traits)
  fixed_names <- payload$metadata$fixed_colnames
  ids <- as.character(raw$breeding_ids %||% payload$ids)

  G0 <- hs_matrix_from_julia(
    raw$genetic_covariance,
    ntraits,
    ntraits,
    "genetic covariance"
  )
  R0 <- hs_matrix_from_julia(
    raw$residual_covariance,
    ntraits,
    ntraits,
    "residual covariance"
  )
  Gcor <- hs_matrix_from_julia(
    raw$genetic_correlation,
    ntraits,
    ntraits,
    "genetic correlation"
  )
  Rcor <- hs_matrix_from_julia(
    raw$residual_correlation,
    ntraits,
    ntraits,
    "residual correlation"
  )
  dimnames(G0) <- dimnames(R0) <- dimnames(Gcor) <- dimnames(Rcor) <-
    list(traits, traits)

  beta <- hs_matrix_from_julia(
    raw$beta,
    length(fixed_names),
    ntraits,
    "fixed effects"
  )
  fixed_effects <- data.frame(
    term = rep(fixed_names, times = ntraits),
    trait = rep(traits, each = length(fixed_names)),
    estimate = as.vector(beta),
    stringsAsFactors = FALSE
  )

  bv <- hs_matrix_from_julia(
    raw$breeding_values,
    length(ids),
    ntraits,
    "breeding values"
  )
  breeding_values <- hs_long_matrix(bv, ids = ids, traits = traits)

  converged <- isTRUE(raw$converged)
  p <- ncol(payload$X)
  n_covariance_parameters <- ntraits * (ntraits + 1L)

  result <- list(
    variance_components = data.frame(
      component = rep(c("genetic", "residual"), each = ntraits),
      trait = rep(traits, times = 2L),
      estimate = c(diag(G0), diag(R0)),
      stringsAsFactors = FALSE
    ),
    heritability = data.frame(
      term = traits,
      trait = traits,
      estimate = as.numeric(raw$heritability),
      stringsAsFactors = FALSE
    ),
    genetic_covariance = G0,
    residual_covariance = R0,
    genetic_correlation = Gcor,
    residual_correlation = Rcor,
    breeding_values = breeding_values,
    fixed_effects = fixed_effects,
    random_effects = list(animal = breeding_values),
    nobs = as.integer(sum(!is.na(payload$Y))),
    converged = converged,
    diagnostics = list(
      target = "multivariate",
      variance_components = "estimated_multivariate_reml",
      optimizer_status = if (converged) "converged" else "not_converged",
      iterations = as.integer(raw$iterations),
      n_traits = ntraits,
      n_records = nrow(payload$Y),
      n_observed_trait_records = sum(!is.na(payload$Y)),
      dense_validation_path = TRUE,
      conditioning_caveat = paste(
        "Experimental dense validation-scale path; the Julia engine inverts",
        "Ainv internally, so deep-inbreeding/high-condition-number pedigrees",
        "remain a twin-side hardening item."
      ),
      genetic_structure = raw$genetic_structure %||% "unstructured"
    )
  )
  if (converged) {
    result$loglik <- as.numeric(raw$loglik)
    result$df <- as.integer(p * ntraits + n_covariance_parameters)
  }
  # Genetic-structure label + number of genetic covariance parameters, for the
  # covariance-structure LRT. Prefer the engine payload field; otherwise derive
  # it from the structure (diagonal = t; unstructured = t(t+1)/2).
  gstruct <- raw$genetic_structure %||% "unstructured"
  result$genetic_structure <- gstruct
  result$n_genetic_params <- if (!is.null(raw$n_genetic_params)) {
    as.integer(raw$n_genetic_params)
  } else if (identical(gstruct, "diagonal")) {
    as.integer(ntraits)
  } else {
    as.integer(ntraits * (ntraits + 1L) / 2L)
  }
  if (!is.null(raw$se_genetic_covariance)) {
    lab <- function(m) {
      m <- as.matrix(m)
      dimnames(m) <- list(traits, traits)
      m
    }
    result$covariance_standard_errors <- list(
      genetic_covariance = lab(raw$se_genetic_covariance),
      residual_covariance = lab(raw$se_residual_covariance),
      genetic_correlation = lab(raw$se_genetic_correlation),
      residual_correlation = lab(raw$se_residual_correlation),
      heritability = stats::setNames(as.numeric(raw$se_heritability), traits)
    )
  }
  result <- hs_attach_multivariate_plot_data(result, raw, traits)
  result
}

hs_validate_multivariate_initial <- function(initial, ntraits) {
  if (is.null(initial)) {
    initial <- list(G0 = diag(1, ntraits), R0 = diag(1, ntraits))
  }
  if (
    !is.list(initial) ||
      is.null(names(initial)) ||
      !all(c("G0", "R0") %in% names(initial))
  ) {
    stop(
      "`initial` for the multivariate target must be a named list with ",
      "`G0` and `R0` covariance matrices.",
      call. = FALSE
    )
  }
  list(
    G0 = hs_validate_initial_covariance(initial$G0, "initial$G0", ntraits),
    R0 = hs_validate_initial_covariance(initial$R0, "initial$R0", ntraits)
  )
}

hs_validate_initial_covariance <- function(x, name, ntraits) {
  x <- as.matrix(x)
  if (!is.numeric(x) || !identical(dim(x), c(ntraits, ntraits))) {
    stop(
      "`",
      name,
      "` must be a numeric ",
      ntraits,
      " x ",
      ntraits,
      " covariance matrix.",
      call. = FALSE
    )
  }
  if (any(!is.finite(x))) {
    stop("`", name, "` must contain only finite values.", call. = FALSE)
  }
  if (!isTRUE(all.equal(x, t(x), tolerance = 1e-8, check.attributes = FALSE))) {
    stop("`", name, "` must be symmetric.", call. = FALSE)
  }
  pd <- tryCatch(
    {
      chol((x + t(x)) / 2)
      TRUE
    },
    error = function(e) FALSE
  )
  if (!isTRUE(pd)) {
    stop("`", name, "` must be positive definite.", call. = FALSE)
  }
  unname(x)
}

hs_matrix_from_julia <- function(x, nrow, ncol, label) {
  out <- as.matrix(x)
  storage.mode(out) <- "double"
  if (!identical(dim(out), c(nrow, ncol))) {
    if (length(out) == nrow * ncol) {
      out <- matrix(as.numeric(out), nrow = nrow, ncol = ncol)
    } else {
      stop(
        "The Julia multivariate result returned ",
        label,
        " with unexpected dimensions.",
        call. = FALSE
      )
    }
  }
  out
}

hs_long_matrix <- function(x, ids, traits) {
  data.frame(
    id = rep(as.character(ids), times = length(traits)),
    trait = rep(as.character(traits), each = length(ids)),
    value = as.vector(x),
    stringsAsFactors = FALSE
  )
}

hs_julia_attach_multivariate_plot_data <- function() {
  JuliaCall::julia_command(paste(
    "if isdefined(HSquared, :genetic_correlation_plot_data);",
    "hsq_gcpd = try;",
    "HSquared.genetic_correlation_plot_data(",
    "hsq_fit.genetic_covariance;",
    "traits = string.(collect(hsq_fit.traits)),",
    "heritabilities = collect(Float64, hsq_fit.heritability));",
    hs_julia_catch_record("genetic_correlation_plot_data"),
    "if hsq_gcpd !== nothing;",
    "hsq_mv_raw[\"genetic_correlation_plot_data\"] = hsq_gcpd;",
    "end;",
    "end;",
    "if isdefined(HSquared, :genetic_pca_plot_data);",
    "hsq_gppd = try;",
    "HSquared.genetic_pca_plot_data(hsq_fit.genetic_covariance);",
    hs_julia_catch_record("genetic_pca_plot_data"),
    "if hsq_gppd !== nothing;",
    "hsq_mv_raw[\"genetic_pca_plot_data\"] = hsq_gppd;",
    "end;",
    "end;"
  ))
  invisible(TRUE)
}

hs_julia_attach_random_regression_plot_data <- function() {
  JuliaCall::julia_command(paste(
    "hsq_rr_plot_ts = collect(range(-1.0, 1.0; length = 25));",
    "if isdefined(HSquared, :rr_genetic_variance_plot_data);",
    "hsq_rr_gvpd = try;",
    "HSquared.rr_genetic_variance_plot_data(",
    "hsq_fit.variance_components.K_g, hsq_rr_plot_ts;",
    "residual = hsq_fit.variance_components.sigma_e2);",
    hs_julia_catch_record("rr_genetic_variance_plot_data"),
    "if hsq_rr_gvpd !== nothing;",
    "hsq_rr_raw[\"rr_genetic_variance_plot_data\"] = hsq_rr_gvpd;",
    "end;",
    "end;",
    "if isdefined(HSquared, :rr_eigenfunctions_plot_data);",
    "hsq_rr_efpd = try;",
    "HSquared.rr_eigenfunctions_plot_data(",
    "hsq_fit.variance_components.K_g, hsq_rr_plot_ts);",
    hs_julia_catch_record("rr_eigenfunctions_plot_data"),
    "if hsq_rr_efpd !== nothing;",
    "hsq_rr_raw[\"rr_eigenfunctions_plot_data\"] = hsq_rr_efpd;",
    "end;",
    "end;",
    "if isdefined(HSquared, :rr_covariance_surface_plot_data);",
    "hsq_rr_sfpd = try;",
    "HSquared.rr_covariance_surface_plot_data(",
    "hsq_fit.variance_components.K_g, hsq_rr_plot_ts);",
    hs_julia_catch_record("rr_covariance_surface_plot_data"),
    "if hsq_rr_sfpd !== nothing;",
    "hsq_rr_raw[\"rr_covariance_surface_plot_data\"] = hsq_rr_sfpd;",
    "end;",
    "end;"
  ))
  invisible(TRUE)
}

# Opt-in, experimental random-regression (reaction-norm) bridge. Mirrors the
# multivariate payload path: assign the univariate response, fixed design, sparse
# record incidence, pedigree, the per-record covariate, and the Legendre order;
# build the n x k Legendre design Phi in Julia from the standardized covariate;
# call the Julia-owned `HSquared.fit_random_regression_reml`; unpack the
# NamedTuple fields into a Dict; normalize to an `hsquared_fit`. The grammar is
# PROVISIONAL (proposed to the twin on HSquared.jl#61, awaiting ack).
hs_fit_julia_random_regression_payload <- function(
  payload,
  project = hs_default_julia_project(),
  iterations = 2000L
) {
  if (!inherits(payload, "hs_bridge_payload")) {
    stop("`payload` must be an internal `hs_bridge_payload`.", call. = FALSE)
  }
  if (is.null(payload$y)) {
    stop(
      "Internal bridge error: the random-regression payload is missing its ",
      "univariate `y` response vector.",
      call. = FALSE
    )
  }
  if (is.null(payload$random_regression)) {
    stop(
      "Internal bridge error: the random-regression target requires an ",
      "`animal(rr(covariate, order = k) | id, ...)` payload.",
      call. = FALSE
    )
  }
  if (is.null(payload$pedigree)) {
    stop(
      "Internal bridge error: the random-regression target requires a ",
      "pedigree animal-model payload.",
      call. = FALSE
    )
  }
  if (!hs_julia_bridge_available(project)) {
    stop(
      "The experimental Julia bridge requires Julia, the `JuliaCall` R ",
      "package, and a local `HSquared.jl` project.",
      call. = FALSE
    )
  }

  rr <- payload$random_regression
  iterations <- hs_validate_iterations(iterations)
  order <- as.integer(rr$order)

  hs_julia_setup(project)
  JuliaCall::julia_assign("hsq_y", payload$y)
  JuliaCall::julia_assign("hsq_X", payload$X)
  hs_julia_assign_sparse_csc("hsq_Z", payload$Z)
  JuliaCall::julia_assign("hsq_id", payload$pedigree$id)
  JuliaCall::julia_assign(
    "hsq_sire",
    hs_parent_for_julia(payload$pedigree$sire)
  )
  JuliaCall::julia_assign("hsq_dam", hs_parent_for_julia(payload$pedigree$dam))
  JuliaCall::julia_assign("hsq_age", as.numeric(rr$values))
  JuliaCall::julia_assign("hsq_order", order)
  JuliaCall::julia_assign("hsq_iterations", iterations)
  hs_julia_fit(
    JuliaCall::julia_command(paste(
      hs_julia_bridge_errors_reset,
      "hsq_ped = HSquared.normalize_pedigree(hsq_id, hsq_sire, hsq_dam);",
      "hsq_Ainv = HSquared.pedigree_inverse(hsq_ped);",
      # Standardize the per-record covariate to [-1, 1] over its observed range and
      # build the n x k normalized-Legendre design Phi (basis convention fixed to
      # Kirkpatrick/Meyer/Schaeffer normalized Legendre on standardized t).
      "hsq_Phi = HSquared.legendre_design(",
      "HSquared.standardize_covariate(hsq_age), hsq_order);",
      "hsq_fit = HSquared.fit_random_regression_reml(",
      "hsq_y, hsq_X, hsq_Phi, hsq_Z, hsq_Ainv;",
      "iterations = hsq_iterations, ids = hsq_ped.ids);",
      "hsq_rr_raw = Dict(",
      "\"K_g\" => Matrix{Float64}(hsq_fit.variance_components.K_g),",
      "\"sigma_e2\" => hsq_fit.variance_components.sigma_e2,",
      "\"beta\" => collect(Float64, hsq_fit.beta),",
      "\"coef_ids\" => string.(collect(hsq_fit.random_coefficients.ids)),",
      "\"coef_values\" => Matrix{Float64}(hsq_fit.random_coefficients.values),",
      "\"loglik\" => hsq_fit.loglik,",
      "\"converged\" => hsq_fit.converged,",
      "\"iterations\" => hsq_fit.iterations,",
      "\"ncoef\" => hsq_fit.basis.ncoef",
      ");"
    )),
    hint = hs_dense_scale_hint
  )
  hs_julia_attach_random_regression_plot_data()

  raw <- JuliaCall::julia_eval("hsq_rr_raw")
  result <- hs_normalize_random_regression_result(raw, payload)
  fit <- hs_new_fit(
    spec = list(
      method = "REML",
      family = list(family = payload$family, link = "identity"),
      target = "random_regression"
    ),
    payload = payload,
    result = result,
    engine = "HSquared.jl"
  )
  hs_julia_surface_bridge_errors(fit)
}

hs_normalize_random_regression_result <- function(raw, payload) {
  rr <- payload$random_regression %||% payload$metadata$random_regression
  fixed_names <- payload$metadata$fixed_colnames
  ids <- as.character(raw$coef_ids %||% payload$ids)
  k <- as.integer(raw$ncoef)
  if (is.na(k) || k < 1L) {
    k <- as.integer(rr$order)
  }
  coef_labels <- paste0("legendre", seq_len(k) - 1L)

  K_g <- hs_matrix_from_julia(raw$K_g, k, k, "coefficient genetic covariance")
  dimnames(K_g) <- list(coef_labels, coef_labels)
  sigma_e2 <- as.numeric(raw$sigma_e2)

  beta <- as.numeric(raw$beta)
  fixed_effects <- data.frame(
    term = fixed_names,
    estimate = beta,
    stringsAsFactors = FALSE
  )

  coefficients <- hs_matrix_from_julia(
    raw$coef_values,
    length(ids),
    k,
    "random-regression coefficients"
  )
  dimnames(coefficients) <- list(ids, coef_labels)
  random_coefficients <- data.frame(
    id = rep(ids, times = k),
    coefficient = rep(coef_labels, each = length(ids)),
    value = as.vector(coefficients),
    stringsAsFactors = FALSE
  )

  converged <- isTRUE(raw$converged)
  p <- ncol(payload$X)
  n_covariance_parameters <- as.integer(k * (k + 1L) / 2L) + 1L

  result <- list(
    variance_components = data.frame(
      component = c(paste0("K_g_", coef_labels), "residual"),
      estimate = c(diag(K_g), sigma_e2),
      stringsAsFactors = FALSE
    ),
    coefficient_covariance = K_g,
    residual_variance = sigma_e2,
    random_coefficients = random_coefficients,
    fixed_effects = fixed_effects,
    nobs = as.integer(length(payload$y)),
    converged = converged,
    # Standardization + basis metadata: extractors recompute the genetic
    # variance / heritability / correlation trajectories in R from K_g and the
    # normalized-Legendre basis, re-standardizing any user-supplied `at =` on the
    # original covariate scale with these recorded bounds.
    random_regression = list(
      covariate = rr$covariate,
      order = k,
      lower = as.numeric(rr$lower),
      upper = as.numeric(rr$upper)
    ),
    diagnostics = list(
      target = "random_regression",
      variance_components = "estimated_random_regression_reml",
      optimizer_status = if (converged) "converged" else "not_converged",
      iterations = as.integer(raw$iterations),
      n_coefficients = k,
      n_records = length(payload$y),
      covariate = rr$covariate,
      covariate_range = c(as.numeric(rr$lower), as.numeric(rr$upper)),
      dense_validation_path = TRUE,
      residual_model = "homogeneous",
      conditioning_caveat = paste(
        "Experimental dense validation-scale path with a HOMOGENEOUS residual",
        "and NO permanent-environment term; both are planned. The Julia engine",
        "inverts Ainv internally, so deep-inbreeding/high-condition-number",
        "pedigrees remain a twin-side hardening item."
      )
    )
  )
  if (converged) {
    result$loglik <- as.numeric(raw$loglik)
    result$df <- as.integer(p + n_covariance_parameters)
  }
  result <- hs_attach_random_regression_plot_data(result, raw)
  result
}

# Normalized Legendre basis row vector phi(t) = [phi_0(t), ..., phi_{k-1}(t)] at a
# standardized covariate t in [-1, 1], mirroring `HSquared.legendre_basis`:
# phi_n(t) = sqrt((2n+1)/2) * P_n(t), with P_n the ordinary Legendre polynomials
# via the Bonnet recurrence. Used by the R-side reaction-norm trajectory
# extractors so they need no live Julia round-trip. Errors -- does not
# clamp -- for |t| > 1 + 1e-10, matching the Julia engine's tolerance exactly
# (#213); within that tolerance t is clamped to [-1, 1] before use, as Julia
# also does.
hs_legendre_basis <- function(t, order) {
  order <- as.integer(order)
  t <- as.numeric(t)
  if (t < -1 - 1e-10 || t > 1 + 1e-10) {
    hs_abort_out_of_range(
      "`t` must be in [-1, 1]; standardize the covariate first ",
      "(see `hs_standardize_covariate()`), got t = ",
      format(t, trim = TRUE),
      "."
    )
  }
  tt <- max(-1, min(1, t))
  p <- numeric(order)
  p[1L] <- 1
  if (order >= 2L) {
    p[2L] <- tt
  }
  if (order >= 3L) {
    for (n in 2:(order - 1L)) {
      p[n + 1L] <- ((2 * n - 1) * tt * p[n] - (n - 1) * p[n - 1L]) / n
    }
  }
  phi <- numeric(order)
  for (n in 0:(order - 1L)) {
    phi[n + 1L] <- sqrt((2 * n + 1) / 2) * p[n + 1L]
  }
  phi
}

# n x order normalized-Legendre design over already-standardized points `ts`,
# mirroring `HSquared.legendre_design`. Errors via `hs_legendre_basis()` if any
# `ts` point is outside [-1 - 1e-10, 1 + 1e-10].
hs_legendre_design <- function(ts, order) {
  do.call(rbind, lapply(ts, hs_legendre_basis, order = order))
}

# Map a raw covariate value/vector onto t in [-1, 1] using the recorded fit
# range, mirroring `HSquared.standardize_covariate`
# (t = 2(a - lower)/(upper - lower) - 1).
hs_standardize_covariate <- function(a, lower, upper) {
  2 * (as.numeric(a) - lower) / (upper - lower) - 1
}

hs_unstandardize_covariate <- function(t, lower, upper) {
  lower + (as.numeric(t) + 1) * (upper - lower) / 2
}

hs_validate_two_effect_initial <- function(initial) {
  if (
    !is.numeric(initial) ||
      !setequal(names(initial), c("sigma_a2", "sigma_c2", "sigma_e2"))
  ) {
    stop(
      "`initial` for the two-effect target must be a named numeric vector ",
      "with `sigma_a2`, `sigma_c2`, and `sigma_e2`.",
      call. = FALSE
    )
  }
  if (any(!is.finite(initial)) || any(initial <= 0)) {
    stop(
      "`initial` variance components must be finite and positive.",
      call. = FALSE
    )
  }
  initial[c("sigma_a2", "sigma_c2", "sigma_e2")]
}

# Opt-in, experimental genomic GREML estimator. Surfaces the Julia-owned
# `HSquared.fit_ai_reml()` REML optimizer on an animal_model_spec built with the
# user-supplied genomic relationship inverse `Ginv` (in place of a pedigree
# Ainv). Reuses the standard result normalizer, relabelling the genetic
# component as "genomic".
hs_fit_julia_genomic_payload <- function(
  payload,
  project = hs_default_julia_project(),
  initial = c(sigma_a2 = 1, sigma_e2 = 1),
  iterations = 100L
) {
  if (!inherits(payload, "hs_bridge_payload")) {
    stop("`payload` must be an internal `hs_bridge_payload`.", call. = FALSE)
  }
  from_markers <- identical(payload$relationship_source, "markers")
  if (is.null(payload$Ginv) && !from_markers) {
    stop(
      "Internal bridge error: the genomic payload is missing its `Ginv`.",
      call. = FALSE
    )
  }
  if (from_markers && is.null(payload$markers)) {
    stop(
      "Internal bridge error: the genomic payload is missing its `markers`.",
      call. = FALSE
    )
  }
  if (!hs_julia_bridge_available(project)) {
    stop(
      "The experimental Julia bridge requires Julia, the `JuliaCall` R ",
      "package, and a local `HSquared.jl` project.",
      call. = FALSE
    )
  }

  initial <- hs_validate_initial_variances(initial)
  iterations <- hs_validate_iterations(iterations)
  hs_julia_setup(project)
  JuliaCall::julia_assign("hsq_y", payload$y)
  JuliaCall::julia_assign("hsq_X", payload$X)
  hs_julia_assign_sparse_csc("hsq_Z", payload$Z)
  JuliaCall::julia_assign("hsq_ids", payload$ids)
  JuliaCall::julia_assign("hsq_initial_sigma_a2", unname(initial[["sigma_a2"]]))
  JuliaCall::julia_assign("hsq_initial_sigma_e2", unname(initial[["sigma_e2"]]))
  JuliaCall::julia_assign("hsq_iterations", iterations)
  if (from_markers) {
    # Build the frozen sample-p VanRaden-1 relationship and its regularized
    # precision in Julia. The engine owns the canonical construction metadata
    # and SHA-256 provenance fingerprints; R carries them unchanged.
    JuliaCall::julia_assign("hsq_markers", payload$markers)
    JuliaCall::julia_assign("hsq_ridge", payload$ridge)
    if (is.null(payload$marker_names)) {
      marker_names_cmd <- "hsq_marker_names = nothing;"
    } else {
      JuliaCall::julia_assign("hsq_marker_names", payload$marker_names)
      marker_names_cmd <- ""
    }
    relinv_cmd <- paste(
      marker_names_cmd,
      "hsq_genomic_construction =",
      "HSquared._genomic_activation_construction(",
      "hsq_markers, hsq_ids; marker_names = hsq_marker_names,",
      "ridge = hsq_ridge);",
      "hsq_Ginvs = sparse(hsq_genomic_construction.Q);",
      "hsq_genomic_provenance = hsq_genomic_construction.provenance;",
      "hsq_genomic_kernel = hsq_genomic_construction.K;"
    )
  } else {
    JuliaCall::julia_assign("hsq_Ginv", payload$Ginv)
    relinv_cmd <- if (identical(payload$relationship, "genomic")) {
      paste(
        "hsq_Ginvs = sparse(hsq_Ginv);",
        "hsq_genomic_provenance =",
        "HSquared._genomic_precision_provenance(hsq_Ginv, hsq_ids);",
        "hsq_genomic_kernel = nothing;"
      )
    } else {
      paste(
        "hsq_Ginvs = sparse(hsq_Ginv);",
        "hsq_genomic_provenance = nothing;",
        "hsq_genomic_kernel = nothing;"
      )
    }
  }
  rel <- payload$relationship
  boundary_eligible <- identical(rel, "genomic") &&
    nrow(payload$Z) == ncol(payload$Z) &&
    nrow(payload$Z) <= 2000L &&
    isTRUE(all.equal(
      as.matrix(payload$Z),
      diag(nrow(payload$Z)),
      tolerance = 1e-12,
      check.attributes = FALSE
    ))
  fit_cmd <- if (boundary_eligible) {
    paste(
      "hsq_boundary_result = HSquared._fit_ai_reml_genomic_boundary(",
      "hsq_spec;",
      "provenance = hsq_genomic_provenance,",
      "kernel = hsq_genomic_kernel,",
      "initial = (sigma_a2 = hsq_initial_sigma_a2,",
      "sigma_e2 = hsq_initial_sigma_e2),",
      "iterations = hsq_iterations);",
      "hsq_fit = hsq_boundary_result.fit;"
    )
  } else {
    paste(
      "hsq_boundary_result = nothing;",
      "hsq_fit = HSquared.fit_ai_reml(",
      "hsq_spec;",
      "initial = (sigma_a2 = hsq_initial_sigma_a2,",
      "sigma_e2 = hsq_initial_sigma_e2),",
      "iterations = hsq_iterations);"
    )
  }
  hs_julia_fit(
    JuliaCall::julia_command(paste(
      hs_julia_bridge_errors_reset,
      relinv_cmd,
      "hsq_spec = HSquared.animal_model_spec(",
      "hsq_y, hsq_X, hsq_Z, hsq_Ginvs;",
      "ids = hsq_ids, method = :REML);",
      fit_cmd,
      "hsq_result = HSquared.result_payload(hsq_fit);"
    )),
    hint = hs_dense_scale_hint
  )
  hs_julia_attach_standard_plot_data()

  raw <- JuliaCall::julia_eval(
    "Dict(String(k) => getfield(hsq_result, k) for k in keys(hsq_result))"
  )
  raw_provenance <- if (identical(payload$relationship, "genomic")) {
    JuliaCall::julia_eval(paste0(
      "Dict(String(k) => (getfield(hsq_genomic_provenance, k) === nothing ",
      "? missing : getfield(hsq_genomic_provenance, k)) ",
      "for k in keys(hsq_genomic_provenance))"
    ))
  } else {
    NULL
  }
  raw_boundary <- if (boundary_eligible) {
    JuliaCall::julia_eval(paste0(
      "Dict(String(k) => (getfield(hsq_boundary_result.boundary, k) === nothing ",
      "? missing : getfield(hsq_boundary_result.boundary, k)) ",
      "for k in keys(hsq_boundary_result.boundary))"
    ))
  } else {
    NULL
  }
  result <- hs_normalize_julia_result(raw, payload)
  result$variance_components$component[
    result$variance_components$component == "animal"
  ] <- rel
  result$heritability$term[result$heritability$term == "animal"] <- rel
  names(result$random_effects)[
    names(result$random_effects) == "animal"
  ] <- rel
  result$diagnostics$variance_components <- paste0(
    "estimated_",
    rel,
    "_ai_reml"
  )
  if (identical(rel, "genomic")) {
    provenance <- hs_normalize_genomic_provenance(raw_provenance)
    payload$relationship_provenance <- provenance
    result$relationship_provenance <- provenance
    result$heritability$component <- "genomic_variance_ratio"
    result$heritability$relationship_scale <- provenance$relationship_scale
    result$heritability$relationship_source <- provenance$relationship_source
    result$heritability$relationship_method <- provenance$relationship_method
    result$heritability$allele_frequency_source <-
      provenance$allele_frequency_source
    result$heritability$ridge <- provenance$ridge
    if (!is.null(raw_boundary)) {
      boundary <- hs_normalize_genomic_boundary(raw_boundary)
      result$genomic_boundary <- boundary
      result$heritability$numerical_estimate <- boundary$numerical_ratio
      if (!is.na(boundary$profile_ratio)) {
        result$heritability$estimate <- boundary$profile_ratio
      }
    }
    # Genomic ratio uncertainty is not yet scale-labelled or separately
    # calibrated. Keep the engine's raw capability out of the public R result
    # until that contract is validated.
    result$heritability_interval <- NULL
    result$heritability_se <- NULL
    if (!is.null(result$genomic_boundary) &&
        result$genomic_boundary$status %in% c("boundary_lower", "boundary_upper")) {
      result$breeding_values <- NULL
      result$breeding_values_plot_data <- NULL
      result$random_effects <- NULL
      result$predictions <- NULL
      result$prediction_error_variance <- NULL
      result$reliability <- NULL
      result$variance_component_se <- NULL
    }
    if (!is.null(result$variance_component_se)) {
      result$variance_component_se$component[
        result$variance_component_se$component == "animal"
      ] <- "genomic"
    }
  }
  fit <- hs_new_fit(
    spec = list(
      method = "REML",
      family = list(family = payload$family, link = "identity"),
      target = rel
    ),
    payload = payload,
    result = result,
    engine = "HSquared.jl"
  )
  hs_julia_surface_bridge_errors(fit)
}

hs_normalize_genomic_boundary <- function(raw) {
  raw <- hs_drop_julia_classes(raw)
  if (is.data.frame(raw)) raw <- as.list(raw)
  if (!is.list(raw)) {
    stop("Internal bridge error: genomic boundary metadata is missing.", call. = FALSE)
  }
  scalar_character <- function(name) {
    x <- raw[[name]]
    if (is.null(x) || !length(x) || all(is.na(x))) NA_character_ else as.character(x[[1L]])
  }
  scalar_numeric <- function(name) {
    x <- raw[[name]]
    if (is.null(x) || !length(x) || all(is.na(x))) NA_real_ else as.numeric(x[[1L]])
  }
  out <- list(
    status = scalar_character("status"),
    reason = scalar_character("reason"),
    profile_ratio = scalar_numeric("profile_ratio"),
    numerical_ratio = scalar_numeric("numerical_ratio"),
    boundary_epsilon = scalar_numeric("boundary_epsilon"),
    profile_loglik = scalar_numeric("profile_loglik"),
    lower_derivative_per_observation = scalar_numeric("lower_derivative_per_observation"),
    upper_derivative_per_observation = scalar_numeric("upper_derivative_per_observation")
  )
  allowed <- c("boundary_lower", "boundary_upper", "interior", "interior_rescued", "boundary_unresolved")
  if (is.na(out$status) || !out$status %in% allowed) {
    stop("Internal bridge error: unknown genomic boundary status.", call. = FALSE)
  }
  if (!identical(out$boundary_epsilon, 1e-7)) {
    stop("Internal bridge error: genomic boundary epsilon drift.", call. = FALSE)
  }
  resolved <- !identical(out$status, "boundary_unresolved")
  required <- unlist(out[c("profile_ratio", "numerical_ratio", "profile_loglik",
    "lower_derivative_per_observation", "upper_derivative_per_observation")])
  if (resolved && any(!is.finite(required))) {
    stop("Internal bridge error: resolved genomic boundary metadata is non-finite.", call. = FALSE)
  }
  if (identical(out$status, "boundary_lower") &&
      (!identical(out$profile_ratio, 0) || !identical(out$numerical_ratio, 1e-7))) {
    stop("Internal bridge error: lower-boundary ratio contract drift.", call. = FALSE)
  }
  if (identical(out$status, "boundary_upper") &&
      (!identical(out$profile_ratio, 1) || !identical(out$numerical_ratio, 1 - 1e-7))) {
    stop("Internal bridge error: upper-boundary ratio contract drift.", call. = FALSE)
  }
  out
}

hs_v07_genomic_boundary_contract <- function() {
  list(
    doc46_commit = "fe96a147",
    doc46_sha256 = "283ab00bab3da925f0ac2916959efacaa7fb711c5da4dce09dd49ea568eef030",
    julia_implementation_commit = "ecc058f380be71058c9cfde373c345ab7a2f6aba",
    boundary_epsilon = 1e-7,
    grid_step = 0.0025,
    derivative_delta = 1e-6,
    kkt_tolerance = 1e-8,
    candidate_id = "v07_genomic_closed_boundary_v1"
  )
}

hs_normalize_genomic_provenance <- function(raw) {
  raw <- hs_drop_julia_classes(raw)
  if (is.data.frame(raw)) {
    raw <- as.list(raw)
  }
  if (!is.list(raw)) {
    stop(
      "Internal bridge error: genomic relationship provenance is missing.",
      call. = FALSE
    )
  }

  value <- function(name, type = c("character", "numeric")) {
    type <- match.arg(type)
    x <- raw[[name]]
    if (is.null(x) || length(x) < 1L || all(is.na(x))) {
      return(if (identical(type, "numeric")) NA_real_ else NA_character_)
    }
    if (identical(type, "numeric")) {
      as.numeric(x[[1L]])
    } else {
      as.character(x[[1L]])
    }
  }

  out <- list(
    relationship_source = value("relationship_source"),
    relationship_method = value("relationship_method"),
    allele_frequency_source = value("allele_frequency_source"),
    ridge = value("ridge", "numeric"),
    scale_denominator = value("scale_denominator", "numeric"),
    relationship_scale = value("relationship_scale"),
    id_order_fingerprint = value("id_order_fingerprint"),
    marker_content_fingerprint = value("marker_content_fingerprint"),
    kernel_fingerprint = value("kernel_fingerprint"),
    precision_fingerprint = value("precision_fingerprint")
  )

  if (
    length(out$relationship_source) != 1L ||
      is.na(out$relationship_source) ||
      !out$relationship_source %in% c("markers", "supplied_Ginv")
  ) {
    stop(
      "Internal bridge error: genomic provenance source must be exactly ",
      "`markers` or `supplied_Ginv`.",
      call. = FALSE
    )
  }

  required <- c("id_order_fingerprint", "precision_fingerprint")
  if (identical(out$relationship_source, "markers")) {
    required <- c(required, "marker_content_fingerprint", "kernel_fingerprint")
  }
  valid_sha256 <- function(x) {
    length(x) == 1L && !is.na(x) && grepl("^[0-9a-f]{64}$", x)
  }
  if (!all(vapply(out[required], valid_sha256, logical(1L)))) {
    stop(
      "Internal bridge error: genomic provenance fingerprints must be ",
      "lowercase SHA-256 values.",
      call. = FALSE
    )
  }
  if (
    identical(out$relationship_source, "markers") &&
      (is.na(out$relationship_method) ||
        is.na(out$allele_frequency_source) ||
        !is.finite(out$ridge) ||
        !is.finite(out$scale_denominator) ||
        !identical(out$relationship_method, "vanraden1") ||
        !identical(out$allele_frequency_source, "sample") ||
        !isTRUE(all.equal(out$ridge, 0.01)) ||
        !is.finite(out$scale_denominator) ||
        out$scale_denominator <= 0 ||
        !identical(out$relationship_scale, "K_lambda"))
  ) {
    stop(
      "Internal bridge error: marker provenance does not match the frozen ",
      "VanRaden-1/sample-p/ridge-0.01 contract.",
      call. = FALSE
    )
  }
  if (
    identical(out$relationship_source, "supplied_Ginv") &&
      (!is.na(out$relationship_method) ||
        !is.na(out$allele_frequency_source) ||
        !is.na(out$ridge) ||
        !is.na(out$scale_denominator) ||
        !is.na(out$marker_content_fingerprint) ||
        !is.na(out$kernel_fingerprint) ||
        !identical(
          out$relationship_scale,
          "inverse_of_supplied_precision"
        ))
  ) {
    stop(
      "Internal bridge error: supplied-Ginv provenance must leave its ",
      "construction method, allele frequencies, ridge, denominator, and ",
      "kernel unknown.",
      call. = FALSE
    )
  }
  out
}

# Opt-in, experimental single-step H^-1 CONSTRUCTION estimator. Builds Ainv + dense
# A from the pedigree and G from the genotyped-subset markers, then fits via the
# Julia-owned `fit_single_step_reml` (which assembles H^-1 = A^-1 + scatter over the
# genotyped rows). Mirrors `fit_ai_reml` on the supplied-Hinv path; experimental,
# dense/validation-scale (docs/design/25).
hs_fit_julia_single_step_construct_payload <- function(
  payload,
  project = hs_default_julia_project(),
  initial = c(sigma_a2 = 1, sigma_e2 = 1),
  iterations = 100L
) {
  if (!inherits(payload, "hs_bridge_payload")) {
    stop("`payload` must be an internal `hs_bridge_payload`.", call. = FALSE)
  }
  if (
    !identical(payload$relationship_source, "construct") ||
      is.null(payload$markers) ||
      is.null(payload$pedigree) ||
      is.null(payload$genotyped_rows)
  ) {
    stop(
      "Internal bridge error: the single-step construction payload is ",
      "incomplete (needs pedigree, markers, and genotyped_rows).",
      call. = FALSE
    )
  }
  if (!hs_julia_bridge_available(project)) {
    stop(
      "The experimental Julia bridge requires Julia, the `JuliaCall` R ",
      "package, and a local `HSquared.jl` project.",
      call. = FALSE
    )
  }

  initial <- hs_validate_initial_variances(initial)
  iterations <- hs_validate_iterations(iterations)
  hs_julia_setup(project)
  JuliaCall::julia_assign("hsq_y", payload$y)
  JuliaCall::julia_assign("hsq_X", payload$X)
  hs_julia_assign_sparse_csc("hsq_Z", payload$Z)
  JuliaCall::julia_assign("hsq_id", payload$pedigree$id)
  JuliaCall::julia_assign(
    "hsq_sire",
    hs_parent_for_julia(payload$pedigree$sire)
  )
  JuliaCall::julia_assign("hsq_dam", hs_parent_for_julia(payload$pedigree$dam))
  JuliaCall::julia_assign("hsq_markers", payload$markers)
  JuliaCall::julia_assign(
    "hsq_grows",
    as.integer(payload$genotyped_rows)
  )
  JuliaCall::julia_assign("hsq_tau", payload$single_step$tau)
  JuliaCall::julia_assign("hsq_omega", payload$single_step$omega)
  JuliaCall::julia_assign("hsq_bw", payload$single_step$blend_weight)
  JuliaCall::julia_assign("hsq_ssridge", payload$single_step$ridge)
  JuliaCall::julia_assign("hsq_initial_sigma_a2", unname(initial[["sigma_a2"]]))
  JuliaCall::julia_assign("hsq_initial_sigma_e2", unname(initial[["sigma_e2"]]))
  JuliaCall::julia_assign("hsq_iterations", iterations)
  hs_julia_fit(
    JuliaCall::julia_command(paste(
      hs_julia_bridge_errors_reset,
      "hsq_ped = HSquared.normalize_pedigree(hsq_id, hsq_sire, hsq_dam);",
      # Guard the genotyped_rows alignment (docs/design/25 section 8): the R-computed
      # genotyped_rows index R's pedigree order, so the engine's normalize_pedigree
      # must preserve that order. Fail loudly rather than fit a misaligned G.
      "collect(String, hsq_ped.ids) == hsq_id ||",
      "error(\"single_step construct: engine pedigree order != R order\");",
      "hsq_Ainv = HSquared.pedigree_inverse(hsq_ped);",
      "hsq_A = HSquared.additive_relationship(hsq_ped);",
      "hsq_G = HSquared.genomic_relationship_matrix(hsq_markers);",
      "hsq_fit = HSquared.fit_single_step_reml(",
      "hsq_y, hsq_X, hsq_Z, hsq_Ainv, hsq_A, hsq_G, hsq_grows;",
      "ids = hsq_ped.ids,",
      "tau = hsq_tau, omega = hsq_omega, blend_weight = hsq_bw, ridge = hsq_ssridge,",
      "initial = (sigma_a2 = hsq_initial_sigma_a2,",
      # hsquared#212: `iterations` was validated and assigned above but never
      # reached the fit call -- silently discarded. Forward it now.
      "sigma_e2 = hsq_initial_sigma_e2), iterations = hsq_iterations);",
      "hsq_result = HSquared.result_payload(hsq_fit);"
    )),
    hint = hs_single_step_ridge_hint
  )
  hs_julia_attach_standard_plot_data()

  raw <- JuliaCall::julia_eval(
    "Dict(String(k) => getfield(hsq_result, k) for k in keys(hsq_result))"
  )
  rel <- payload$relationship
  result <- hs_normalize_julia_result(raw, payload)
  result$variance_components$component[
    result$variance_components$component == "animal"
  ] <- rel
  result$heritability$term[result$heritability$term == "animal"] <- rel
  names(result$random_effects)[
    names(result$random_effects) == "animal"
  ] <- rel
  result$diagnostics$variance_components <- "estimated_single_step_construct_ai_reml"
  fit <- hs_new_fit(
    spec = list(
      method = "REML",
      family = list(family = payload$family, link = "identity"),
      target = "single_step_construct"
    ),
    payload = payload,
    result = result,
    engine = "HSquared.jl"
  )
  hs_julia_surface_bridge_errors(fit)
}

# Opt-in, experimental supplied-Gamma H^Gamma single-step bridge. This mirrors
# the ordinary construction helper above, but delegates H construction to the
# Julia-owned metafounder precision path. Gamma is supplied by the user; this
# bridge does not estimate Gamma or expose metafounder-specific extractors.
hs_fit_julia_metafounder_single_step_payload <- function(
  payload,
  project = hs_default_julia_project(),
  initial = c(sigma_a2 = 1, sigma_e2 = 1),
  iterations = 100L
) {
  if (!inherits(payload, "hs_bridge_payload")) {
    stop("`payload` must be an internal `hs_bridge_payload`.", call. = FALSE)
  }
  if (
    !identical(payload$relationship_source, "metafounder_single_step") ||
      is.null(payload$markers) ||
      is.null(payload$pedigree) ||
      is.null(payload$genotyped_rows) ||
      is.null(payload$group_of) ||
      is.null(payload$Gamma)
  ) {
    stop(
      "Internal bridge error: the metafounder single-step payload is ",
      "incomplete (needs pedigree, markers, genotyped_rows, group_of, and ",
      "Gamma).",
      call. = FALSE
    )
  }
  if (!hs_julia_bridge_available(project)) {
    stop(
      "The experimental Julia bridge requires Julia, the `JuliaCall` R ",
      "package, and a local `HSquared.jl` project.",
      call. = FALSE
    )
  }

  initial <- hs_validate_initial_variances(initial)
  iterations <- hs_validate_iterations(iterations)
  hs_julia_setup(project)
  JuliaCall::julia_assign("hsq_y", payload$y)
  JuliaCall::julia_assign("hsq_X", payload$X)
  hs_julia_assign_sparse_csc("hsq_Z", payload$Z)
  JuliaCall::julia_assign("hsq_id", payload$pedigree$id)
  JuliaCall::julia_assign(
    "hsq_sire",
    hs_parent_for_julia(payload$pedigree$sire)
  )
  JuliaCall::julia_assign("hsq_dam", hs_parent_for_julia(payload$pedigree$dam))
  JuliaCall::julia_assign("hsq_markers", payload$markers)
  JuliaCall::julia_assign(
    "hsq_grows",
    as.integer(payload$genotyped_rows)
  )
  JuliaCall::julia_assign("hsq_group_of", unname(payload$group_of))
  JuliaCall::julia_assign("hsq_Gamma_vec", as.numeric(payload$Gamma))
  JuliaCall::julia_assign("hsq_Gamma_n", as.integer(nrow(payload$Gamma)))
  JuliaCall::julia_assign("hsq_tau", payload$single_step$tau)
  JuliaCall::julia_assign("hsq_omega", payload$single_step$omega)
  JuliaCall::julia_assign("hsq_bw", payload$single_step$blend_weight)
  JuliaCall::julia_assign("hsq_ssridge", payload$single_step$ridge)
  JuliaCall::julia_assign("hsq_initial_sigma_a2", unname(initial[["sigma_a2"]]))
  JuliaCall::julia_assign("hsq_initial_sigma_e2", unname(initial[["sigma_e2"]]))
  JuliaCall::julia_assign("hsq_iterations", iterations)
  hs_julia_fit(
    JuliaCall::julia_command(paste(
      hs_julia_bridge_errors_reset,
      "hsq_ped = HSquared.normalize_pedigree(hsq_id, hsq_sire, hsq_dam);",
      "collect(String, hsq_ped.ids) == hsq_id ||",
      "error(\"metafounder single_step: engine pedigree order != R order\");",
      "hsq_G = HSquared.genomic_relationship_matrix(hsq_markers);",
      "hsq_Gamma = reshape(collect(Float64, hsq_Gamma_vec),",
      "Int(hsq_Gamma_n), Int(hsq_Gamma_n));",
      "hsq_fit = HSquared.fit_metafounder_single_step_reml(",
      "hsq_y, hsq_X, hsq_Z, hsq_ped, hsq_group_of, hsq_Gamma, hsq_G, hsq_grows;",
      "ids = hsq_ped.ids,",
      "tau = hsq_tau, omega = hsq_omega, blend_weight = hsq_bw, ridge = hsq_ssridge,",
      "initial = (sigma_a2 = hsq_initial_sigma_a2,",
      # hsquared#212: `iterations` was validated and assigned above but never
      # reached the fit call -- silently discarded. Forward it now.
      "sigma_e2 = hsq_initial_sigma_e2), iterations = hsq_iterations);",
      "hsq_result = HSquared.result_payload(hsq_fit);"
    )),
    hint = hs_single_step_ridge_hint
  )
  hs_julia_attach_standard_plot_data()

  raw <- JuliaCall::julia_eval(
    "Dict(String(k) => getfield(hsq_result, k) for k in keys(hsq_result))"
  )
  rel <- payload$relationship
  result <- hs_normalize_julia_result(raw, payload)
  result$variance_components$component[
    result$variance_components$component == "animal"
  ] <- rel
  result$heritability$term[result$heritability$term == "animal"] <- rel
  names(result$random_effects)[
    names(result$random_effects) == "animal"
  ] <- rel
  result$diagnostics$target <- "metafounder_single_step"
  result$diagnostics$variance_components <-
    "estimated_metafounder_single_step_ai_reml"
  result$diagnostics$gamma_source <- "supplied"
  fit <- hs_new_fit(
    spec = list(
      method = "REML",
      family = list(family = payload$family, link = "identity"),
      target = "metafounder_single_step"
    ),
    payload = payload,
    result = result,
    engine = "HSquared.jl"
  )
  hs_julia_surface_bridge_errors(fit)
}

hs_fit_julia_snp_blup_payload <- function(
  payload,
  project = hs_default_julia_project(),
  variance_components = NULL
) {
  if (!inherits(payload, "hs_bridge_payload")) {
    stop("`payload` must be an internal `hs_bridge_payload`.", call. = FALSE)
  }
  if (
    !identical(payload$relationship_source, "markers") ||
      is.null(payload$markers)
  ) {
    stop(
      "Internal bridge error: SNP-BLUP requires a marker-matrix genomic ",
      "payload.",
      call. = FALSE
    )
  }
  if (!hs_julia_bridge_available(project)) {
    stop(
      "The experimental Julia bridge requires Julia, the `JuliaCall` R ",
      "package, and a local `HSquared.jl` project.",
      call. = FALSE
    )
  }

  variance_components <- hs_validate_snp_blup_variances(variance_components)
  sigma_g2 <- unname(variance_components[["sigma_g2"]])
  sigma_e2 <- unname(variance_components[["sigma_e2"]])

  markers_ind <- payload$markers
  # Per-record marker design: each record carries its individual's genotype.
  # `Z` maps records to the `ids` order that the `markers` rows are in, so
  # `Z %*% markers` aligns the marker rows with the response `y`.
  markers_rec <- as.matrix(payload$Z %*% markers_ind)

  hs_julia_setup(project)
  JuliaCall::julia_assign("hsq_y", payload$y)
  JuliaCall::julia_assign("hsq_X", payload$X)
  JuliaCall::julia_assign("hsq_markers_rec", markers_rec)
  JuliaCall::julia_assign("hsq_markers_ind", unname(markers_ind))
  JuliaCall::julia_assign("hsq_sigma_g2", sigma_g2)
  JuliaCall::julia_assign("hsq_sigma_e2", sigma_e2)
  hs_julia_fit(
    JuliaCall::julia_command(paste(
      "hsq_snp = HSquared.fit_snp_blup(",
      "hsq_y, hsq_X, hsq_markers_rec, hsq_sigma_g2, hsq_sigma_e2);",
      # Per-individual GEBV at the same allele-frequency centering as the fit.
      "hsq_Wind = HSquared.centered_markers(",
      "hsq_markers_ind; allele_frequencies = hsq_snp.p).W;",
      "hsq_gebv_ind = hsq_Wind * hsq_snp.marker_effects;",
      "hsq_fitted = hsq_X * hsq_snp.beta .+ hsq_snp.gebv;",
      "hsq_snp_raw = Dict(",
      "\"marker_effects\" => hsq_snp.marker_effects,",
      "\"gebv\" => hsq_gebv_ind,",
      "\"beta\" => hsq_snp.beta,",
      "\"p\" => hsq_snp.p,",
      "\"fitted\" => hsq_fitted,",
      "\"k\" => hsq_snp.k,",
      "\"nobs\" => length(hsq_y)",
      ");"
    )),
    hint = hs_dense_scale_hint
  )

  raw <- JuliaCall::julia_eval("hsq_snp_raw")
  result <- hs_normalize_julia_snp_blup_result(
    raw,
    payload,
    variance_components
  )
  hs_new_fit(
    spec = list(
      method = payload$method,
      family = list(family = payload$family, link = "identity"),
      target = "snp_blup"
    ),
    payload = payload,
    result = result,
    engine = "HSquared.jl"
  )
}

# Opt-in, experimental REML-estimated SNP-BLUP. Unlike the supplied-variance
# `hs_fit_julia_snp_blup_payload`, this estimates the genomic and residual
# variance components from the markers by REML (Julia-owned
# `HSquared.fit_snp_blup_reml`), so `genomic(1 | id, markers = M)` no longer needs
# the user to supply `sigma_g2`/`sigma_e2`. Experimental, dense/validation-scale
# (mirrors the engine row V2-SNPBLUP, partial).
hs_fit_julia_snp_blup_reml_payload <- function(
  payload,
  project = hs_default_julia_project()
) {
  if (!inherits(payload, "hs_bridge_payload")) {
    stop("`payload` must be an internal `hs_bridge_payload`.", call. = FALSE)
  }
  if (
    !identical(payload$relationship_source, "markers") ||
      is.null(payload$markers)
  ) {
    stop(
      "Internal bridge error: SNP-BLUP requires a marker-matrix genomic payload.",
      call. = FALSE
    )
  }
  if (!hs_julia_bridge_available(project)) {
    stop(
      "The experimental Julia bridge requires Julia, the `JuliaCall` R ",
      "package, and a local `HSquared.jl` project.",
      call. = FALSE
    )
  }

  markers_ind <- payload$markers
  markers_rec <- as.matrix(payload$Z %*% markers_ind)

  hs_julia_setup(project)
  JuliaCall::julia_assign("hsq_y", payload$y)
  JuliaCall::julia_assign("hsq_X", payload$X)
  JuliaCall::julia_assign("hsq_markers_rec", markers_rec)
  JuliaCall::julia_assign("hsq_markers_ind", unname(markers_ind))
  hs_julia_fit(
    JuliaCall::julia_command(paste(
      "hsq_snp = HSquared.fit_snp_blup_reml(hsq_y, hsq_X, hsq_markers_rec);",
      # Per-individual GEBV at the same allele-frequency centering as the fit.
      "hsq_Wind = HSquared.centered_markers(",
      "hsq_markers_ind; allele_frequencies = hsq_snp.p).W;",
      "hsq_gebv_ind = hsq_Wind * hsq_snp.marker_effects;",
      "hsq_fitted = hsq_X * hsq_snp.beta .+ hsq_snp.gebv;",
      "hsq_snp_raw = Dict(",
      "\"marker_effects\" => hsq_snp.marker_effects,",
      "\"gebv\" => hsq_gebv_ind,",
      "\"beta\" => hsq_snp.beta,",
      "\"p\" => hsq_snp.p,",
      "\"fitted\" => hsq_fitted,",
      "\"k\" => hsq_snp.k,",
      "\"sigma_g2\" => hsq_snp.sigma_g2,",
      "\"sigma_e2\" => hsq_snp.sigma_e2,",
      "\"loglik\" => hsq_snp.loglik,",
      "\"converged\" => hsq_snp.converged,",
      "\"nobs\" => length(hsq_y)",
      ");"
    )),
    hint = hs_dense_scale_hint
  )

  raw <- JuliaCall::julia_eval("hsq_snp_raw")
  estimated_vc <- c(
    sigma_g2 = as.numeric(raw$sigma_g2),
    sigma_e2 = as.numeric(raw$sigma_e2)
  )
  result <- hs_normalize_julia_snp_blup_result(
    raw,
    payload,
    estimated_vc,
    provenance = "estimated_snp_blup_reml",
    converged = isTRUE(raw$converged),
    loglik = as.numeric(raw$loglik)
  )
  # df for AIC/BIC: the fixed effects + the two REML-estimated variance
  # components (sigma_g2, sigma_e2). The marker effects are random (BLUP), not
  # free parameters. (The supplied-variance path estimates no VCs, so no df.)
  result$df <- as.integer(ncol(payload$X) + 2L)
  hs_new_fit(
    spec = list(
      method = "SNP-BLUP-REML",
      family = list(family = payload$family, link = "identity"),
      target = "snp_blup"
    ),
    payload = payload,
    result = result,
    engine = "HSquared.jl"
  )
}

hs_normalize_julia_snp_blup_result <- function(
  raw,
  payload,
  variance_components,
  provenance = "supplied",
  converged = TRUE,
  loglik = NULL
) {
  sigma_g2 <- unname(variance_components[["sigma_g2"]])
  sigma_e2 <- unname(variance_components[["sigma_e2"]])

  fixed_effects <- as.numeric(raw$beta)
  fixed_names <- payload$metadata$fixed_colnames
  if (length(fixed_effects) == length(fixed_names)) {
    names(fixed_effects) <- fixed_names
  }

  ids <- as.character(payload$ids)
  genomic_bv <- data.frame(id = ids, value = as.numeric(raw$gebv))

  effects <- as.numeric(raw$marker_effects)
  marker_labels <- payload$marker_names
  if (is.null(marker_labels) || length(marker_labels) != length(effects)) {
    marker_labels <- as.character(seq_along(effects))
  }
  markers_rec <- if (!is.null(payload$Z)) {
    as.matrix(payload$Z %*% payload$markers)
  } else {
    payload$markers
  }
  marker_allele_frequencies <- if (!is.null(raw$p)) as.numeric(raw$p) else NULL

  list(
    variance_components = data.frame(
      component = c("genomic", "residual"),
      estimate = c(sigma_g2, sigma_e2)
    ),
    heritability = data.frame(
      term = "genomic",
      estimate = sigma_g2 / (sigma_g2 + sigma_e2)
    ),
    breeding_values = genomic_bv,
    fixed_effects = fixed_effects,
    marker_effects = data.frame(
      marker = as.character(marker_labels),
      effect = effects
    ),
    marker_allele_frequencies = marker_allele_frequencies,
    marker_variance_explained = hs_marker_variance_explained_from_snp_blup(
      effects = effects,
      markers = markers_rec,
      marker_labels = marker_labels,
      allele_frequencies = marker_allele_frequencies
    ),
    random_effects = list(genomic = genomic_bv),
    predictions = data.frame(.fitted = as.numeric(raw$fitted)),
    nobs = as.integer(raw$nobs),
    loglik = loglik,
    diagnostics = list(
      target = "snp_blup",
      variance_components = provenance,
      optimizer_status = if (identical(provenance, "supplied")) {
        "not_run"
      } else if (isTRUE(converged)) {
        "converged"
      } else {
        "not_converged"
      },
      n_markers = length(effects)
    ),
    converged = isTRUE(converged)
  )
}

hs_marker_variance_explained_from_snp_blup <- function(
  effects,
  markers,
  marker_labels = NULL,
  allele_frequencies = NULL
) {
  effects <- as.numeric(effects)
  markers <- as.matrix(markers)

  if (!is.numeric(markers) || ncol(markers) != length(effects)) {
    stop(
      "Internal bridge error: marker effects and marker matrix columns are ",
      "not aligned.",
      call. = FALSE
    )
  }
  if (nrow(markers) < 1L || ncol(markers) < 1L) {
    stop("Internal bridge error: marker matrix is empty.", call. = FALSE)
  }
  if (any(!is.finite(effects)) || any(!is.finite(markers))) {
    stop(
      "Internal bridge error: marker variance explained requires finite ",
      "marker effects and marker dosages.",
      call. = FALSE
    )
  }
  if (!is.null(allele_frequencies)) {
    allele_frequencies <- as.numeric(allele_frequencies)
    if (
      length(allele_frequencies) != ncol(markers) ||
        any(!is.finite(allele_frequencies))
    ) {
      stop(
        "Internal bridge error: allele frequencies must align with marker ",
        "columns.",
        call. = FALSE
      )
    }
    if (any(allele_frequencies < 0 | allele_frequencies > 1)) {
      stop(
        "Internal bridge error: allele frequencies must lie in [0, 1].",
        call. = FALSE
      )
    }
  }

  if (is.null(marker_labels) || length(marker_labels) != length(effects)) {
    marker_labels <- colnames(markers)
  }
  if (is.null(marker_labels) || length(marker_labels) != length(effects)) {
    marker_labels <- as.character(seq_along(effects))
  }

  center <- if (is.null(allele_frequencies)) {
    colMeans(markers)
  } else {
    2 * allele_frequencies
  }
  centered <- sweep(markers, 2L, center, check.margin = FALSE)
  centered_marker_variance <- colMeans(centered^2)
  contribution <- centered_marker_variance * effects^2
  total <- sum(contribution)
  proportion <- if (isTRUE(total > 0)) {
    contribution / total
  } else {
    rep(NA_real_, length(contribution))
  }

  data.frame(
    marker = as.character(marker_labels),
    effect = effects,
    centered_marker_variance = as.numeric(centered_marker_variance),
    contribution = as.numeric(contribution),
    proportion = as.numeric(proportion),
    stringsAsFactors = FALSE
  )
}

hs_julia_setup <- function(project) {
  project <- normalizePath(project, winslash = "/", mustWork = TRUE)
  if (
    isTRUE(hs_julia_bridge_state$initialized) &&
      identical(hs_julia_bridge_state$project, project)
  ) {
    return(invisible(TRUE))
  }

  JuliaCall::julia_setup(installJulia = FALSE, verbose = FALSE)
  JuliaCall::julia_assign("hsq_project", project)
  JuliaCall::julia_command(
    "using Pkg; Pkg.activate(hsq_project); using HSquared; using SparseArrays;"
  )
  hs_julia_bridge_state$initialized <- TRUE
  hs_julia_bridge_state$project <- project
  invisible(TRUE)
}

hs_julia_assign_payload <- function(payload, initial) {
  JuliaCall::julia_assign("hsq_y", payload$y)
  JuliaCall::julia_assign("hsq_X", payload$X)
  hs_julia_assign_sparse_csc("hsq_Z", payload$Z)
  JuliaCall::julia_assign("hsq_id", payload$pedigree$id)
  JuliaCall::julia_assign(
    "hsq_sire",
    hs_parent_for_julia(payload$pedigree$sire)
  )
  JuliaCall::julia_assign("hsq_dam", hs_parent_for_julia(payload$pedigree$dam))
  JuliaCall::julia_assign("hsq_method", payload$method)
  JuliaCall::julia_assign("hsq_initial_sigma_a2", unname(initial[["sigma_a2"]]))
  JuliaCall::julia_assign("hsq_initial_sigma_e2", unname(initial[["sigma_e2"]]))
  invisible(TRUE)
}

hs_julia_assign_sparse_csc <- function(name, x) {
  slots <- hs_sparse_csc_slots(x)
  JuliaCall::julia_assign(paste0(name, "_nrow"), slots$nrow)
  JuliaCall::julia_assign(paste0(name, "_ncol"), slots$ncol)
  JuliaCall::julia_assign(paste0(name, "_colptr"), slots$colptr)
  JuliaCall::julia_assign(paste0(name, "_rowval"), slots$rowval)
  JuliaCall::julia_assign(paste0(name, "_nzval"), slots$nzval)
  JuliaCall::julia_command(paste0(
    name,
    " = HSquared.sparse_csc_matrix(",
    name,
    "_nrow, ",
    name,
    "_ncol, ",
    name,
    "_colptr, ",
    name,
    "_rowval, ",
    name,
    "_nzval; index_base = :zero);"
  ))
  invisible(TRUE)
}

hs_sparse_csc_slots <- function(x) {
  if (!inherits(x, "dgCMatrix")) {
    stop("`x` must be a `Matrix::dgCMatrix` object.", call. = FALSE)
  }

  list(
    nrow = as.integer(nrow(x)),
    ncol = as.integer(ncol(x)),
    colptr = as.integer(x@p),
    rowval = as.integer(x@i),
    nzval = as.numeric(x@x)
  )
}

hs_parent_for_julia <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- "0"
  x
}

# R-side NA -> Julia NaN sentinel for multivariate Y.
# Keeps the R payload's NA cells intact; only the Julia assign copy is rewritten.
# Needed because JuliaCall 0.17.6 marshals NA_real_ to Missing, not NaN.
hs_y_matrix_for_julia <- function(Y) {
  if (!is.matrix(Y)) {
    stop("`Y` must be a matrix.", call. = FALSE)
  }
  out <- Y
  storage.mode(out) <- "double"
  out[is.na(out)] <- NaN
  out
}

# hsquared#225: `initial` for the non-Gaussian target (`fit_laplace_reml()`'s
# single-variance-component families -- poisson/bernoulli/binomial) is a list
# with `sigma_a2` only, matching the engine's `initial = (sigma_a2 = ...,)`
# NamedTuple shape. `NULL` (the default) means "let Julia centre its own
# log-scale search bracket, `log(sigma_a2) +/- 6`, at its hard-coded
# sigma_a2 = 1.0".
hs_validate_nongaussian_initial <- function(initial) {
  if (is.null(initial)) {
    return(NULL)
  }
  if (is.null(names(initial)) || !"sigma_a2" %in% names(initial)) {
    stop(
      "`initial` for the non-Gaussian target must be a list with `sigma_a2`.",
      call. = FALSE
    )
  }
  out <- initial[["sigma_a2"]]
  out <- suppressWarnings(as.numeric(out))
  if (length(out) != 1L || !is.finite(out) || out <= 0) {
    stop(
      "`initial$sigma_a2` must be a single positive finite value.",
      call. = FALSE
    )
  }
  out
}

# hsquared#225: `restart_check` opts into `fit_laplace_reml()`'s two-start
# restart (HSquared.jl#327), which refits once from a bumped start and flags
# `boundary = TRUE` when the estimate moves with the start. Default FALSE
# matches the engine's own default.
hs_validate_restart_check <- function(restart_check) {
  if (
    !is.logical(restart_check) || length(restart_check) != 1L ||
      is.na(restart_check)
  ) {
    stop(
      "`engine_control$restart_check` must be a single TRUE/FALSE value.",
      call. = FALSE
    )
  }
  isTRUE(restart_check)
}

hs_validate_initial_variances <- function(initial) {
  if (
    is.null(names(initial)) ||
      !all(c("sigma_a2", "sigma_e2") %in% names(initial))
  ) {
    stop(
      "`initial` must include `sigma_a2` and `sigma_e2`.",
      call. = FALSE
    )
  }
  out <- as.numeric(initial[c("sigma_a2", "sigma_e2")])
  names(out) <- c("sigma_a2", "sigma_e2")
  if (any(!is.finite(out)) || any(out <= 0)) {
    stop(
      "`initial` variance values must be positive and finite.",
      call. = FALSE
    )
  }
  out
}

# hsquared#212: `initial` for the multi-effect (K independent random-effect
# blocks + residual) target is a plain positive numeric vector of length
# K + 1, in block order (animal first, then each i.i.d. block, then the
# residual) -- matching `HSquared.fit_multi_effect_reml()`'s own `initial`
# argument shape (a vector, not a named list; there is no per-block name to
# attach since K is caller-determined).
hs_validate_multi_effect_initial <- function(initial, k) {
  if (!is.numeric(initial) || length(initial) != k) {
    stop(
      "`initial` for the multi-effect target must be a numeric vector of ",
      "length ",
      k,
      " (one value per random-effect block, in formula order, plus the ",
      "residual).",
      call. = FALSE
    )
  }
  if (any(!is.finite(initial)) || any(initial <= 0)) {
    stop(
      "`initial` variance components must be finite and positive.",
      call. = FALSE
    )
  }
  as.numeric(initial)
}

# hsquared#212: `initial` for the direct-maternal (correlated 2x2 G_dm) target
# is a list with `G_dm` (a 2x2 positive-definite matrix) and `sigma_e2` (a
# positive scalar) -- matching `HSquared.fit_direct_maternal_reml()`'s own
# `initial = (G_dm = ..., sigma_e2 = ...)` NamedTuple shape.
hs_validate_direct_maternal_initial <- function(initial) {
  if (!is.list(initial) || !all(c("G_dm", "sigma_e2") %in% names(initial))) {
    stop(
      "`initial` for the direct-maternal target must be a list with ",
      "`G_dm` (a 2x2 matrix) and `sigma_e2` (a positive scalar).",
      call. = FALSE
    )
  }
  G_dm <- as.matrix(initial[["G_dm"]])
  if (!identical(dim(G_dm), c(2L, 2L)) || any(!is.finite(G_dm))) {
    stop("`initial$G_dm` must be a finite 2x2 matrix.", call. = FALSE)
  }
  storage.mode(G_dm) <- "double"
  if (
    !isTRUE(all.equal(G_dm, t(G_dm))) ||
      eigen(G_dm, only.values = TRUE)$values[2L] <= 0
  ) {
    stop(
      "`initial$G_dm` must be a symmetric positive-definite 2x2 matrix.",
      call. = FALSE
    )
  }
  sigma_e2 <- as.numeric(initial[["sigma_e2"]])
  if (length(sigma_e2) != 1L || !is.finite(sigma_e2) || sigma_e2 <= 0) {
    stop(
      "`initial$sigma_e2` must be a single positive finite value.",
      call. = FALSE
    )
  }
  list(G_dm = G_dm, sigma_e2 = sigma_e2)
}

# hsquared#212: `engine_control` keys silently discarded by a target LOOK
# accepted (no error, no warning) but never reach the Julia call. This
# declares, per `target`, which keys the bridge actually forwards, and errors
# -- naming the key and the target -- when a caller supplies one that is not
# honoured. `target` (the routing key itself) and `julia_project` (read by
# every target's dispatch to resolve the Julia project) are always allowed and
# are not repeated in the table below. Call once per `hsquared()` Julia
# dispatch, before any `hs_fit_julia_*_payload` builder runs.
hs_engine_control_honoured_keys <- list(
  fit_animal_model = c("initial", "max_dense_cells"),
  henderson_mme = "variance_components",
  metafounder = "variance_components",
  sparse_reml = c("initial", "iterations"),
  ai_reml = c("initial", "iterations", "em_warmup"),
  repeatability = c("initial", "iterations", "max_dense_cells", "scale_method"),
  two_effect = c("initial", "iterations"),
  # multi_effect: initial/iterations are honoured on the `scale_method =
  # "dense"` (default) route only; `scale_method = "auto"` does not forward
  # them (HSquared.jl#343, a known remaining gap, not fixed here).
  multi_effect = c("initial", "iterations", "scale_method"),
  direct_maternal = c("initial", "iterations"),
  genomic = c("initial", "iterations"),
  single_step = c("initial", "iterations"),
  single_step_construct = c("initial", "iterations"),
  metafounder_single_step = c("initial", "iterations"),
  snp_blup = "variance_components",
  relmat = c("initial", "iterations"),
  precision = c("initial", "iterations"),
  multivariate = c("initial", "iterations", "genetic_structure", "rank"),
  random_regression = "iterations",
  # initial (hsquared#225): a list with `sigma_a2`, the centre of the
  # engine's log-scale search bracket, `log(sigma_a2) +/- 6` (default centre
  # sigma_a2 = 1.0). restart_check (hsquared#225): opt-in two-start
  # refit that flags `boundary = TRUE` when the estimate moves with the start
  # (HSquared.jl#327).
  nongaussian = c("marginal", "iterations", "initial", "restart_check")
)

hs_engine_control_forwarding <- function(control, target) {
  honoured <- c(
    hs_engine_control_honoured_keys[[target]],
    "target",
    "julia_project"
  )
  supplied <- names(control$engine_control)
  unsupported <- setdiff(supplied, honoured)
  if (length(unsupported) > 0L) {
    forwarded <- setdiff(honoured, c("target", "julia_project"))
    hs_abort_unsupported_syntax(
      "`engine_control` key",
      if (length(unsupported) > 1L) "s " else " ",
      paste(sprintf("`%s`", unsupported), collapse = ", "),
      if (length(unsupported) > 1L) " are" else " is",
      " not honoured by `target = \"",
      target,
      "\"`. Supported keys (besides `target`/`julia_project`): ",
      if (length(forwarded) == 0L) {
        "(none)"
      } else {
        paste(sprintf("`%s`", forwarded), collapse = ", ")
      },
      "."
    )
  }
  forwarded <- setdiff(honoured, c("target", "julia_project"))
  control$engine_control[intersect(supplied, forwarded)]
}

hs_validate_iterations <- function(iterations) {
  iterations <- suppressWarnings(as.integer(iterations))
  if (length(iterations) != 1L || is.na(iterations) || iterations <= 0L) {
    stop("`iterations` must be a single positive integer.", call. = FALSE)
  }
  iterations
}

# em_warmup: opt-in EM-REML warm-start iterations before the AI/Newton step (engine
# `fit_ai_reml`, V1-AI-REML). 0 (default) = off / byte-identical to the pre-warm-start path.
hs_validate_em_warmup <- function(em_warmup) {
  em_warmup <- suppressWarnings(as.integer(em_warmup))
  if (length(em_warmup) != 1L || is.na(em_warmup) || em_warmup < 0L) {
    stop("`em_warmup` must be a single non-negative integer.", call. = FALSE)
  }
  em_warmup
}

# max_dense_cells: the R-facing lever for the engine's dense-validation size
# guard (hsquared#214, #217). Bounds `nobs^2 + nanimals^2` on the dense
# fitters (the default `animal()` route via `fit_animal_model()` /
# `fit_variance_components()`, and `target = "repeatability"`); default
# 1e6 mirrors the engine's own `DEFAULT_MAX_DENSE_CELLS` unchanged.
hs_validate_max_dense_cells <- function(max_dense_cells) {
  max_dense_cells <- suppressWarnings(as.numeric(max_dense_cells))
  if (
    length(max_dense_cells) != 1L ||
      is.na(max_dense_cells) ||
      max_dense_cells != as.integer(max_dense_cells) ||
      max_dense_cells <= 0
  ) {
    stop(
      "`engine_control$max_dense_cells` must be a single positive integer.",
      call. = FALSE
    )
  }
  max_dense_cells
}

# Production fence (0.9): allowed `engine_control$target` values below are
# validation-scale live routes only — not production sparse fitting. Default
# `engine = "fit"` auto-selects ai_reml / genomic / multivariate without listing
# them here. payload_v2 block routing is separate (direct_maternal, multi_effect).
# FA/lowrank stay blocked in hs_validate_genetic_structure_control(). PATH_ONLY
# interval smoke (C1-ext) is not a target here. See design-45 DRAFT.
hs_validate_julia_target <- function(target) {
  if (!is.character(target) || length(target) != 1L || is.na(target)) {
    stop(
      "`engine_control$target` must be a single string.",
      call. = FALSE
    )
  }
  # D2: matrix-free REML is engine-only (V1-MATFREE-REML). Do not add an R
  # target; name the closest live large-scale siblings instead of a bare
  # allowlist dump.
  if (
    target %in%
      c("matrix_free", "matrix_free_reml", "matfree", "matrix_free_mc_em_reml")
  ) {
    hs_abort_unsupported_syntax(
      "`engine_control$target = \"",
      target,
      "\"` is not an R bridge target. HSquared.jl `fit_matrix_free_reml` ",
      "(V1-MATFREE-REML) stays engine-only for experimental 0.9.0. Closest ",
      "live R paths: `target = \"ai_reml\"` / `\"sparse_reml\"` for the ",
      "univariate animal model, or `target = \"multi_effect\"` / ",
      "`\"repeatability\"` with `scale_method = \"auto\"` for the related ",
      "experimental large-scale Monte-Carlo route. public_covered_count ",
      "stays 7.",
      call. = FALSE
    )
  }
  if (
    !target %in%
      c(
        "fit_animal_model",
        "henderson_mme",
        "metafounder",
        "sparse_reml",
        "ai_reml",
        "repeatability",
        "two_effect",
        "multi_effect",
        "genomic",
        "single_step",
        "single_step_construct",
        "metafounder_single_step",
        "snp_blup",
        "relmat",
        "precision",
        "multivariate",
        "random_regression",
        "nongaussian",
        "direct_maternal"
      )
  ) {
    hs_abort_unsupported_syntax(
      "`engine_control$target` must be one of \"fit_animal_model\", ",
      "\"henderson_mme\", \"sparse_reml\", \"ai_reml\", \"repeatability\", ",
      "\"metafounder\", \"two_effect\", \"multi_effect\", \"genomic\", ",
      "\"single_step\", ",
      "\"single_step_construct\", \"metafounder_single_step\", \"snp_blup\", ",
      "\"relmat\", \"precision\", \"multivariate\", ",
      "\"random_regression\", \"nongaussian\", or \"direct_maternal\".",
      call. = FALSE
    )
  }
  target
}

hs_validate_genetic_structure_control <- function(control, target) {
  value <- hs_engine_control_value(control, "genetic_structure", NULL)
  if (is.null(value)) {
    return("unstructured")
  }
  if (!is.character(value) || length(value) != 1L || is.na(value)) {
    stop(
      "`engine_control$genetic_structure` must be a single string.",
      call. = FALSE
    )
  }
  allowed <- c("unstructured", "diagonal", "lowrank", "factor_analytic")
  if (!value %in% allowed) {
    hs_abort_unsupported_syntax(
      "`engine_control$genetic_structure` must be one of \"unstructured\", ",
      "\"diagonal\", \"lowrank\", or \"factor_analytic\"."
    )
  }
  if (!identical(target, "multivariate")) {
    hs_abort_unsupported_syntax(
      "`engine_control$genetic_structure` is only planned for the ",
      "`target = \"multivariate\"` bridge. Remove `genetic_structure`, or use ",
      "`target = \"multivariate\"` with a `cbind(...)` response."
    )
  }
  if (identical(value, "factor_analytic")) {
    hs_abort_unsupported_syntax(
      "`genetic_structure = \"factor_analytic\"` is planned on the R ",
      "surface and not activated on the R bridge. Julia V4-FA is ",
      "engine-covered (HSquared.jl 60895208 / #300); that is not an ",
      "R-public factor-analytic fit. The rotation convention is already ",
      "ratified (rotation-invariant functionals only, never loadings). ",
      "Use `genetic_structure = \"unstructured\"` or `\"diagonal\"`. ",
      "public_covered_count stays 7."
    )
  }
  if (identical(value, "lowrank")) {
    hs_abort_unsupported_syntax(
      "`genetic_structure = \"lowrank\"` is planned, not activated on ",
      "the R bridge. The opt-in multivariate path estimates ",
      "`\"unstructured\"` or `\"diagonal\"` G0; use one of those. ",
      "public_covered_count stays 7."
    )
  }
  # "unstructured" (default) and "diagonal" are both reachable. "diagonal" has
  # no loadings and no rotation ambiguity (it is just per-trait genetic
  # variances with zero genetic covariances), so it is honesty-clean to surface
  # ahead of lowrank/factor_analytic.
  rank <- hs_engine_control_value(control, "rank", NULL)
  if (!is.null(rank)) {
    if (
      !is.numeric(rank) ||
        length(rank) != 1L ||
        is.na(rank) ||
        !is.finite(rank) ||
        rank < 1L ||
        rank != as.integer(rank)
    ) {
      stop(
        "`engine_control$rank` must be a single positive integer.",
        call. = FALSE
      )
    }
    hs_abort_unsupported_syntax(
      "`engine_control$rank` is reserved for future `lowrank` and ",
      "`factor_analytic` structured covariance controls. The current ",
      "multivariate bridge estimates unstructured or diagonal G0 with ",
      "unstructured R0 only; remove `rank` until low-rank or ",
      "factor-analytic support is available."
    )
  }
  value
}

# All opt-in engine targets that can fit a given non-default random effect. A
# genomic marker primary fits either by GREML on the built relationship
# (`genomic`) or as a marker-effect model (`snp_blup`; supplied-variance, or
# REML-estimated when variances are omitted).
hs_effect_targets <- function(type) {
  switch(
    type,
    permanent = "repeatability",
    common_env = "two_effect",
    # maternal_genetic supports both the INDEPENDENT two-effect target (default
    # suggestion) and the CORRELATED direct-maternal target (opt-in Phase 4).
    maternal_genetic = c("two_effect", "direct_maternal"),
    iid_effects = "multi_effect",
    metafounder = "metafounder",
    genomic = c("genomic", "snp_blup"),
    single_step = c(
      "single_step",
      "single_step_construct",
      "metafounder_single_step"
    ),
    # Experimental supplied-relationship primary effects: each fits only through
    # its own opt-in supplied-relationship-inverse target.
    relmat = "relmat",
    precision = "precision",
    stop("Unknown random effect type: ", type, call. = FALSE)
  )
}

# Map a parsed second random effect to the opt-in engine target that fits it.
# Map a parsed non-default random effect (a second effect, or the genomic
# primary effect) to the opt-in engine target that fits it.
hs_second_effect_target <- function(type) {
  switch(
    type,
    permanent = "repeatability",
    common_env = "two_effect",
    # Routing map only: two_effect stays the single suggested *name* for
    # callers that interpolate one target. Default-path user copy must use
    # hs_maternal_genetic_default_path_message() so the covered
    # direct_maternal sibling is named. Do not auto-route.
    maternal_genetic = "two_effect",
    iid_effects = "multi_effect",
    metafounder = "metafounder",
    genomic = "genomic",
    single_step = "single_step",
    relmat = "relmat",
    precision = "precision",
    stop("Unknown random effect type: ", type, call. = FALSE)
  )
}

# Default-path copy for maternal_genetic(). Names the covered correlated
# sibling first, then the experimental independent two_effect path. Does
# not change routing: neither target is auto-selected on engine = "fit".
hs_maternal_genetic_default_path_message <- function() {
  paste0(
    "`maternal_genetic()` is not on the default path. ",
    "Two opt-in targets fit this formula; neither is auto-routed.\n\n",
    "Closest working call (covered correlated 2x2 G / Willham triple, ",
    "not a scalar h2):\n\n",
    "  control = hs_control(engine = \"julia\", engine_control = list(",
    "target = \"direct_maternal\"))\n\n",
    "Experimental alternative (independent maternal):\n\n",
    "  control = hs_control(engine = \"julia\", engine_control = list(",
    "target = \"two_effect\"))\n\n",
    "See formula_status() and the current-limits article."
  )
}

hs_abort_maternal_genetic_default_path <- function() {
  hs_abort_unsupported_syntax(
    hs_maternal_genetic_default_path_message(),
    call. = FALSE
  )
}

# Extra clause for a wrong-target maternal_genetic() formula on engine = "julia".
hs_maternal_genetic_allowed_targets_note <- function() {
  paste0(
    " `target = \"direct_maternal\"` is covered (correlated 2x2 G / Willham ",
    "triple, not a scalar h2); `target = \"two_effect\"` is experimental ",
    "(independent maternal)."
  )
}

hs_validate_supplied_variances <- function(
  variance_components,
  target = "henderson_mme"
) {
  if (is.null(variance_components)) {
    stop(
      "`engine_control$variance_components` is required when ",
      "`target = \"",
      target,
      "\"`.",
      call. = FALSE
    )
  }
  if (
    is.null(names(variance_components)) ||
      !all(c("sigma_a2", "sigma_e2") %in% names(variance_components))
  ) {
    stop(
      "`engine_control$variance_components` must include `sigma_a2` and ",
      "`sigma_e2`.",
      call. = FALSE
    )
  }
  out <- as.numeric(variance_components[c("sigma_a2", "sigma_e2")])
  names(out) <- c("sigma_a2", "sigma_e2")
  if (any(!is.finite(out)) || any(out <= 0)) {
    stop(
      "`engine_control$variance_components` values must be positive and ",
      "finite.",
      call. = FALSE
    )
  }
  out
}

hs_validate_snp_blup_variances <- function(variance_components) {
  if (is.null(variance_components)) {
    stop(
      "`engine_control$variance_components` is required when ",
      "`target = \"snp_blup\"` (named `sigma_g2` and `sigma_e2`).",
      call. = FALSE
    )
  }
  if (
    is.null(names(variance_components)) ||
      !all(c("sigma_g2", "sigma_e2") %in% names(variance_components))
  ) {
    stop(
      "`engine_control$variance_components` must include `sigma_g2` (genomic) ",
      "and `sigma_e2` (residual) for `target = \"snp_blup\"`.",
      call. = FALSE
    )
  }
  out <- as.numeric(variance_components[c("sigma_g2", "sigma_e2")])
  names(out) <- c("sigma_g2", "sigma_e2")
  if (any(!is.finite(out)) || any(out <= 0)) {
    stop(
      "`engine_control$variance_components` values must be positive and ",
      "finite.",
      call. = FALSE
    )
  }
  out
}

hs_plot_data_try <- function(expr) {
  tryCatch(expr, error = function(e) NULL)
}

hs_plot_data_list <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  x <- hs_drop_julia_classes(x)
  if (is.data.frame(x)) {
    x <- as.list(x)
  }
  if (!is.list(x)) {
    return(NULL)
  }
  x
}

hs_plot_data_character <- function(x, n = NULL) {
  if (is.null(x)) {
    return(NULL)
  }
  out <- as.character(x)
  if (!is.null(n) && length(out) != n) {
    return(NULL)
  }
  out
}

hs_plot_data_numeric <- function(x, n = NULL) {
  if (is.null(x)) {
    return(NULL)
  }
  out <- as.numeric(x)
  if (!is.null(n) && length(out) != n) {
    return(NULL)
  }
  out
}

hs_plot_data_scalar <- function(x) {
  if (is.null(x) || length(x) < 1L) {
    return(NULL)
  }
  x[[1L]]
}

hs_plot_data_matrix <- function(x, nr = NULL, nc = NULL) {
  if (is.null(x)) {
    return(NULL)
  }
  out <- as.matrix(x)
  storage.mode(out) <- "double"
  if (!is.null(nr) && !is.null(nc)) {
    if (!identical(dim(out), c(nr, nc))) {
      if (length(out) != nr * nc) {
        return(NULL)
      }
      out <- matrix(as.numeric(out), nrow = nr, ncol = nc)
    }
  } else if (!is.null(nr) && nrow(out) != nr) {
    if (length(out) %% nr != 0L) {
      return(NULL)
    }
    out <- matrix(as.numeric(out), nrow = nr)
  }
  out
}

hs_normalize_variance_components_plot_data <- function(pd) {
  pd <- hs_plot_data_list(pd)
  if (is.null(pd)) {
    return(NULL)
  }
  term <- hs_plot_data_character(pd$term)
  estimate <- hs_plot_data_numeric(pd$estimate)
  if (is.null(term) || is.null(estimate) || length(term) != length(estimate)) {
    return(NULL)
  }
  n <- length(term)
  out <- list(term = term, estimate = estimate)
  lo <- hs_plot_data_numeric(pd$lo, n)
  hi <- hs_plot_data_numeric(pd$hi, n)
  panel <- hs_plot_data_character(pd$panel, n)
  if (!is.null(lo)) {
    out$lo <- lo
  }
  if (!is.null(hi)) {
    out$hi <- hi
  }
  if (!is.null(panel)) {
    out$panel <- panel
  }
  interval_status <- hs_plot_data_scalar(pd$interval_status)
  interval_method <- hs_plot_data_scalar(pd$interval_method)
  if (!is.null(interval_status)) {
    out$interval_status <- as.character(interval_status)
  }
  if (!is.null(interval_method)) {
    out$interval_method <- as.character(interval_method)
  }
  out
}

hs_normalize_breeding_values_plot_data <- function(pd) {
  pd <- hs_plot_data_list(pd)
  if (is.null(pd)) {
    return(NULL)
  }
  id <- hs_plot_data_character(pd$id %||% pd$ids)
  value <- hs_plot_data_numeric(
    pd$value %||% pd$values %||% pd$breeding_value %||% pd$breeding_values
  )
  if (is.null(id) || is.null(value) || length(id) != length(value)) {
    return(NULL)
  }
  n <- length(id)
  out <- list(id = id, value = value)
  trait <- hs_plot_data_character(pd$trait, n)
  pev <- hs_plot_data_numeric(pd$pev %||% pd$prediction_error_variance, n)
  if (
    !is.null(trait) &&
      !(length(unique(trait)) == 1L && unique(trait) %in% c("1", "trait_1"))
  ) {
    out$trait <- trait
  }
  if (!is.null(pev)) {
    out$pev <- pev
  }
  pev_scale <- hs_plot_data_scalar(pd$pev_scale)
  if (!is.null(pev_scale)) {
    out$pev_scale <- as.character(pev_scale)
  }
  out
}

hs_attach_standard_plot_data <- function(result, raw) {
  vcpd <- hs_plot_data_try(
    hs_normalize_variance_components_plot_data(
      raw$variance_components_plot_data
    )
  )
  if (!is.null(vcpd)) {
    result$variance_components_plot_data <- vcpd
  }
  bvpd <- hs_plot_data_try(
    hs_normalize_breeding_values_plot_data(raw$breeding_values_plot_data)
  )
  if (!is.null(bvpd)) {
    result$breeding_values_plot_data <- bvpd
  }
  result
}

hs_normalize_genetic_correlation_plot_data <- function(pd, traits = NULL) {
  pd <- hs_plot_data_list(pd)
  if (is.null(pd) || is.null(pd$genetic_correlations)) {
    return(NULL)
  }
  rg <- hs_plot_data_matrix(pd$genetic_correlations)
  if (is.null(rg) || nrow(rg) != ncol(rg)) {
    return(NULL)
  }
  pd_traits <- hs_plot_data_character(pd$traits)
  traits <- pd_traits %||% traits
  if (is.null(traits) || length(traits) != nrow(rg)) {
    traits <- paste0("trait_", seq_len(nrow(rg)))
  }
  dimnames(rg) <- list(traits, traits)
  out <- list(
    traits = as.character(traits),
    genetic_correlations = rg,
    rotation_invariant = isTRUE(hs_plot_data_scalar(pd$rotation_invariant))
  )
  h2 <- hs_plot_data_numeric(pd$heritabilities, length(traits))
  if (!is.null(h2)) {
    out$heritabilities <- h2
  }
  out
}

hs_normalize_genetic_pca_plot_data <- function(pd) {
  pd <- hs_plot_data_list(pd)
  if (is.null(pd)) {
    return(NULL)
  }
  eigenvalues <- hs_plot_data_numeric(pd$eigenvalues)
  if (is.null(eigenvalues)) {
    return(NULL)
  }
  n <- length(eigenvalues)
  out <- list(
    eigenvalues = eigenvalues,
    rotation_invariant = isTRUE(hs_plot_data_scalar(pd$rotation_invariant)),
    is_eigenstructure_not_loadings = isTRUE(hs_plot_data_scalar(
      pd$is_eigenstructure_not_loadings
    ))
  )
  variance_explained <- hs_plot_data_numeric(pd$variance_explained, n)
  axis_labels <- hs_plot_data_character(pd$axis_labels, n)
  if (!is.null(variance_explained)) {
    out$variance_explained <- variance_explained
  }
  if (!is.null(axis_labels)) {
    out$axis_labels <- axis_labels
  }
  out
}

hs_attach_multivariate_plot_data <- function(result, raw, traits = NULL) {
  gcpd <- hs_plot_data_try(
    hs_normalize_genetic_correlation_plot_data(
      raw$genetic_correlation_plot_data,
      traits
    )
  )
  if (!is.null(gcpd)) {
    result$genetic_correlation_plot_data <- gcpd
  }
  gppd <- hs_plot_data_try(
    hs_normalize_genetic_pca_plot_data(raw$genetic_pca_plot_data)
  )
  if (!is.null(gppd)) {
    result$genetic_pca_plot_data <- gppd
  }
  result
}

hs_rr_payload_covariate <- function(pd, rr) {
  cov <- hs_plot_data_numeric(pd$covariate)
  if (is.null(cov)) {
    return(NULL)
  }
  hs_unstandardize_covariate(cov, rr$lower, rr$upper)
}

hs_normalize_rr_genetic_variance_plot_data <- function(pd, rr) {
  pd <- hs_plot_data_list(pd)
  if (is.null(pd)) {
    return(NULL)
  }
  cov <- hs_rr_payload_covariate(pd, rr)
  value <- hs_plot_data_numeric(pd$value %||% pd$genetic_variance)
  if (is.null(cov) || is.null(value) || length(cov) != length(value)) {
    return(NULL)
  }
  n <- length(cov)
  out <- list(covariate = cov, value = value, genetic_variance = value)
  h2 <- hs_plot_data_numeric(pd$heritability, n)
  if (!is.null(h2)) {
    out$heritability <- h2
  }
  basis_order <- hs_plot_data_scalar(pd$basis_order)
  supplied <- hs_plot_data_scalar(pd$supplied)
  if (!is.null(basis_order)) {
    out$basis_order <- as.integer(basis_order)
  }
  if (!is.null(supplied)) {
    out$supplied <- isTRUE(supplied)
  }
  out
}

hs_normalize_rr_eigenfunctions_plot_data <- function(pd, rr) {
  pd <- hs_plot_data_list(pd)
  if (is.null(pd)) {
    return(NULL)
  }
  cov <- hs_rr_payload_covariate(pd, rr)
  if (is.null(cov)) {
    return(NULL)
  }
  eigenfunctions <- hs_plot_data_matrix(pd$eigenfunctions, nr = length(cov))
  if (is.null(eigenfunctions)) {
    return(NULL)
  }
  k <- ncol(eigenfunctions)
  out <- list(
    covariate = cov,
    eigenfunctions = eigenfunctions,
    rotation_invariant = isTRUE(hs_plot_data_scalar(pd$rotation_invariant))
  )
  variance_explained <- hs_plot_data_numeric(pd$variance_explained, k)
  axis <- hs_plot_data_numeric(pd$axis, k)
  if (!is.null(variance_explained)) {
    out$variance_explained <- variance_explained
  }
  if (!is.null(axis)) {
    out$axis <- as.integer(axis)
  }
  out
}

hs_normalize_rr_covariance_surface_plot_data <- function(pd, rr) {
  pd <- hs_plot_data_list(pd)
  if (is.null(pd)) {
    return(NULL)
  }
  cov <- hs_rr_payload_covariate(pd, rr)
  if (is.null(cov)) {
    return(NULL)
  }
  surface <- hs_plot_data_matrix(pd$surface, nr = length(cov), nc = length(cov))
  if (is.null(surface)) {
    return(NULL)
  }
  list(
    covariate = cov,
    surface = surface,
    is_correlation = isTRUE(hs_plot_data_scalar(pd$is_correlation))
  )
}

hs_attach_random_regression_plot_data <- function(result, raw) {
  rr <- result$random_regression
  gvpd <- hs_plot_data_try(
    hs_normalize_rr_genetic_variance_plot_data(
      raw$rr_genetic_variance_plot_data,
      rr
    )
  )
  if (!is.null(gvpd)) {
    result$rr_genetic_variance_plot_data <- gvpd
  }
  efpd <- hs_plot_data_try(
    hs_normalize_rr_eigenfunctions_plot_data(
      raw$rr_eigenfunctions_plot_data,
      rr
    )
  )
  if (!is.null(efpd)) {
    result$rr_eigenfunctions_plot_data <- efpd
  }
  sfpd <- hs_plot_data_try(
    hs_normalize_rr_covariance_surface_plot_data(
      raw$rr_covariance_surface_plot_data,
      rr
    )
  )
  if (!is.null(sfpd)) {
    result$rr_covariance_surface_plot_data <- sfpd
  }
  result
}

hs_normalize_julia_result <- function(raw, payload) {
  breeding_values <- raw$breeding_values
  animal <- raw$random_effects$animal

  fixed_effects <- as.numeric(raw$fixed_effects)
  fixed_names <- payload$metadata$fixed_colnames
  if (length(fixed_effects) == length(fixed_names)) {
    names(fixed_effects) <- fixed_names
  }

  result <- list(
    variance_components = data.frame(
      component = c("animal", "residual"),
      estimate = c(
        raw$variance_components$sigma_a2,
        raw$variance_components$sigma_e2
      )
    ),
    heritability = data.frame(
      term = "animal",
      estimate = raw$heritability
    ),
    breeding_values = data.frame(
      id = as.character(breeding_values$ids),
      value = as.numeric(breeding_values$values)
    ),
    fixed_effects = fixed_effects,
    random_effects = list(
      animal = data.frame(
        id = as.character(animal$ids),
        value = as.numeric(animal$values)
      )
    ),
    loglik = as.numeric(raw$loglik),
    df = as.integer(raw$df),
    nobs = as.integer(raw$nobs),
    predictions = data.frame(.fitted = as.numeric(raw$predictions)),
    diagnostics = hs_drop_julia_classes(raw$diagnostics),
    converged = isTRUE(raw$converged)
  )

  if (!is.null(raw$prediction_error_variance)) {
    result$prediction_error_variance <- hs_julia_id_values(
      raw$prediction_error_variance
    )
  }
  if (!is.null(raw$reliability)) {
    result$reliability <- hs_julia_id_values(raw$reliability)
  }
  if (!is.null(raw$heritability_interval)) {
    result$heritability_interval <- hs_normalize_heritability_interval(
      raw$heritability_interval
    )
  }
  if (!is.null(raw$variance_component_se)) {
    result$variance_component_se <- hs_normalize_variance_component_se(
      raw$variance_component_se
    )
  }
  if (!is.null(raw$heritability_se)) {
    result$heritability_se <- as.numeric(raw$heritability_se)
  }

  result <- hs_attach_standard_plot_data(result, raw)
  result
}

# Normalize the engine's heritability_interval NamedTuple
# (heritability, lower, upper, level, [se], method) into a one-row data frame.
# `se` is absent for the profile method.
hs_normalize_heritability_interval <- function(hi) {
  data.frame(
    estimate = as.numeric(hi$heritability),
    lower = as.numeric(hi$lower),
    upper = as.numeric(hi$upper),
    level = as.numeric(hi$level),
    se = if (!is.null(hi$se)) as.numeric(hi$se) else NA_real_,
    method = as.character(hi$method),
    stringsAsFactors = FALSE
  )
}

# Normalize the engine's variance_component_standard_errors NamedTuple
# (sigma_a2, sigma_e2) into a data frame of component standard errors.
hs_normalize_variance_component_se <- function(vcse) {
  data.frame(
    component = c("animal", "residual"),
    se = c(as.numeric(vcse$sigma_a2), as.numeric(vcse$sigma_e2)),
    stringsAsFactors = FALSE
  )
}

# Normalize the engine's repeatability_interval NamedTuple
# (repeatability, lower, upper, level, se) into a one-row data frame. This is the
# delta-method CI for the repeatability coefficient t = (Va + Vpe)/Vp (engine row
# V3-REPEAT-REML, partial), not for h2; the engine offers only the delta method
# for t, so there is no `method` column.
hs_normalize_repeatability_interval <- function(ri) {
  data.frame(
    estimate = as.numeric(ri$repeatability),
    lower = as.numeric(ri$lower),
    upper = as.numeric(ri$upper),
    level = as.numeric(ri$level),
    se = as.numeric(ri$se),
    stringsAsFactors = FALSE
  )
}

hs_normalize_julia_henderson_mme_result <- function(
  raw,
  payload,
  variance_components
) {
  fixed_effects <- as.numeric(raw$fixed_effects)
  fixed_names <- payload$metadata$fixed_colnames
  if (length(fixed_effects) == length(fixed_names)) {
    names(fixed_effects) <- fixed_names
  }

  animal <- data.frame(
    id = as.character(raw$animal_ids),
    value = as.numeric(raw$animal_effects)
  )

  result <- list(
    variance_components = data.frame(
      component = c("animal", "residual"),
      estimate = c(
        unname(variance_components[["sigma_a2"]]),
        unname(variance_components[["sigma_e2"]])
      )
    ),
    heritability = data.frame(
      term = "animal",
      estimate = unname(variance_components[["sigma_a2"]]) /
        sum(unname(variance_components))
    ),
    breeding_values = animal,
    fixed_effects = fixed_effects,
    random_effects = list(animal = animal),
    predictions = data.frame(.fitted = as.numeric(raw$fitted)),
    nobs = as.integer(raw$nobs),
    diagnostics = list(
      target = "henderson_mme",
      variance_components = "supplied",
      optimizer_status = "not_run"
    ),
    converged = TRUE
  )

  if (!is.null(raw$prediction_error_variance)) {
    result$prediction_error_variance <- hs_julia_id_values(
      raw$prediction_error_variance
    )
  }
  if (!is.null(raw$reliability)) {
    result$reliability <- hs_julia_id_values(raw$reliability)
  }

  result
}

hs_julia_id_values <- function(x) {
  data.frame(
    id = as.character(x$ids),
    value = as.numeric(x$values)
  )
}

hs_drop_julia_classes <- function(x) {
  if (is.list(x)) {
    x <- lapply(x, hs_drop_julia_classes)
    class(x) <- NULL
  }
  x
}
