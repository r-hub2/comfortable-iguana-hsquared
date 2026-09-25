test_that("hs_control stores validated defaults", {
  control <- hs_control()

  expect_s3_class(control, "hs_control")
  expect_equal(control$engine, "fit")
  expect_equal(control$backend, "auto")
  expect_equal(control$accelerator, "auto")
  expect_equal(control$precision, "float64")
  expect_equal(control$save, "minimal")
  expect_equal(control$engine_control, list())
})

test_that("hs_control validates engine selection", {
  control <- hs_control(
    engine = "julia",
    engine_control = list(
      initial = c(sigma_a2 = 1, sigma_e2 = 1),
      target = "henderson_mme",
      variance_components = c(sigma_a2 = 1.2, sigma_e2 = 0.8)
    )
  )

  expect_equal(control$engine, "julia")
  expect_equal(control$engine_control$initial, c(sigma_a2 = 1, sigma_e2 = 1))
  expect_equal(control$engine_control$target, "henderson_mme")
  expect_equal(
    control$engine_control$variance_components,
    c(sigma_a2 = 1.2, sigma_e2 = 0.8)
  )

  expect_error(
    hs_control(engine = "not-an-engine"),
    "'arg' should be one of",
    fixed = TRUE
  )
})

test_that("hs_control preserves planned backend vocabulary", {
  expect_equal(hs_control(backend = "threads")$backend, "threads")
  expect_equal(hs_control(backend = "cuda")$backend, "cuda")
  expect_equal(hs_control(backend = "amdgpu")$backend, "amdgpu")
  expect_equal(hs_control(backend = "metal")$backend, "metal")
  expect_equal(hs_control(backend = "oneapi")$backend, "oneapi")
  expect_equal(hs_control(accelerator = "gpu")$accelerator, "gpu")
  expect_equal(hs_control(accelerator = "metal")$accelerator, "metal")
  expect_equal(hs_control(accelerator = "amdgpu")$accelerator, "amdgpu")
  expect_equal(hs_control(accelerator = "oneapi")$accelerator, "oneapi")

  expect_error(hs_control(backend = "tpu"), "'arg' should be one of")
  expect_error(hs_control(accelerator = "tpu"), "'arg' should be one of")
})

test_that("backend_info separates control vocabulary from execution", {
  info <- backend_info(hs_control(backend = "metal"))

  expect_s3_class(info, "hs_backend_info")
  expect_equal(
    info$backend,
    c("cpu", "threads", "cuda", "amdgpu", "metal", "oneapi")
  )
  expect_true(info$requested[info$backend == "metal"])
  expect_true(all(info$selectable))
  expect_false(any(info$execution_available))
  expect_true(all(info$status == "planned"))

  gpu_info <- backend_info(hs_control(accelerator = "gpu"))
  expect_true(all(
    gpu_info$requested[
      gpu_info$backend %in%
        c(
          "cuda",
          "amdgpu",
          "metal",
          "oneapi"
        )
    ]
  ))

  expect_error(
    backend_info(control = list()),
    "`control` must be created by `hs_control\\(\\)`."
  )
})

test_that("formula_status separates parsed, reserved, and planned grammar", {
  status <- formula_status()

  expect_s3_class(status, "hs_formula_status")
  expect_equal(nrow(status), 36L)
  expect_true("term" %in% names(status))
  expect_true("syntax_status" %in% names(status))
  expect_true("fitting_status" %in% names(status))
  expect_equal(
    status$syntax_status[status$term == "animal(1 | id, pedigree = ped)"],
    "parsed"
  )
  expect_equal(
    status$syntax_status[
      status$term == "animal(1 | id) with data = hs_data(..., pedigree = ped)"
    ],
    "parsed"
  )
  expect_true(all(
    status$fitting_status[status$syntax_status != "parsed"] == "not available"
  ))
  expect_true("permanent(1 | id)" %in% status$term)
  expect_true("genomic(1 | id, Ginv = Ginv)" %in% status$term)
  rr_term <- "animal(rr(covariate, order = 2) | id, pedigree = ped)"
  expect_true(rr_term %in% status$term)
  expect_equal(status$syntax_status[status$term == rr_term], "parsed")
  expect_equal(
    status$fitting_status[status$term == rr_term],
    "fitted (opt-in random-regression)"
  )
  expect_equal(
    status$fitting_status[
      status$term == "cbind(trait1, trait2) ~ animal(1 | id, pedigree = ped)"
    ],
    "fitted (default route, covered multivariate)"
  )
  expect_true(all(
    c(
      "animal(trait | id, pedigree = ped, cov = us())",
      "animal(trait | id, pedigree = ped, cov = diag())",
      "animal(trait | id, pedigree = ped, cov = lowrank(K = 2))",
      "animal(trait | id, pedigree = ped, cov = fa(K = 2))"
    ) %in%
      status$term
  ))
  expect_equal(
    status$syntax_status[
      status$term == "animal(trait | id, pedigree = ped, cov = us())"
    ],
    "planned"
  )
  expect_equal(
    status$syntax_status[
      status$term == "animal(trait | id, pedigree = ped, cov = diag())"
    ],
    "planned"
  )
  expect_equal(
    status$syntax_status[
      status$term == "animal(trait | id, pedigree = ped, cov = lowrank(K = 2))"
    ],
    "planned"
  )
  expect_equal(
    status$syntax_status[
      status$term == "animal(trait | id, pedigree = ped, cov = fa(K = 2))"
    ],
    "reserved"
  )
  expect_match(
    status$current_behavior[
      status$term == "animal(trait | id, pedigree = ped, cov = fa(K = 2))"
    ],
    "60895208",
    fixed = TRUE
  )
  expect_match(
    status$current_behavior[
      status$term == "animal(trait | id, pedigree = ped, cov = fa(K = 2))"
    ],
    "public_covered_count stays 7",
    fixed = TRUE
  )
  expect_true(all(
    c(
      "missing = miss_control(response = \"include\")",
      "mi(x) with missing = miss_control(predictor = \"model\")"
    ) %in%
      status$term
  ))
  expect_equal(
    status$syntax_status[
      status$term == "mi(x) with missing = miss_control(predictor = \"model\")"
    ],
    "planned"
  )
  expect_match(
    status$current_behavior[
      status$term == "mi(x) with missing = miss_control(predictor = \"model\")"
    ],
    "no mi\\(\\) marker is exported"
  )
  expect_true(any(status$syntax_status == "planned"))
  expect_match(capture.output(print(status))[[1L]], "<hs_formula_status>")
  expect_match(
    paste(capture.output(print(status)), collapse = "\n"),
    "permanent\\(1 \\| id\\)"
  )
  expect_match(
    paste(capture.output(print(status)), collapse = "\n"),
    "reserved/planned:",
    fixed = TRUE
  )
  expect_true(any(
    status$term == "single_step(1 | id, pedigree = ped, markers = M)" &
      status$fitting_status == "fitted (opt-in single-step construction)"
  ))
  single_step_construct_row <- status[
    status$term == "single_step(1 | id, pedigree = ped, markers = M)",
  ]
  expect_match(
    single_step_construct_row$current_behavior,
    "HSquared.jl #295",
    fixed = TRUE
  )
  expect_match(
    single_step_construct_row$current_behavior,
    "public_covered_count stays 7",
    fixed = TRUE
  )
  expect_match(
    single_step_construct_row$current_behavior,
    "not a default-route promotion",
    fixed = TRUE
  )
  single_step_bundle_row <- status[
    status$term ==
      paste0(
        "single_step(1 | id) with data = hs_data(..., pedigree = ped, ",
        "genotypes = M)"
      ),
  ]
  expect_equal(single_step_bundle_row$syntax_status, "parsed")
  expect_equal(
    single_step_bundle_row$fitting_status,
    "fitted (opt-in single-step bundle construction)"
  )
  expect_match(
    single_step_bundle_row$current_behavior,
    "genotypes",
    fixed = TRUE
  )
  expect_match(
    single_step_bundle_row$current_behavior,
    "target = \"single_step_construct\"",
    fixed = TRUE
  )
  expect_true(any(
    grepl("Gamma = Gamma", status$term, fixed = TRUE) &
      status$syntax_status == "parsed" &
      status$fitting_status == "fitted (opt-in supplied-Gamma H^Gamma)"
  ))
  metafounder_row <- status[
    status$term ==
      "metafounder(1 | id, pedigree = ped, group = group, Gamma = Gamma)",
  ]
  expect_equal(metafounder_row$syntax_status, "parsed")
  expect_equal(
    metafounder_row$fitting_status,
    "fitted (opt-in supplied-Gamma metafounder)"
  )
  expect_match(metafounder_row$current_behavior, "Gamma", fixed = TRUE)
  expect_match(
    metafounder_row$current_behavior,
    "supplied, not estimated",
    fixed = TRUE
  )
  subset <- status[
    status$category == "multivariate and factor analytic",
    c("term", "syntax_status", "fitting_status")
  ]
  expect_s3_class(subset, "hs_formula_status")
  expect_output(print(subset), "cov = lowrank")
})

test_that("formula_status helper vectors stay the same length", {
  n <- lengths(list(
    term = hsquared:::hs_formula_status_terms(),
    category = hsquared:::hs_formula_status_categories(),
    phase = hsquared:::hs_formula_status_phases(),
    syntax_status = hsquared:::hs_formula_status_syntax(),
    fitting_status = hsquared:::hs_formula_status_fitting(),
    current_behavior = hsquared:::hs_formula_status_behavior()
  ))
  expect_equal(length(unique(n)), 1L)
  expect_equal(unname(unique(n)), nrow(formula_status()))
})

test_that("formula_status print header is derived from the printed rows", {
  status <- formula_status()
  printed <- paste(capture.output(print(status)), collapse = "\n")

  expect_match(printed, "default:")
  expect_match(printed, "animal\\(\\)")
  expect_match(printed, "cbind\\(\\)")
  expect_match(printed, "opt-in:")
  expect_match(printed, "animal\\(rr\\(\\)\\)")
  expect_match(printed, "single_step\\(\\)")
  expect_match(printed, "metafounder\\(\\)")
  expect_match(printed, "relmat\\(\\)")
  expect_match(printed, "precision\\(\\)")
  expect_match(printed, "\\(1 \\| group\\)")
  expect_match(printed, "\\$current_behavior")
  expect_false(grepl(
    "permanent/common_env/maternal_genetic/genomic fit opt-in",
    printed,
    fixed = TRUE
  ))
  expect_false(grepl("parsed today:", printed, fixed = TRUE))

  rr <- status[
    status$term == "animal(rr(covariate, order = 2) | id, pedigree = ped)",
  ]
  class(rr) <- unique(c("hs_formula_status", class(rr)))
  rr_print <- paste(capture.output(print(rr)), collapse = "\n")
  expect_match(rr_print, "animal\\(rr\\(\\)\\)")
  expect_false(grepl("cbind()", rr_print, fixed = TRUE))
  expect_false(grepl("default:", rr_print, fixed = TRUE))

  planned <- status[status$syntax_status == "planned", ]
  class(planned) <- unique(c("hs_formula_status", class(planned)))
  planned_print <- paste(capture.output(print(planned)), collapse = "\n")
  expect_match(planned_print, "reserved/planned:")
  expect_false(grepl("default:", planned_print, fixed = TRUE))
  expect_false(grepl("opt-in:", planned_print, fixed = TRUE))
})

test_that("validation_status separates evidence from planned validation", {
  status <- validation_status()

  expect_s3_class(status, "hs_validation_status")
  expect_equal(nrow(status), 21L)
  expect_equal(
    status$status[
      status$capability ==
        "experimental supplied-relationship estimator (opt-in: genomic, single-step)"
    ],
    "covered"
  )
  supplied_row <- status[
    status$capability ==
      "experimental supplied-relationship estimator (opt-in: genomic, single-step)",
  ]
  expect_match(supplied_row$evidence, "single_step_construct", fixed = TRUE)
  expect_match(supplied_row$claim_boundary, "genomic GREML", fixed = TRUE)
  expect_match(
    supplied_row$claim_boundary,
    "public_covered_count is 7",
    fixed = TRUE
  )
  expect_equal(
    status$status[
      status$capability == "experimental repeatability estimator (opt-in)"
    ],
    "partial"
  )
  expect_equal(
    status$status[
      status$capability == "experimental multivariate REML estimator (opt-in)"
    ],
    "covered"
  )
  # Phase 1 (2026-07-01): common-environment two-effect leg COVERED (public_covered_count 1->2).
  # Phase 2-R (2026-07-01): generalized to arbitrary-N independent (1|g) effects, same row,
  # COVERED (public_covered_count 2->3; engine V3-NEFFECT-REML gate + sommer comparator + exact
  # live R-Julia parity). Maternal leg (A2=pedigree) stays experimental (label reflects the split).
  expect_equal(
    status$status[
      status$capability ==
        "two-effect / arbitrary-N independent-effect estimator (opt-in; covered: common-env + (1|g) iid / A2=I; experimental: maternal / A2=pedigree)"
    ],
    "covered"
  )
  expect_true(all(
    c("capability", "phase", "status", "evidence", "claim_boundary") %in%
      names(status)
  ))
  expect_equal(
    status$status[
      status$capability == "supplied-variance Henderson MME fixture"
    ],
    "partial"
  )
  expect_equal(
    status$status[status$capability == "sparse REML likelihood identity"],
    "partial"
  )
  expect_equal(
    status$status[
      status$capability == "Mrode-style supplied-variance outputs"
    ],
    "partial"
  )
  expect_equal(
    status$status[
      status$capability == "experimental sparse REML estimator (opt-in)"
    ],
    "partial"
  )
  default_row <- status[
    status$capability ==
      "univariate Gaussian animal-model fit (default path, AI-REML)",
  ]
  expect_equal(default_row$status, "covered")
  expect_match(default_row$claim_boundary, "REML = FALSE", fixed = TRUE)
  expect_match(default_row$claim_boundary, "default fit", fixed = TRUE)
  expect_match(
    default_row$claim_boundary,
    'engine = "validate"',
    fixed = TRUE
  )
  expect_false(grepl(
    "rejected on the fit path",
    default_row$claim_boundary,
    fixed = TRUE
  ))
  expect_equal(
    status$status[
      status$capability ==
        "external published-REML recovery (gryphon, R reference)"
    ],
    "covered"
  )
  expect_equal(
    status$status[
      status$capability ==
        "known-truth DGP variance-component recovery (R reference)"
    ],
    "covered"
  )
  expect_equal(
    status$status[status$capability == "ASReml comparison policy"],
    "planned"
  )
  multivariate_row <- status[
    status$capability == "experimental multivariate REML estimator (opt-in)",
  ]
  expect_equal(multivariate_row$status, "covered")
  expect_match(
    multivariate_row$evidence,
    "100-replicate cold-start t=2 known-truth recovery study",
    fixed = TRUE
  )
  expect_match(
    multivariate_row$evidence,
    "full-unstructured residual sommer comparator",
    fixed = TRUE
  )
  expect_match(
    multivariate_row$evidence,
    "published Mrode Example 5.1 multiple-trait supplied-G0/R0 BLUP/MME",
    fixed = TRUE
  )
  expect_match(
    multivariate_row$evidence,
    "Bayesian MCMCglmm agreement probe",
    fixed = TRUE
  )
  expect_match(
    multivariate_row$evidence,
    "explicitly not a same-estimand REML comparator",
    fixed = TRUE
  )
  expect_match(
    multivariate_row$claim_boundary,
    "This R public surface is COVERED at",
    fixed = TRUE
  )
  expect_match(
    multivariate_row$claim_boundary,
    "published Mrode-style supplied-variance BLUP/MME anchor",
    fixed = TRUE
  )
  expect_match(
    multivariate_row$claim_boundary,
    "Twin C8 broader-DGP confirm is banked",
    fixed = TRUE
  )
  expect_match(
    multivariate_row$claim_boundary,
    "14/16 pass",
    fixed = TRUE
  )
  expect_match(
    multivariate_row$claim_boundary,
    "MV-5 is SUPERSEDED",
    fixed = TRUE
  )
  expect_match(
    multivariate_row$claim_boundary,
    "hs_require_suggests",
    fixed = TRUE
  )
  expect_match(
    multivariate_row$claim_boundary,
    "scoped to k = 2 unstructured",
    fixed = TRUE
  )
  expect_match(
    multivariate_row$claim_boundary,
    "k >= 3 stays parseable-and-fittable-but-experimental",
    fixed = TRUE
  )
  expect_match(
    multivariate_row$claim_boundary,
    "\"diagonal\" stays experimental",
    fixed = TRUE
  )
  expect_match(
    multivariate_row$claim_boundary,
    "MCMCglmm leg is not same-estimand REML parity",
    fixed = TRUE
  )
  expect_match(
    status$claim_boundary[
      status$capability == "CPU/GPU backend comparison"
    ],
    "no backend execution",
    fixed = TRUE
  )
  expect_match(capture.output(print(status))[[1L]], "<hs_validation_status>")
  expect_match(
    paste(capture.output(print(status)), collapse = "\n"),
    "supplied-variance Henderson MME fixture"
  )
})

test_that("capability ids stay stable and labels carry current wording", {
  status <- validation_status()

  expect_true("capability_label" %in% names(status))
  expect_type(status$capability_label, "character")
  expect_length(status$capability_label, nrow(status))
  expect_false(anyNA(status$capability_label))
  expect_false(any(!nzchar(status$capability_label)))

  # Every override must name a live capability id. Without this, renaming an id
  # would silently strand its alias and the label would revert to the stale
  # wording with no test failing.
  overrides <- hs_validation_status_label_overrides()
  expect_true(all(names(overrides) %in% status$capability))

  # Rows with no override are unaliased: label is the id verbatim.
  unaliased <- !(status$capability %in% names(overrides))
  expect_identical(
    status$capability_label[unaliased],
    status$capability[unaliased]
  )
})

test_that("the multivariate row is looked up by its historical opt-in id", {
  status <- validation_status()

  # Dated evidence cites this id verbatim - docs/dev-log/comparator-runs/
  # 2026-06-21-multivariate-tool-availability.md and
  # 2026-09-01-blupf90-tool-unavailability.md. Lookup by it must keep working.
  historical_id <- "experimental multivariate REML estimator (opt-in)"
  row <- status[status$capability == historical_id, ]
  expect_equal(nrow(row), 1L)
  expect_equal(row$status, "covered")

  # MV-4 made the route default, so the reader-facing label drops "opt-in"
  # while the id keeps it.
  expect_false(grepl("opt-in", row$capability_label, fixed = TRUE))
  expect_match(row$capability_label, "default route", fixed = TRUE)

  # The printed table shows the label, not the stale id.
  printed <- paste(capture.output(print(status)), collapse = "\n")
  expect_match(
    printed,
    "multivariate REML estimator (default route)",
    fixed = TRUE
  )
})

test_that("hs_control validates engine_control", {
  expect_error(
    hs_control(engine_control = "not-a-list"),
    "`engine_control` must be a list.",
    fixed = TRUE
  )

  expect_error(
    hs_control(engine_control = list(1)),
    "`engine_control` must be a named list.",
    fixed = TRUE
  )
})

test_that("hsquared validates basic call shape", {
  expect_error(hsquared(), "`formula` is required.", fixed = TRUE)

  expect_error(hsquared(y ~ x), "`data` is required.", fixed = TRUE)

  expect_error(
    hsquared(y ~ x, data = data.frame(y = 1, x = 1), control = list()),
    "`control` must be created by `hs_control()`.",
    fixed = TRUE
  )
})

test_that("the validate engine validates and returns the spec without fitting", {
  ped <- data.frame(
    id = c("a", "b", "c", "d"),
    sire = c(NA, NA, "a", "a"),
    dam = c(NA, NA, "b", "c")
  )
  dat <- data.frame(
    y = c(1, 2, 3),
    sex = c("f", "m", "f"),
    age = c(1, 2, 3),
    id = c("a", "c", "d")
  )

  # `engine = "validate"` confirms the contract with a message and returns the
  # validated spec invisibly (it no longer stops), so it can be inspected.
  expect_message(
    spec <- hsquared(
      y ~ sex + age + animal(1 | id, pedigree = ped),
      data = dat,
      control = hs_control(engine = "validate")
    ),
    "Validated the v0.1 animal-model contract",
    fixed = TRUE
  )
  expect_type(spec, "list")
  expect_match(spec$bridge$target, "fit_animal_model", fixed = TRUE)
})

test_that("the default engine fits, and errors clearly without the Julia engine", {
  ped <- data.frame(
    id = c("a", "b", "c"),
    sire = c(NA, NA, "a"),
    dam = c(NA, NA, "b")
  )
  dat <- data.frame(
    y = c(1, 2, 3),
    sex = c("m", "f", "m"),
    id = c("a", "b", "c")
  )

  # Default engine = "fit"; with no Julia engine available it errors with
  # actionable install guidance rather than silently doing nothing.
  expect_error(
    hsquared(
      y ~ sex + animal(1 | id, pedigree = ped),
      data = dat,
      control = hs_control(
        engine = "fit",
        engine_control = list(julia_project = tempfile())
      )
    ),
    "requires the HSquared.jl Julia",
    fixed = TRUE
  )
})

test_that("the default fit path rejects REML = FALSE rather than mislabeling ML", {
  ped <- data.frame(
    id = c("a", "b", "c"),
    sire = c(NA, NA, "a"),
    dam = c(NA, NA, "b")
  )
  dat <- data.frame(
    y = c(1, 2, 3),
    sex = c("m", "f", "m"),
    id = c("a", "b", "c")
  )

  # The default fit path estimates variance components by REML only. A
  # `REML = FALSE` request must be rejected honestly (ML is not implemented),
  # not silently run as REML and returned mislabeled as "ML". This is a pure
  # request-validity error, so it fires before any Julia-engine check.
  expect_error(
    hsquared(
      y ~ sex + animal(1 | id, pedigree = ped),
      data = dat,
      family = stats::gaussian(),
      REML = FALSE
    ),
    "ML estimation",
    fixed = TRUE
  )
})
