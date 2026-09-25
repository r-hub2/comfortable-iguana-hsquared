# design-41 §3 #8 — Hinv-cell R↔engine parity receipt.
# Covered cell only: supplied-Hinv from G = A22 + 0.05 I (design-56 §A.2).
# Not VanRaden markers=. Not a covered flip. Count stays 7.

hs_ss_hinv_fixture_dir <- function() {
  testthat::test_path("fixtures", "ss_hinv_parity")
}

hs_ss_hinv_read_labeled <- function(path) {
  raw <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  ids <- names(raw)[-1L]
  mat <- as.matrix(raw[, -1L, drop = FALSE])
  storage.mode(mat) <- "double"
  dimnames(mat) <- list(as.character(raw[[1L]]), ids)
  mat
}

hs_ss_hinv_julia_project <- function() {
  env <- Sys.getenv("HSQUARED_JULIA_PROJECT", unset = "")
  if (nzchar(env) && file.exists(file.path(env, "Project.toml"))) {
    return(normalizePath(env, winslash = "/", mustWork = FALSE))
  }
  candidates <- c(
    path.expand("~/local-scratch/lanes/HSquared.jl-08-ss-20260903"),
    file.path(
      dirname(normalizePath(testthat::test_path("..", ".."), mustWork = FALSE)),
      "HSquared.jl"
    ),
    "/Users/z3437171/Dropbox/Github Local/HSquared.jl"
  )
  for (p in candidates) {
    if (file.exists(file.path(p, "Project.toml"))) {
      return(normalizePath(p, winslash = "/", mustWork = FALSE))
    }
  }
  hsquared:::hs_default_julia_project()
}

test_that("Hinv-cell fixture pins the teaching-kernel packet (Julia-free)", {
  dir <- hs_ss_hinv_fixture_dir()
  sha <- readLines(file.path(dir, "ENGINE_PACKET_SHA.txt"), warn = FALSE)
  expect_equal(trimws(sha[[1L]]), "0b03d67e")

  meta <- utils::read.csv(
    file.path(dir, "metadata.csv"),
    stringsAsFactors = FALSE
  )
  kv <- stats::setNames(meta$value, meta$key)
  expect_equal(unname(kv[["n"]]), "6")
  expect_equal(unname(kv[["g_shift"]]), "0.05")
  expect_equal(unname(kv[["tau"]]), "1")
  expect_equal(unname(kv[["omega"]]), "1")
  expect_equal(unname(kv[["blend_weight"]]), "0")
  expect_equal(unname(kv[["ridge"]]), "0")
  expect_equal(unname(kv[["mrode_ch11_anchor"]]), "NO_ANCHOR")
  expect_equal(unname(kv[["not_a_covered_flip"]]), "true")

  hinv <- hs_ss_hinv_read_labeled(file.path(dir, "engine_hinv.csv"))
  expect_equal(dim(hinv), c(6L, 6L))
  expect_equal(rownames(hinv), c("s1", "d1", "d2", "o1", "o2", "o3"))
  expect_equal(colnames(hinv), rownames(hinv))
  expect_true(isTRUE(all.equal(hinv, t(hinv), tolerance = 1e-12)))

  g <- hs_ss_hinv_read_labeled(file.path(dir, "G.csv"))
  expect_equal(dim(g), c(3L, 3L))
  expect_equal(rownames(g), c("o1", "o2", "o3"))
  # teaching kernel: G = A22 + 0.05 I, not VanRaden
  a <- hs_ss_hinv_read_labeled(file.path(dir, "A.csv"))
  a22 <- a[c("o1", "o2", "o3"), c("o1", "o2", "o3")]
  expect_equal(g, a22 + diag(0.05, 3L), tolerance = 1e-12)
})

test_that("supplied-Hinv teaching-kernel cell matches the engine [live]", {
  testthat::skip_on_cran()
  project <- hs_ss_hinv_julia_project()
  testthat::skip_if_not(
    hsquared:::hs_julia_bridge_available(project),
    "JuliaCall, Julia, and a local HSquared.jl are required for Hinv-cell parity."
  )
  hsquared:::hs_julia_setup(project)

  dir <- hs_ss_hinv_fixture_dir()
  pinned_hinv <- hs_ss_hinv_read_labeled(file.path(dir, "engine_hinv.csv"))
  ids <- rownames(pinned_hinv)
  ped <- data.frame(
    id = ids,
    sire = c(NA, NA, NA, "s1", "s1", "s1"),
    dam = c(NA, NA, NA, "d1", "d1", "d2"),
    stringsAsFactors = FALSE
  )

  JuliaCall::julia_assign("hsq_ids", ids)
  JuliaCall::julia_assign("hsq_sire", c("0", "0", "0", "s1", "s1", "s1"))
  JuliaCall::julia_assign("hsq_dam", c("0", "0", "0", "d1", "d1", "d2"))
  JuliaCall::julia_assign("hsq_grows", as.integer(c(4L, 5L, 6L)))
  JuliaCall::julia_command(paste(
    "using LinearAlgebra;",
    "hsq_ped = HSquared.normalize_pedigree(hsq_ids, hsq_sire, hsq_dam);",
    "hsq_Ainv = Matrix(HSquared.pedigree_inverse(hsq_ped));",
    "hsq_A = Matrix(HSquared.additive_relationship(hsq_ped));",
    "hsq_G = hsq_A[hsq_grows, hsq_grows] + 0.05 * Matrix{Float64}(I, 3, 3);",
    "hsq_Hinv = Matrix(HSquared.single_step_inverse(",
    "hsq_Ainv, hsq_A, hsq_G, hsq_grows;",
    "tau = 1.0, omega = 1.0, blend_weight = 0.0, ridge = 0.0));"
  ))
  live_hinv <- JuliaCall::julia_eval("hsq_Hinv")
  dimnames(live_hinv) <- list(ids, ids)
  expect_equal(live_hinv, pinned_hinv, tolerance = 1e-10)

  # Same y on both sides. Three records / animal so the REML cell is identified
  # enough for a parity comparison. This is NOT the n=6 recovery-smoke draw and
  # is NOT a recovery-to-truth gate.
  set.seed(20260903L)
  rec <- rep(ids, each = 3L)
  dat <- data.frame(
    y = 2 +
      as.numeric(factor(rec, levels = ids)) * 0.15 +
      stats::rnorm(18L, 0, 0.35),
    id = rec,
    stringsAsFactors = FALSE
  )

  fit_r <- hsquared(
    y ~ single_step(1 | id, Hinv = pinned_hinv),
    data = dat,
    family = stats::gaussian(),
    control = hs_control(
      engine = "julia",
      engine_control = list(
        target = "single_step",
        julia_project = project,
        initial = c(sigma_a2 = 0.5, sigma_e2 = 1.0),
        iterations = 200L
      )
    )
  )
  expect_s3_class(fit_r, "hsquared_fit")
  expect_equal(fit_r$spec$target, "single_step")
  expect_false(identical(fit_r$spec$target, "single_step_construct"))

  vc_r <- variance_components(fit_r)
  sa_r <- vc_r$estimate[vc_r$component == "single_step"]
  se_r <- vc_r$estimate[vc_r$component == "residual"]
  expect_length(sa_r, 1L)
  expect_length(se_r, 1L)
  expect_true(is.finite(sa_r) && is.finite(se_r))

  bv_r <- breeding_values(fit_r)
  expect_setequal(as.character(bv_r$id), ids)

  y <- as.numeric(dat$y)
  X <- matrix(1, nrow = length(y), ncol = 1L)
  Z <- matrix(0, nrow = length(y), ncol = length(ids))
  Z[cbind(seq_along(y), match(dat$id, ids))] <- 1
  JuliaCall::julia_assign("hsq_y", y)
  JuliaCall::julia_assign("hsq_X", X)
  JuliaCall::julia_assign("hsq_Z", Z)
  JuliaCall::julia_command(paste(
    "hsq_fit = HSquared.fit_single_step_reml(",
    "hsq_y, hsq_X, hsq_Z, hsq_Ainv, hsq_A, hsq_G, hsq_grows;",
    "tau = 1.0, omega = 1.0, blend_weight = 0.0, ridge = 0.0,",
    "initial = (sigma_a2 = 0.5, sigma_e2 = 1.0), ids = hsq_ids);",
    "hsq_bv = HSquared.breeding_values(hsq_fit);"
  ))
  sa_j <- JuliaCall::julia_eval("hsq_fit.variance_components.sigma_a2")
  se_j <- JuliaCall::julia_eval("hsq_fit.variance_components.sigma_e2")
  bv_ids <- as.character(JuliaCall::julia_eval("collect(hsq_bv.ids)"))
  bv_vals <- as.numeric(JuliaCall::julia_eval("collect(hsq_bv.values)"))

  expect_equal(as.numeric(sa_r), as.numeric(sa_j), tolerance = 1e-8)
  expect_equal(as.numeric(se_r), as.numeric(se_j), tolerance = 1e-8)
  expect_equal(sort(as.character(bv_r$id)), sort(bv_ids))
  ord <- match(ids, as.character(bv_r$id))
  j_ord <- match(ids, bv_ids)
  expect_equal(
    as.numeric(bv_r$value[ord]),
    as.numeric(bv_vals[j_ord]),
    tolerance = 1e-8
  )
})
