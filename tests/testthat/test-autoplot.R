test_that("theme_hsquared returns a ggplot2 theme", {
  th <- theme_hsquared()
  expect_s3_class(th, "theme")
})

mock_uni_fit <- function() {
  structure(
    list(
      result = list(
        variance_components = data.frame(
          component = c("animal", "residual"),
          estimate = c(0.6, 0.9)
        ),
        variance_component_se = data.frame(
          component = c("animal", "residual"),
          se = c(0.12, 0.10)
        ),
        heritability = data.frame(term = "animal", estimate = 0.4),
        heritability_se = 0.07,
        breeding_values = data.frame(
          id = letters[1:6],
          value = c(-1, -0.4, -0.1, 0.2, 0.5, 1)
        ),
        prediction_error_variance = data.frame(
          id = letters[1:6],
          value = c(0.2, 0.25, 0.22, 0.21, 0.24, 0.2)
        )
      )
    ),
    class = "hsquared_fit"
  )
}

mock_mv_fit <- function() {
  g <- matrix(
    c(1, 0.3, 0.3, 1),
    2,
    2,
    dimnames = list(c("t1", "t2"), c("t1", "t2"))
  )
  structure(
    list(result = list(genetic_correlation = g)),
    class = "hsquared_fit"
  )
}

mock_mv_fit_h2 <- function(h2 = c(t1 = 0.05, t2 = 0.4)) {
  f <- mock_mv_fit()
  f$result$heritability <- data.frame(
    term = names(h2),
    estimate = as.numeric(h2),
    stringsAsFactors = FALSE
  )
  f
}

# A fit carrying the engine `genetic_correlation_plot_data` payload. This mock
# exercises the auto-detect branch without requiring a live Julia bridge. The
# payload correlation (0.55) differs from the `genetic_correlation` field (0.3)
# to prove consumption.
mock_mv_fit_plotdata <- function() {
  structure(
    list(
      result = list(
        genetic_correlation = matrix(
          c(1, 0.3, 0.3, 1),
          2,
          2,
          dimnames = list(c("t1", "t2"), c("t1", "t2"))
        ),
        genetic_correlation_plot_data = list(
          traits = c("t1", "t2"),
          genetic_correlations = matrix(c(1, 0.55, 0.55, 1), 2, 2),
          heritabilities = c(0.05, 0.4),
          rotation_invariant = TRUE
        )
      )
    ),
    class = "hsquared_fit"
  )
}

# A fit carrying the engine Set-B `variance_components_plot_data` payload (twin
# PR #95). This exercises the auto-detect branch without a live Julia bridge.
mock_uni_fit_vcpd <- function(h2_hi = 0.54) {
  structure(
    list(
      result = list(
        variance_components_plot_data = list(
          term = c("sigma_a2", "sigma_e2", "h2"),
          estimate = c(0.6, 0.9, 0.40),
          lo = c(0.36, 0.70, 0.26),
          hi = c(0.84, 1.10, h2_hi),
          panel = c(
            "variance components",
            "variance components",
            "heritability"
          ),
          level = 0.95,
          interval_method = "asymptotic_reml",
          interval_status = "experimental_asymptotic",
          supplied = FALSE
        )
      )
    ),
    class = "hsquared_fit"
  )
}

mock_gwas <- function() {
  structure(
    data.frame(
      marker = paste0("m", 1:20),
      effect = rnorm(20),
      se = runif(20, 0.1, 0.3),
      z = rnorm(20),
      chisq = rchisq(20, 1),
      p_value = runif(20),
      bonferroni_p = runif(20),
      bh_qvalue = runif(20),
      lod = runif(20, 0, 3)
    ),
    class = c("hs_gwas", "data.frame")
  )
}

test_that("autoplot.hsquared_fit variance returns a ggplot with a zero line", {
  p <- autoplot(mock_uni_fit(), "variance")
  expect_s3_class(p, "ggplot")
  expect_true(any(vapply(
    p$layers,
    function(l) inherits(l$geom, "GeomVline"),
    logical(1)
  )))
})

test_that("variance autoplot consumes the engine variance_components_plot_data", {
  p <- autoplot(mock_uni_fit_vcpd(), "variance")
  expect_s3_class(p, "ggplot")
  # engine terms/panels are used directly (sigma_a2/sigma_e2/h2)
  expect_true(all(
    c("sigma_a2", "sigma_e2", "h2") %in% as.character(p$data$term)
  ))
  expect_equal(p$data$estimate[as.character(p$data$term) == "sigma_a2"], 0.6)
  expect_equal(p$data$lo[as.character(p$data$term) == "h2"], 0.26)
  expect_true(grepl("experimental", p$labels$subtitle))
  expect_equal(
    attr(p, "hsquared_meta")$interval_status,
    "experimental_asymptotic"
  )
})

test_that("variance payload annotates an h2 interval crossing [0,1]", {
  p <- autoplot(mock_uni_fit_vcpd(h2_hi = 1.05), "variance")
  expect_true(grepl("[0,1] boundary", p$labels$subtitle, fixed = TRUE))
})

test_that("variance boundary note is h2-only (a variance whisker crossing is not flagged)", {
  # default mock has sigma_e2 hi = 1.10 (> 1) on the variance panel, h2 in [0,1].
  p <- autoplot(mock_uni_fit_vcpd(h2_hi = 0.54), "variance")
  expect_true(any(p$data$hi[p$data$panel == "variance components"] > 1))
  expect_false(grepl("[0,1] boundary", p$labels$subtitle, fixed = TRUE))
})

test_that("variance payload with only term + estimate degrades to points only", {
  fit <- structure(
    list(
      result = list(
        variance_components_plot_data = list(
          term = c("sigma_a2", "sigma_e2"),
          estimate = c(0.6, 0.9)
        )
      )
    ),
    class = "hsquared_fit"
  )
  p <- autoplot(fit, "variance")
  expect_true(all(is.na(p$data$lo)))
  expect_equal(attr(p, "hsquared_meta")$interval_status, "none")
})

test_that("variance recompute path annotates an h2 CI crossing [0,1]", {
  # no variance_components_plot_data -> the recompute branch; h2 = 0.9, se = 0.30
  # gives a +/- 1.96 SE whisker that crosses 1.
  fit <- mock_uni_fit()
  fit$result$heritability <- data.frame(term = "animal", estimate = 0.9)
  fit$result$heritability_se <- 0.30
  p <- autoplot(fit, "variance")
  h_row <- p$data[p$data$panel == "heritability", ]
  expect_equal(h_row$hi, 0.9 + 1.96 * 0.30, tolerance = 1e-9)
  expect_true(grepl("[0,1] boundary", p$labels$subtitle, fixed = TRUE))
})

test_that("variance payload with NaN intervals draws points only", {
  fit <- mock_uni_fit_vcpd()
  fit$result$variance_components_plot_data$lo <- rep(NaN, 3)
  fit$result$variance_components_plot_data$hi <- rep(NaN, 3)
  fit$result$variance_components_plot_data$interval_status <- "none"
  p <- autoplot(fit, "variance")
  expect_true(all(is.na(p$data$lo)))
  expect_equal(attr(p, "hsquared_meta")$interval_status, "none")
})

test_that("variance payload takes precedence over the recompute fields", {
  fit <- mock_uni_fit() # has variance_components (animal/residual)
  fit$result$variance_components_plot_data <-
    mock_uni_fit_vcpd()$result$variance_components_plot_data
  p <- autoplot(fit, "variance")
  # payload terms win (sigma_a2), not the recompute component name (animal)
  expect_true("sigma_a2" %in% as.character(p$data$term))
  expect_false("animal" %in% as.character(p$data$term))
})

test_that("autoplot.hsquared_fit breeding_values returns a ggplot", {
  p <- autoplot(mock_uni_fit(), "breeding_values")
  expect_s3_class(p, "ggplot")
  # the PEV band is present
  expect_true(any(vapply(
    p$layers,
    function(l) inherits(l$geom, "GeomRibbon"),
    logical(1)
  )))
})

test_that("breeding_values autoplot consumes the engine breeding_values_plot_data", {
  # a fit carrying ONLY the engine payload (no extractable breeding_values), so
  # autoplot must use the payload.
  fit <- structure(
    list(
      result = list(
        breeding_values_plot_data = list(
          id = letters[1:5],
          trait = rep(1L, 5),
          value = c(-0.8, -0.2, 0.1, 0.4, 0.9),
          pev = c(0.20, 0.22, 0.19, 0.21, 0.18),
          pev_scale = "validation"
        )
      )
    ),
    class = "hsquared_fit"
  )
  p <- autoplot(fit, "breeding_values")
  expect_s3_class(p, "ggplot")
  # the plotted values come from the payload (sorted EBV caterpillar)
  expect_setequal(p$data$value, c(-0.8, -0.2, 0.1, 0.4, 0.9))
  # the PEV band is drawn from the payload pev
  expect_true(all(is.finite(p$data$lo)) && all(is.finite(p$data$hi)))
  expect_true(any(vapply(
    p$layers,
    function(l) inherits(l$geom, "GeomRibbon"),
    logical(1)
  )))
})

test_that("hs_breeding_values_from_payload is rename-robust and guards bad input", {
  pl <- list(
    ids = c("a", "b"),
    values = c(1, 2),
    prediction_error_variance = c(0.1, 0.2)
  )
  out <- hsquared:::hs_breeding_values_from_payload(pl)
  expect_equal(out$id, c("a", "b"))
  expect_equal(out$value, c(1, 2))
  expect_equal(out$pev, c(0.1, 0.2))
  # NULL / mismatched payloads fall back (return NULL)
  expect_null(hsquared:::hs_breeding_values_from_payload(NULL))
  expect_null(hsquared:::hs_breeding_values_from_payload(
    list(id = c("a", "b"), value = 1)
  ))
})

test_that("autoplot.hsquared_fit g_matrix returns a ggplot for multivariate", {
  p <- autoplot(mock_mv_fit(), "g_matrix")
  expect_s3_class(p, "ggplot")
})

test_that("g_matrix errors on a univariate fit (no correlation matrix)", {
  expect_error(
    autoplot(mock_uni_fit(), "g_matrix"),
    "multivariate fit",
    fixed = TRUE
  )
})

# A multivariate fit carrying a genetic covariance G (eigen_G reads it).
mock_mv_fit_gcov <- function() {
  g <- matrix(c(2.0, 0.6, 0.6, 0.5), 2, 2)
  structure(
    list(result = list(genetic_covariance = g)),
    class = "hsquared_fit"
  )
}

test_that("g_geometry draws a rotation-invariant eigenvalue scree", {
  p <- autoplot(mock_mv_fit_gcov(), "g_geometry")
  expect_s3_class(p, "ggplot")
  ev <- sort(eigen(matrix(c(2, 0.6, 0.6, 0.5), 2, 2))$values, decreasing = TRUE)
  expect_equal(
    sort(p$data$eigenvalue, decreasing = TRUE),
    ev,
    tolerance = 1e-10
  )
  expect_equal(sum(p$data$variance_explained), 1, tolerance = 1e-10)
  # §3 binding rule: g_geometry must be rotation_invariant; no loadings drawn.
  expect_equal(attr(p, "hsquared_meta")$rotation_status, "rotation_invariant")
  expect_equal(attr(p, "hsquared_meta")$type, "g_geometry")
})

test_that("g_geometry errors on a fit without a genetic covariance", {
  expect_error(
    autoplot(mock_uni_fit(), "g_geometry"),
    "g_geometry",
    fixed = TRUE
  )
})

test_that("g_geometry consumes the engine genetic_pca_plot_data payload", {
  fit <- structure(
    list(
      result = list(
        genetic_pca_plot_data = list(
          eigenvalues = c(2.4, 0.1),
          variance_explained = c(0.96, 0.04),
          axis_labels = c("PC1", "PC2"),
          rotation_invariant = TRUE,
          is_eigenstructure_not_loadings = TRUE
        )
      )
    ),
    class = "hsquared_fit"
  )
  p <- autoplot(fit, "g_geometry")
  expect_equal(p$data$eigenvalue, c(2.4, 0.1))
  expect_true(any(grepl("96%", p$data$label)))
})

test_that("g_geometry ignores a payload that is not rotation-invariant", {
  fit <- mock_mv_fit_gcov()
  fit$result$genetic_pca_plot_data <- list(
    eigenvalues = c(9, 9),
    rotation_invariant = FALSE
  )
  p <- autoplot(fit, "g_geometry")
  # falls back to eigen_G recompute -> not the bogus c(9, 9).
  expect_false(all(p$data$eigenvalue == 9))
})

test_that("g_geometry ignores a payload flagged as loadings (not eigenstructure)", {
  fit <- mock_mv_fit_gcov()
  fit$result$genetic_pca_plot_data <- list(
    eigenvalues = c(9, 9),
    rotation_invariant = TRUE,
    is_eigenstructure_not_loadings = FALSE # §3-enforced: must not be drawn
  )
  p <- autoplot(fit, "g_geometry")
  expect_false(all(p$data$eigenvalue == 9)) # fell back to eigen_G recompute
})

test_that("g_geometry flags a non-PSD payload and omits variance-share labels", {
  fit <- structure(
    list(
      result = list(
        genetic_pca_plot_data = list(
          eigenvalues = c(1.73, -0.23),
          rotation_invariant = TRUE,
          is_eigenstructure_not_loadings = TRUE
        )
      )
    ),
    class = "hsquared_fit"
  )
  p <- autoplot(fit, "g_geometry")
  expect_true(any(p$data$eigenvalue < 0)) # negative bar drawn (honest)
  expect_true(all(p$data$label == "")) # % labels suppressed
  expect_true(grepl("non-positive-definite", p$labels$subtitle))
  expect_match(attr(p, "hsquared_meta")$notes, "non-positive-definite")
})

test_that("g_geometry recomputes variance_explained and labels on length mismatch", {
  fit <- structure(
    list(
      result = list(
        genetic_pca_plot_data = list(
          eigenvalues = c(2, 1, 1),
          variance_explained = c(0.9), # wrong length -> recompute
          axis_labels = c("only one"), # wrong length -> PC labels
          rotation_invariant = TRUE,
          is_eigenstructure_not_loadings = TRUE
        )
      )
    ),
    class = "hsquared_fit"
  )
  p <- autoplot(fit, "g_geometry")
  expect_equal(sum(p$data$variance_explained), 1, tolerance = 1e-10)
  expect_equal(levels(p$data$axis), c("PC1", "PC2", "PC3"))
})

test_that("g_geometry yields NA variance shares for an all-zero eigenstructure", {
  fit <- structure(
    list(
      result = list(
        genetic_pca_plot_data = list(
          eigenvalues = c(0, 0),
          rotation_invariant = TRUE,
          is_eigenstructure_not_loadings = TRUE
        )
      )
    ),
    class = "hsquared_fit"
  )
  p <- autoplot(fit, "g_geometry")
  expect_true(all(is.na(p$data$variance_explained)))
  expect_true(all(p$data$label == ""))
})

test_that("g_geometry payload and recompute agree on identical G", {
  g <- matrix(c(2.0, 0.6, 0.6, 0.5), 2, 2)
  ev <- sort(eigen(g, symmetric = TRUE)$values, decreasing = TRUE)
  fit_recompute <- mock_mv_fit_gcov() # genetic_covariance == g
  fit_payload <- structure(
    list(
      result = list(
        genetic_pca_plot_data = list(
          eigenvalues = ev,
          variance_explained = ev / sum(ev),
          axis_labels = c("PC1", "PC2"),
          rotation_invariant = TRUE,
          is_eigenstructure_not_loadings = TRUE
        )
      )
    ),
    class = "hsquared_fit"
  )
  d_r <- autoplot(fit_recompute, "g_geometry")$data
  d_p <- autoplot(fit_payload, "g_geometry")$data
  expect_equal(d_p$eigenvalue, d_r$eigenvalue, tolerance = 1e-10)
  expect_equal(
    d_p$variance_explained,
    d_r$variance_explained,
    tolerance = 1e-10
  )
  expect_equal(d_p$label, d_r$label)
})

test_that("g_matrix flags off-diagonal cells involving a low-h2 trait", {
  p <- autoplot(mock_mv_fit_h2(c(t1 = 0.05, t2 = 0.4)), "g_matrix")
  expect_s3_class(p, "ggplot")
  off <- p$data[as.character(p$data$row) != as.character(p$data$col), ]
  dia <- p$data[as.character(p$data$row) == as.character(p$data$col), ]
  # t1 h2 = 0.05 < 0.1 -> every off-diagonal cell is imprecise; diagonal is not.
  expect_true(all(off$low))
  expect_false(any(dia$low))
  expect_true(grepl("imprecise", p$labels$subtitle))
})

test_that("g_matrix renders the flag as a glyph, not a literal escape", {
  p <- autoplot(mock_mv_fit_h2(c(t1 = 0.05, t2 = 0.4)), "g_matrix")
  dagger <- intToUtf8(0x2020)
  # the dagger glyph is present on flagged labels and in the subtitle ...
  expect_true(any(grepl(dagger, p$data$label, fixed = TRUE)))
  expect_true(grepl(dagger, p$labels$subtitle, fixed = TRUE))
  # ... and the literal escape text never leaks (guards the backslash-escape bug).
  expect_false(any(grepl("u2020", p$data$label)))
  expect_false(grepl("u2020", p$labels$subtitle))
  expect_false(grepl("u00b2", p$labels$subtitle))
})

test_that("g_matrix does not flag when all traits clear low_h2", {
  p <- autoplot(mock_mv_fit_h2(c(t1 = 0.5, t2 = 0.6)), "g_matrix")
  expect_false(any(p$data$low))
  expect_false(grepl("imprecise", p$labels$subtitle))
})

test_that("g_matrix low_h2 threshold is configurable", {
  fit <- mock_mv_fit_h2(c(t1 = 0.2, t2 = 0.6))
  # the default 0.1 does not flag (0.2 >= 0.1) ...
  expect_false(any(autoplot(fit, "g_matrix")$data$low))
  # ... but a higher threshold does, and the threshold shows in the subtitle.
  p <- autoplot(fit, "g_matrix", low_h2 = 0.3)
  off <- p$data[as.character(p$data$row) != as.character(p$data$col), ]
  expect_true(all(off$low))
  expect_true(grepl("0.3", p$labels$subtitle, fixed = TRUE))
})

test_that("g_matrix handles a single NA heritability without false flags", {
  # t1 NA, t2 clears -> the off-diagonal involves no finite low-h2 trait: no flag.
  p_na <- autoplot(mock_mv_fit_h2(c(t1 = NA, t2 = 0.4)), "g_matrix")
  expect_false(any(p_na$data$low))
  # t1 low, t2 NA -> the off-diagonal involves t1 (finite, low): flagged.
  p_low <- autoplot(mock_mv_fit_h2(c(t1 = 0.05, t2 = NA)), "g_matrix")
  off <- p_low$data[
    as.character(p_low$data$row) != as.character(p_low$data$col),
  ]
  expect_true(all(off$low))
})

test_that("g_matrix degrades gracefully without heritabilities (no flag)", {
  p <- autoplot(mock_mv_fit(), "g_matrix")
  expect_s3_class(p, "ggplot")
  expect_false(any(p$data$low))
  # recompute path uses the genetic_correlation field value (0.3).
  off <- p$data[as.character(p$data$row) != as.character(p$data$col), ]
  expect_true(all(abs(off$value - 0.3) < 1e-9))
})

test_that("g_matrix consumes the engine genetic_correlation_plot_data payload", {
  p <- autoplot(mock_mv_fit_plotdata(), "g_matrix")
  off <- p$data[as.character(p$data$row) != as.character(p$data$col), ]
  # payload correlation 0.55 is used, not the 0.3 in result$genetic_correlation
  expect_true(all(abs(off$value - 0.55) < 1e-9))
  # payload heritabilities (t1 = 0.05) flag the imprecise off-diagonal cells
  expect_true(all(off$low))
})

test_that("g_matrix ignores a payload that is not rotation-invariant", {
  fit <- mock_mv_fit_plotdata()
  fit$result$genetic_correlation_plot_data$rotation_invariant <- FALSE
  p <- autoplot(fit, "g_matrix")
  off <- p$data[as.character(p$data$row) != as.character(p$data$col), ]
  # falls back to the recompute path -> the 0.3 genetic_correlation, not 0.55.
  expect_true(all(abs(off$value - 0.3) < 1e-9))
})

test_that("g_matrix drops mismatched-length payload heritabilities (no flag)", {
  fit <- mock_mv_fit_plotdata()
  # 3 heritabilities for a 2-trait G -> the guard drops h2 rather than mis-flag.
  fit$result$genetic_correlation_plot_data$heritabilities <- c(0.02, 0.5, 0.5)
  p <- autoplot(fit, "g_matrix")
  expect_false(any(p$data$low))
})

test_that("g_matrix handles a payload with NULL traits (engine default labels)", {
  fit <- mock_mv_fit_plotdata()
  fit$result$genetic_correlation_plot_data$traits <- NULL
  fit$result$genetic_correlation_plot_data$genetic_correlations <-
    unname(fit$result$genetic_correlation_plot_data$genetic_correlations)
  p <- autoplot(fit, "g_matrix")
  expect_s3_class(p, "ggplot")
  # default labels match the engine preparer convention (trait_1, trait_2).
  expect_true(all(c("trait_1", "trait_2") %in% as.character(p$data$row)))
})

test_that("g_matrix payload and recompute paths agree on identical G", {
  # Julia-free consumer parity: the same correlation via the payload vs the
  # recompute path must yield identical tidy frames (value, low, label).
  corr <- matrix(c(1, 0.42, 0.42, 1), 2, 2)
  h2 <- c(t1 = 0.05, t2 = 0.4)
  fit_recompute <- mock_mv_fit_h2(h2)
  fit_recompute$result$genetic_correlation <- matrix(
    c(1, 0.42, 0.42, 1),
    2,
    2,
    dimnames = list(c("t1", "t2"), c("t1", "t2"))
  )
  fit_payload <- structure(
    list(
      result = list(
        genetic_correlation_plot_data = list(
          traits = c("t1", "t2"),
          genetic_correlations = corr,
          heritabilities = as.numeric(h2),
          rotation_invariant = TRUE
        )
      )
    ),
    class = "hsquared_fit"
  )
  d_recompute <- autoplot(fit_recompute, "g_matrix")$data
  d_payload <- autoplot(fit_payload, "g_matrix")$data
  expect_equal(d_payload$value, d_recompute$value)
  expect_equal(d_payload$low, d_recompute$low)
  expect_equal(d_payload$label, d_recompute$label)
})

test_that("reaction_norm consumes the engine rr_genetic_variance_plot_data payload", {
  # payload-only fit: if the code tried to recompute it would error (not an RR
  # fit), so a successful figure using the payload values proves consumption.
  pd <- list(
    covariate = c(1, 2, 3),
    value = c(0.5, 0.8, 0.6),
    heritability = c(0.2, 0.3, 0.25)
  )
  fit <- structure(
    list(result = list(rr_genetic_variance_plot_data = pd)),
    class = "hsquared_fit"
  )
  p <- autoplot(fit, "reaction_norm")
  expect_s3_class(p, "ggplot")
  gv <- p$data[p$data$panel == "genetic variance", ]
  expect_equal(gv$value, c(0.5, 0.8, 0.6))
  expect_equal(attr(p, "hsquared_meta")$rotation_status, "rotation_invariant")
})

test_that("reaction_norm payload consumption is rename-robust", {
  base <- list(covariate = c(1, 2, 3), heritability = c(0.2, 0.3, 0.25))
  fit_v <- structure(
    list(
      result = list(
        rr_genetic_variance_plot_data = c(base, list(value = c(0.5, 0.8, 0.6)))
      )
    ),
    class = "hsquared_fit"
  )
  fit_g <- structure(
    list(
      result = list(
        rr_genetic_variance_plot_data = c(
          base,
          list(genetic_variance = c(0.5, 0.8, 0.6))
        )
      )
    ),
    class = "hsquared_fit"
  )
  gv_v <- autoplot(fit_v, "reaction_norm")$data
  gv_g <- autoplot(fit_g, "reaction_norm")$data
  expect_equal(gv_v$value, gv_g$value)
  expect_equal(
    gv_v$value[gv_v$panel == "genetic variance"],
    c(0.5, 0.8, 0.6)
  )
})

test_that("autoplot.hs_gwas returns a Manhattan ggplot with a Bonferroni line", {
  p <- autoplot(mock_gwas())
  expect_s3_class(p, "ggplot")
  expect_true(any(vapply(
    p$layers,
    function(l) inherits(l$geom, "GeomHline"),
    logical(1)
  )))
})

test_that("rr_eigenfunctions autoplot consumes the engine payload", {
  ef_mat <- matrix(c(0.5, 0.6, 0.7, -0.1, 0.0, 0.1), nrow = 3) # 3 cov x 2 axes
  fit <- structure(
    list(
      result = list(
        rr_eigenfunctions_plot_data = list(
          covariate = c(1, 2, 3),
          eigenfunctions = ef_mat,
          variance_explained = c(0.8, 0.2),
          rotation_invariant = TRUE
        )
      )
    ),
    class = "hsquared_fit"
  )
  p <- autoplot(fit, "rr_eigenfunctions")
  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), 6L) # 3 covariate x 2 axes
  expect_equal(sort(p$data$value), sort(as.numeric(ef_mat)))
  expect_equal(attr(p, "hsquared_meta")$type, "rr_eigenfunctions")
  expect_equal(attr(p, "hsquared_meta")$rotation_status, "rotation_invariant")
  expect_true(any(grepl("80%", levels(p$data$axis_label), fixed = TRUE)))
})

test_that("rr_surface autoplot consumes the engine payload", {
  surf <- matrix(c(0.4, 0.2, 0.1, 0.2, 0.5, 0.2, 0.1, 0.2, 0.6), 3, 3)
  fit <- structure(
    list(
      result = list(
        rr_covariance_surface_plot_data = list(
          covariate = c(1, 2, 3),
          surface = surf,
          is_correlation = FALSE,
          rotation_invariant = TRUE
        )
      )
    ),
    class = "hsquared_fit"
  )
  p <- autoplot(fit, "rr_surface")
  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), 9L) # 3 x 3 grid
  expect_equal(sort(p$data$value), sort(as.numeric(surf)))
  expect_equal(attr(p, "hsquared_meta")$type, "rr_surface")
  expect_equal(attr(p, "hsquared_meta")$rotation_status, "rotation_invariant")
  expect_true(grepl("covariance", p$labels$title))
})

test_that("autoplot.hs_gwas qq returns a ggplot with a y=x null and lambda_GC", {
  set.seed(1)
  p <- autoplot(mock_gwas(), "qq")
  expect_s3_class(p, "ggplot")
  # the y = x null reference line is present
  expect_true(any(vapply(
    p$layers,
    function(l) inherits(l$geom, "GeomAbline"),
    logical(1)
  )))
  # QQ data: observed/expected -log10(p), one point per marker, sorted-aligned
  expect_equal(nrow(p$data), 20L)
  expect_true(all(c("expected", "observed") %in% names(p$data)))
  expect_false(is.unsorted(rev(p$data$expected))) # expected descending
  m <- attr(p, "hsquared_meta")
  expect_equal(m$type, "qq")
  expect_equal(m$interval_status, "uncalibrated")
  expect_true(grepl("lambda_GC", p$labels$subtitle))
})

test_that("autoplot.hs_gwas notes a relatedness-uncorrected single-marker scan", {
  g <- mock_gwas()
  attr(g, "scan_method") <- "single"
  expect_true(grepl("UNcorrected", autoplot(g, "manhattan")$labels$subtitle))
  expect_true(grepl("UNcorrected", autoplot(g, "qq")$labels$subtitle))
  # the default (mixed) scan carries no such note
  expect_false(grepl("UNcorrected", autoplot(mock_gwas())$labels$subtitle))
})

test_that("autoplot.hs_gwas notes a LOCO genomic scan", {
  g <- mock_gwas()
  attr(g, "scan_method") <- "loco"
  # the note carries the load-bearing honesty half (the pedigree-VC caveat),
  # not just the bare word "LOCO"
  expect_true(grepl(
    "LOCO genomic correction \\(pedigree-estimated VCs\\)",
    autoplot(g, "manhattan")$labels$subtitle
  ))
  expect_true(grepl(
    "pedigree-estimated VCs",
    autoplot(g, "qq")$labels$subtitle
  ))
  # the default (mixed) scan carries no LOCO note
  expect_false(grepl("LOCO", autoplot(mock_gwas())$labels$subtitle))
})

test_that("autoplot.hs_gwas defaults to the Manhattan", {
  expect_equal(attr(autoplot(mock_gwas()), "hsquared_meta")$type, "manhattan")
  expect_equal(
    attr(autoplot(mock_gwas(), "manhattan"), "hsquared_meta")$type,
    "manhattan"
  )
})

test_that("hs_recovery_forest returns a ggplot and validates columns", {
  rec <- data.frame(
    target = c("a", "b"),
    bias = c(0.01, -0.02),
    mcse = c(0.02, 0.03)
  )
  p <- hs_recovery_forest(rec)
  expect_s3_class(p, "ggplot")
  expect_error(
    hs_recovery_forest(data.frame(x = 1)),
    "must have",
    fixed = TRUE
  )
})

test_that("variance autoplot errors when there are no variance components", {
  bad <- structure(list(result = list()), class = "hsquared_fit")
  expect_error(autoplot(bad, "variance"), "no variance-component")
})

test_that("autoplot figures carry the hsquared_meta honest-status attribute", {
  m_g <- attr(autoplot(mock_mv_fit(), "g_matrix"), "hsquared_meta")
  expect_equal(m_g$type, "g_matrix")
  expect_equal(m_g$rotation_status, "rotation_invariant")

  m_gw <- attr(autoplot(mock_gwas()), "hsquared_meta")
  expect_equal(m_gw$type, "manhattan")
  expect_equal(m_gw$interval_status, "uncalibrated")

  m_v <- attr(autoplot(mock_uni_fit(), "variance"), "hsquared_meta")
  expect_equal(m_v$type, "variance")
  expect_equal(m_v$interval_status, "experimental_asymptotic")
})
