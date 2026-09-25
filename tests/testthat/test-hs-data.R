test_that("hs_data stores phenotype, pedigree, and genotype ID maps", {
  phenotypes <- data.frame(id = c("a", "b", "b"), y = c(1, 2, 3))
  pedigree <- data.frame(
    id = c("a", "b", "c"),
    sire = c(NA, NA, "a"),
    dam = c(NA, NA, "b")
  )
  genotypes <- matrix(
    1,
    nrow = 2,
    ncol = 3,
    dimnames = list(c("b", "d"), c("m1", "m2", "m3"))
  )

  data <- hs_data(
    phenotypes = phenotypes,
    pedigree = pedigree,
    genotypes = genotypes,
    markers = data.frame(
      marker = c("m1", "m2", "m3"),
      chr = c("1", "1", "2"),
      pos = c(10, 20, 5)
    )
  )

  expect_s3_class(data, "hs_data")
  expect_s3_class(data$marker_spec, "hs_marker_map_spec")
  expect_s3_class(data$genotype_marker_spec, "hs_genotype_marker_spec")
  expect_equal(data$id_map$phenotype_ids, c("a", "b"))
  expect_equal(data$id_map$pedigree_ids, c("a", "b", "c"))
  expect_equal(data$id_map$genotype_ids, c("b", "d"))
  expect_equal(data$marker_spec$marker_ids, c("m1", "m2", "m3"))
  expect_equal(
    data$marker_spec$columns,
    list(
      marker = 1L,
      chromosome = 2L,
      position = 3L
    )
  )
  expect_equal(data$genotype_marker_spec$marker_ids, c("m1", "m2", "m3"))
  expect_equal(data$genotype_marker_spec$marker_map_index, c(1L, 2L, 3L))
  expect_equal(data$id_map$phenotypes_without_pedigree, character())
  expect_equal(data$id_map$phenotypes_without_genotypes, "a")
  expect_equal(data$id_map$genotypes_without_phenotypes, "d")
  expect_s3_class(summary(data), "summary_hs_data")

  overlap <- summary(data)$id_overlap
  expect_equal(
    overlap$metric,
    c(
      "phenotype_ids",
      "pedigree_ids",
      "genotype_ids",
      "expression_ids",
      "phenotypes_without_pedigree",
      "phenotypes_without_genotypes",
      "genotypes_without_phenotypes",
      "phenotypes_without_expression",
      "expression_without_phenotypes"
    )
  )
  expect_equal(
    overlap$count,
    c(2L, 3L, 2L, 0L, 0L, 1L, 1L, 2L, 0L)
  )

  genotype_status <- summary(data)$genotype_status
  expect_equal(
    genotype_status$metric,
    c(
      "genotype_rows",
      "genotype_ids",
      "genotype_marker_columns",
      "named_genotype_marker_columns",
      "unnamed_genotype_marker_columns",
      "duplicate_genotype_marker_columns",
      "missing_genotype_values",
      "component_type"
    )
  )
  expect_equal(
    genotype_status$value,
    c("2", "2", "3", "3", "0", "0", "0", "matrix")
  )

  pedigree_status <- summary(data)$pedigree_status
  expect_equal(
    pedigree_status$metric,
    c(
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
    )
  )
  expect_equal(
    pedigree_status$count,
    c(3L, 3L, 2L, 1L, 2L, 1L, 1L, 1L, 0L, 0L, 0L, 0L)
  )

  marker_status <- summary(data)$marker_status
  expect_equal(
    marker_status$metric,
    c(
      "marker_map_markers",
      "genotype_marker_columns",
      "aligned_marker_columns",
      "chromosomes",
      "position_min",
      "position_max",
      "alignment"
    )
  )
  expect_equal(
    marker_status$value,
    c("3", "3", "3", "2", "5", "20", "checked")
  )
})

test_that("summary.hs_data reports partial marker diagnostics", {
  phenotypes <- data.frame(id = "a", y = 1)

  marker_only <- hs_data(
    phenotypes = phenotypes,
    markers = data.frame(marker = "m1", chr = "1", pos = 10)
  )
  expect_equal(
    summary(marker_only)$marker_status$value,
    c("1", "0", "0", "1", "10", "10", "not_checked_no_genotypes")
  )

  genotype_only <- hs_data(
    phenotypes = phenotypes,
    genotypes = data.frame(id = "a", m1 = 0, m2 = 1)
  )
  expect_equal(
    summary(genotype_only)$genotype_status$value,
    c("1", "1", "2", "2", "0", "0", "0", "data.frame")
  )
  expect_equal(
    summary(genotype_only)$marker_status$value,
    c(
      "0",
      "2",
      "0",
      "not_available",
      "not_available",
      "not_available",
      "not_checked_no_marker_map"
    )
  )

  unnamed_matrix <- hs_data(
    phenotypes = phenotypes,
    genotypes = matrix(
      0,
      nrow = 1,
      ncol = 2,
      dimnames = list("a", NULL)
    )
  )
  expect_equal(summary(unnamed_matrix)$marker_status$value[2], "2")
  expect_equal(
    summary(unnamed_matrix)$genotype_status$value,
    c("1", "1", "2", "0", "2", "0", "0", "matrix")
  )
  expect_equal(
    summary(unnamed_matrix)$marker_status$value[7],
    "not_checked_no_marker_map"
  )

  phenotype_only <- hs_data(phenotypes = phenotypes)
  expect_null(summary(phenotype_only)$pedigree_status)
  expect_null(summary(phenotype_only)$genotype_status)
  expect_null(summary(phenotype_only)$marker_status)
})

test_that("hs_data reports duplicate and missing genotype diagnostics", {
  phenotypes <- data.frame(id = c("a", "b"), y = c(1, 2))
  genotypes <- data.frame(
    id = c("a", "b"),
    m1 = c(0, NA),
    m1 = c(1, 2),
    check.names = FALSE
  )

  data <- hs_data(
    phenotypes = phenotypes,
    genotypes = genotypes
  )

  expect_equal(
    summary(data)$genotype_status$value,
    c("2", "2", "2", "2", "0", "1", "1", "data.frame")
  )
})

test_that("summary.hs_data reports pedigree warning diagnostics", {
  data <- hs_data(
    phenotypes = data.frame(id = "a", y = 1),
    pedigree = data.frame(
      id = c("a", "b", "b", "c"),
      sire = c(NA, "ghost", NA, "a"),
      dam = c(NA, "a", NA, "a")
    )
  )

  status <- summary(data)$pedigree_status
  expect_equal(status$count[status$metric == "pedigree_rows"], 4L)
  expect_equal(status$count[status$metric == "pedigree_ids"], 3L)
  expect_equal(status$count[status$metric == "missing_known_parent_ids"], 1L)
  expect_equal(status$count[status$metric == "duplicate_pedigree_ids"], 1L)
  expect_equal(status$count[status$metric == "same_known_parent_rows"], 1L)
})

test_that("data_status exposes hs_data diagnostics directly", {
  phenotypes <- data.frame(id = c("a", "b"), y = c(1, 2))
  pedigree <- data.frame(
    id = c("a", "b"),
    sire = c(NA, NA),
    dam = c(NA, NA)
  )
  genotypes <- matrix(
    0,
    nrow = 2,
    ncol = 2,
    dimnames = list(c("a", "b"), c("m1", "m2"))
  )
  data <- hs_data(
    phenotypes = phenotypes,
    pedigree = pedigree,
    genotypes = genotypes,
    markers = data.frame(
      marker = c("m1", "m2"),
      chr = c("1", "2"),
      pos = c(10, 20)
    )
  )

  status <- data_status(data)

  expect_s3_class(status, "hs_data_status")
  expect_equal(
    status$components,
    c("phenotypes", "pedigree", "genotypes", "markers")
  )
  expect_equal(status$id_overlap$count[[1L]], 2L)
  expect_equal(status$pedigree_status$count[[5L]], 2L)
  expect_equal(status$genotype_status$value[[3L]], "2")
  expect_equal(status$marker_status$value[[7L]], "checked")
  expect_match(capture.output(print(status))[[1L]], "<hs_data_status>")
})

test_that("hs_data reports environment key diagnostics", {
  phenotypes <- data.frame(
    id = c("a", "b", "c"),
    env = c("E1", "E1", "E3"),
    y = c(1, 2, 3)
  )
  environment <- data.frame(
    env = c("E1", "E2", "E2"),
    temperature = c(18, 20, 21)
  )

  data <- hs_data(
    phenotypes = phenotypes,
    environment = environment,
    environment_id = "env"
  )

  expect_s3_class(data$environment_spec, "hs_environment_spec")
  expect_equal(data$environment_spec$phenotype_environment_ids, c("E1", "E3"))
  expect_equal(data$environment_spec$environment_ids, c("E1", "E2"))
  expect_equal(data$environment_spec$phenotypes_without_environment, "E3")
  expect_equal(data$environment_spec$environment_without_phenotypes, "E2")
  expect_equal(data$environment_spec$duplicate_environment_ids, "E2")

  environment_status <- summary(data)$environment_status
  expect_equal(
    environment_status$metric,
    c(
      "environment_rows",
      "environment_key",
      "environment_ids",
      "phenotype_environment_ids",
      "phenotype_environment_ids_with_metadata",
      "environment_only_ids",
      "phenotype_environment_ids_without_metadata",
      "duplicate_environment_ids"
    )
  )
  expect_equal(
    environment_status$value,
    c("3", "env", "2", "2", "1", "1", "1", "1")
  )

  status <- data_status(data)
  expect_equal(status$environment_status$value[[2L]], "env")
  expect_equal(status$environment_status$value[[7L]], "1")
})

test_that("hs_data reports unkeyed environment tables without overchecking", {
  data <- hs_data(
    phenotypes = data.frame(id = "a", y = 1),
    environment = data.frame(env = "E1", rainfall = 4)
  )

  environment_status <- summary(data)$environment_status
  expect_equal(environment_status$value[[1L]], "1")
  expect_equal(
    environment_status$value[[2L]],
    "not_checked_no_environment_id"
  )
  expect_equal(environment_status$value[[3L]], "not_available")
})

test_that("hs_data accepts expression IDs from an ID column", {
  phenotypes <- data.frame(sample = c("s1", "s2"), y = c(1, 2))
  expression <- data.frame(sample = c("s2", "s3"), gene1 = c(10, 20))

  data <- hs_data(
    phenotypes = phenotypes,
    expression = expression,
    id = "sample"
  )

  expect_equal(data$id_map$expression_ids, c("s2", "s3"))
  expect_equal(data$id_map$phenotypes_without_expression, "s1")
  expect_equal(data$id_map$expression_without_phenotypes, "s3")
  expect_equal(summary(data)$expression_status$value[[7L]], "data.frame")
})

test_that("hs_data reports expression feature diagnostics", {
  phenotypes <- data.frame(id = c("s1", "s2"), y = c(1, 2))
  expression <- matrix(
    1:6,
    nrow = 2,
    dimnames = list(c("s1", "s2"), c("gene1", "gene2", "gene2"))
  )

  data <- hs_data(
    phenotypes = phenotypes,
    expression = expression
  )

  expression_status <- summary(data)$expression_status
  expect_equal(
    expression_status$metric,
    c(
      "expression_rows",
      "expression_ids",
      "expression_features",
      "named_expression_features",
      "unnamed_expression_features",
      "duplicate_expression_features",
      "component_type"
    )
  )
  expect_equal(
    expression_status$value,
    c("2", "2", "3", "3", "0", "1", "matrix")
  )

  status <- data_status(data)
  expect_equal(status$expression_status$value[[3L]], "3")
  expect_equal(status$expression_status$value[[6L]], "1")
})

test_that("hs_data reports unnamed expression matrix features", {
  data <- hs_data(
    phenotypes = data.frame(id = "s1", y = 1),
    expression = matrix(
      1,
      nrow = 1,
      ncol = 2,
      dimnames = list("s1", NULL)
    )
  )

  expression_status <- summary(data)$expression_status
  expect_equal(
    expression_status$value,
    c("1", "1", "2", "0", "2", "0", "matrix")
  )

  partly_named <- hs_data(
    phenotypes = data.frame(id = "s1", y = 1),
    expression = matrix(
      1,
      nrow = 1,
      ncol = 2,
      dimnames = list("s1", c("gene1", NA))
    )
  )
  expect_equal(
    summary(partly_named)$expression_status$value,
    c("1", "1", "2", "1", "1", "0", "matrix")
  )
})

test_that("hs_data reports expression annotation diagnostics", {
  phenotypes <- data.frame(id = c("s1", "s2"), y = c(1, 2))
  expression <- matrix(
    c(1, 2, 3, 4),
    nrow = 2,
    dimnames = list(c("s1", "s2"), c("gene1", "gene3"))
  )
  annotation <- data.frame(
    gene_id = c("gene1", "gene2", "gene2"),
    chromosome = c("1", "1", "2")
  )

  data <- hs_data(
    phenotypes = phenotypes,
    expression = expression,
    annotation = annotation,
    annotation_id = "gene_id"
  )

  expect_s3_class(data$annotation_spec, "hs_annotation_spec")
  expect_equal(data$annotation_spec$annotation_features, c("gene1", "gene2"))
  expect_equal(data$annotation_spec$expression_features, c("gene1", "gene3"))
  expect_equal(data$annotation_spec$expression_without_annotation, "gene3")
  expect_equal(data$annotation_spec$annotation_without_expression, "gene2")
  expect_equal(data$annotation_spec$duplicate_annotation_features, "gene2")

  annotation_status <- summary(data)$annotation_status
  expect_equal(
    annotation_status$metric,
    c(
      "annotation_rows",
      "annotation_key",
      "annotation_features",
      "expression_features",
      "expression_features_with_annotation",
      "annotation_only_features",
      "expression_features_without_annotation",
      "duplicate_annotation_features"
    )
  )
  expect_equal(
    annotation_status$value,
    c("3", "gene_id", "2", "2", "1", "1", "1", "1")
  )

  status <- data_status(data)
  expect_equal(status$annotation_status$value[[2L]], "gene_id")
  expect_equal(status$annotation_status$value[[7L]], "1")
})

test_that("hs_data reports unkeyed annotation tables without overchecking", {
  data <- hs_data(
    phenotypes = data.frame(id = "s1", y = 1),
    annotation = data.frame(gene = "gene1", chr = "1")
  )

  annotation_status <- summary(data)$annotation_status
  expect_equal(annotation_status$value[[1L]], "1")
  expect_equal(
    annotation_status$value[[2L]],
    "not_checked_no_annotation_id"
  )
  expect_equal(annotation_status$value[[3L]], "not_available")
})

test_that("hs_data rejects missing phenotype IDs and incomplete pedigree IDs", {
  expect_error(
    hs_data(data.frame(id = c("a", NA), y = c(1, 2))),
    "`phenotypes` IDs cannot be missing",
    fixed = TRUE
  )

  expect_error(
    hs_data(
      phenotypes = data.frame(id = c("a", "b"), y = c(1, 2)),
      pedigree = data.frame(id = "a", sire = NA, dam = NA)
    ),
    "not present in `pedigree`: b",
    fixed = TRUE
  )
})

test_that("hs_data rejects unsupported component shapes", {
  expect_error(
    hs_data(
      phenotypes = data.frame(id = "a", y = 1),
      genotypes = matrix(1, nrow = 1)
    ),
    "matrix must have individual IDs as row names",
    fixed = TRUE
  )

  expect_error(
    hs_data(
      phenotypes = data.frame(id = "a", y = 1),
      genotypes = data.frame(m1 = 1)
    ),
    "must contain ID column `id` or non-missing row names",
    fixed = TRUE
  )

  expect_error(
    hs_data(
      phenotypes = data.frame(id = "a", y = 1),
      markers = matrix(1)
    ),
    "`markers` must be a data frame",
    fixed = TRUE
  )

  expect_error(
    hs_data(
      phenotypes = data.frame(id = "a", y = 1),
      annotation = matrix(1)
    ),
    "`annotation` must be a data frame",
    fixed = TRUE
  )

  expect_error(
    hs_data(
      phenotypes = data.frame(id = "a", y = 1),
      environment = matrix(1)
    ),
    "`environment` must be a data frame",
    fixed = TRUE
  )
})

test_that("hs_data validates annotation key inputs", {
  expect_error(
    hs_data(
      phenotypes = data.frame(id = "s1", y = 1),
      annotation_id = "gene_id"
    ),
    "`annotation_id` can be supplied only when `annotation` is supplied",
    fixed = TRUE
  )

  expect_error(
    hs_data(
      phenotypes = data.frame(id = "s1", y = 1),
      annotation = data.frame(gene_id = "gene1"),
      annotation_id = ""
    ),
    "`annotation_id` must be one non-empty column name",
    fixed = TRUE
  )

  expect_error(
    hs_data(
      phenotypes = data.frame(id = "s1", y = 1),
      annotation = data.frame(feature = "gene1"),
      annotation_id = "gene_id"
    ),
    "`annotation_id` column `gene_id` was not found in `annotation`",
    fixed = TRUE
  )

  expect_error(
    hs_data(
      phenotypes = data.frame(id = "s1", y = 1),
      expression = matrix(1, nrow = 1, dimnames = list("s1", NULL)),
      annotation = data.frame(gene_id = "gene1"),
      annotation_id = "gene_id"
    ),
    "`expression` matrix must have feature IDs as column names",
    fixed = TRUE
  )

  expect_error(
    hs_data(
      phenotypes = data.frame(id = "s1", y = 1),
      expression = data.frame(id = "s1"),
      annotation = data.frame(gene_id = "gene1"),
      annotation_id = "gene_id"
    ),
    "`expression` must contain at least one feature column",
    fixed = TRUE
  )

  expect_error(
    hs_data(
      phenotypes = data.frame(id = "s1", y = 1),
      annotation = data.frame(gene_id = NA),
      annotation_id = "gene_id"
    ),
    "`annotation` column `gene_id` cannot contain missing or empty values",
    fixed = TRUE
  )
})

test_that("hs_data validates environment key inputs", {
  expect_error(
    hs_data(
      phenotypes = data.frame(id = "a", y = 1),
      environment_id = "env"
    ),
    "`environment_id` can be supplied only when `environment` is supplied",
    fixed = TRUE
  )

  expect_error(
    hs_data(
      phenotypes = data.frame(id = "a", y = 1),
      environment = data.frame(env = "E1"),
      environment_id = ""
    ),
    "`environment_id` must be one non-empty column name",
    fixed = TRUE
  )

  expect_error(
    hs_data(
      phenotypes = data.frame(id = "a", y = 1),
      environment = data.frame(env = "E1"),
      environment_id = "env"
    ),
    "`environment_id` column `env` was not found in `phenotypes`",
    fixed = TRUE
  )

  expect_error(
    hs_data(
      phenotypes = data.frame(id = "a", env = "E1", y = 1),
      environment = data.frame(site = "E1"),
      environment_id = "env"
    ),
    "`environment_id` column `env` was not found in `environment`",
    fixed = TRUE
  )

  expect_error(
    hs_data(
      phenotypes = data.frame(id = "a", env = NA, y = 1),
      environment = data.frame(env = "E1"),
      environment_id = "env"
    ),
    "`phenotypes` column `env` cannot contain missing or empty values",
    fixed = TRUE
  )
})

test_that("hs_data validates marker-map columns and marker positions", {
  phenotypes <- data.frame(id = "a", y = 1)

  expect_error(
    hs_data(
      phenotypes = phenotypes,
      markers = data.frame(marker = "m1", chr = "1")
    ),
    "Missing: position.",
    fixed = TRUE
  )

  expect_error(
    hs_data(
      phenotypes = phenotypes,
      markers = data.frame(
        snp = c("m1", "m1"),
        chromosome = c("1", "1"),
        bp = c(10, 20)
      )
    ),
    "`markers` contains duplicate marker IDs.",
    fixed = TRUE
  )

  expect_error(
    hs_data(
      phenotypes = phenotypes,
      markers = data.frame(
        id = "m1",
        chrom = "1",
        base_pair = -1
      )
    ),
    "finite non-negative numeric positions",
    fixed = TRUE
  )

  expect_error(
    hs_data(
      phenotypes = phenotypes,
      markers = data.frame(
        id = "m1",
        chr = "",
        pos = 1
      )
    ),
    "chromosome column cannot contain missing or empty values",
    fixed = TRUE
  )
})

test_that("hs_data validates genotype marker columns against marker maps", {
  phenotypes <- data.frame(id = "a", y = 1)

  expect_error(
    hs_data(
      phenotypes = phenotypes,
      genotypes = matrix(
        0,
        nrow = 1,
        ncol = 2,
        dimnames = list("a", c("m1", "m_extra"))
      ),
      markers = data.frame(
        marker = c("m1", "m2"),
        chr = c("1", "1"),
        pos = c(10, 20)
      )
    ),
    "missing from `markers`: m_extra; missing from `genotypes`: m2",
    fixed = TRUE
  )

  expect_error(
    hs_data(
      phenotypes = phenotypes,
      genotypes = matrix(
        0,
        nrow = 1,
        ncol = 1,
        dimnames = list("a", NULL)
      ),
      markers = data.frame(marker = "m1", chr = "1", pos = 10)
    ),
    "must have marker IDs as column names",
    fixed = TRUE
  )

  expect_error(
    hs_data(
      phenotypes = phenotypes,
      genotypes = data.frame(id = "a"),
      markers = data.frame(marker = "m1", chr = "1", pos = 10)
    ),
    "at least one marker column",
    fixed = TRUE
  )

  data <- hs_data(
    phenotypes = phenotypes,
    genotypes = data.frame(id = "a", m2 = 1, m1 = 0),
    markers = data.frame(
      marker = c("m1", "m2"),
      chr = c("1", "1"),
      pos = c(10, 20)
    )
  )

  expect_equal(data$genotype_marker_spec$marker_ids, c("m2", "m1"))
  expect_equal(data$genotype_marker_spec$marker_map_index, c(2L, 1L))
})

test_that("hsquared can validate the v0.1 formula from an hs_data bundle", {
  pedigree <- data.frame(
    id = c("a", "b"),
    sire = c(NA, NA),
    dam = c(NA, NA)
  )
  phenotypes <- data.frame(y = c(1, 2), id = c("a", "b"))
  bundle <- hs_data(
    phenotypes = phenotypes,
    pedigree = pedigree
  )

  expect_message(
    spec <- hsquared(
      y ~ animal(1 | id, pedigree = pedigree),
      data = bundle,
      control = hs_control(engine = "validate")
    ),
    "Validated the v0.1 animal-model contract",
    fixed = TRUE
  )
  expect_match(spec$bridge$target, "fit_animal_model", fixed = TRUE)
})
