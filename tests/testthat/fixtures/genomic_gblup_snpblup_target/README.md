# Genomic GBLUP / SNP-BLUP target fixture (#49)

Julia-native supplied-variance genomic target for the #49 external-comparator
gate. The fixture pins a small VanRaden method-1 genomic relationship with
supplied allele frequencies, so `G` is positive definite and the precision-route
GBLUP (`fit_gblup`) can be compared directly with the equivalent marker-route
SNP-BLUP (`fit_snp_blup`) without ridge regularization.

This is a serialized Julia target, not external comparator evidence. It exists
so the R lane or an external executable/package can fit the same estimand and
record a documented tolerance later.

## Model

- 4 individuals (`g1`-`g4`), 6 markers (`m1`-`m6`).
- Fixed effect: intercept only.
- Supplied variance components: `sigma_g2 = 2`, `sigma_e2 = 1`.
- VanRaden method 1 with supplied allele frequencies in
  `allele_frequencies.csv`.
- `G` is positive definite; `Ginv = inv(G)` is serialized for precision-route
  comparators.

## Files

- `phenotypes.csv` - individual IDs and response values.
- `markers.csv` - genotype dosages, rows aligned to `phenotypes.csv`.
- `allele_frequencies.csv` - supplied marker allele frequencies.
- `expected_genomic_relationship.csv` - target `G`.
- `expected_genomic_precision.csv` - target `Ginv`.
- `expected_beta.csv` - intercept estimate.
- `expected_gebv.csv` - GBLUP and SNP-BLUP genomic breeding values.
- `expected_marker_effects.csv` - SNP-BLUP marker effects.
- `expected_metadata.csv` - variance components, VanRaden scale, and route
  agreement metadata.
- `generate.jl` - reproducible generator.

## Regenerate

```sh
julia --project=. test/fixtures/genomic_gblup_snpblup_target/generate.jl
```

## Comparator protocol

An external comparator should use the same `y`, intercept-only `X`, marker matrix,
supplied allele frequencies, and supplied variance components. Precision-route
comparators can consume `expected_genomic_precision.csv`; marker-effect
comparators can consume `markers.csv` plus `allele_frequencies.csv`.

Record package/executable versions, input translation details, fitted beta,
GEBVs, marker effects when available, and the max absolute deviation from the
stored targets. Do not treat this fixture as a covered validation row until an
independent same-estimand comparator run is recorded.
