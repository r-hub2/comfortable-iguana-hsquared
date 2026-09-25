# Phase 4 Multi-Trait Parity Fixture

This directory stores a deterministic two-trait animal-model fixture for future
R-lane sommer/ASReml/BLUPF90 parity checks.

It is not a published external comparator and does not promote any validation
row to covered status. The expected values are Julia REML targets from
`fit_multivariate_reml` on this fixture, recorded so another engine can compare
against the same data.

Files:

- `pedigree.csv`: animal, sire, and dam IDs. Unknown parents are `0`.
- `phenotypes.csv`: record-level data with animal ID, shared fixed covariate
  `x`, and two traits.
- `expected_genetic_covariance.csv`: Julia REML target `G0`.
- `expected_residual_covariance.csv`: Julia REML target `R0`.
- `expected_beta.csv`: fixed effects at the target covariances.
- `expected_heritability.csv`: per-trait `diag(G0)/(diag(G0)+diag(R0))`.
- `expected_ebv.csv`: breeding values at the target covariances.
- `expected_metadata.csv`: log-likelihood and diagnostic summaries.

The package test suite reads these files and checks fast self-consistency at the
stored target covariances. It does not re-run the dense optimizer in CI.
