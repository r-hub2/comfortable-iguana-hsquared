# structured_covariance_parity — diagonal-G multivariate REML target

A deterministic two-trait Julia REML target for the **`:diagonal`** genetic
covariance structure, serialized for the R lane (`hsquared`) to run
sommer / ASReml / BLUPF90 diagonal-G comparators against, and to wire the
**diagonal-vs-unstructured LRT** (issue #42 / #47).

This is the rotation-free subset of the structured-covariance bridge: a diagonal
genetic covariance `G0 = diag(σ²g_1, …, σ²g_t)` has **no loadings and no rotation
ambiguity**, so it is safe to surface across the bridge ahead of the lowrank/fa
rotation/interpretation convention (which stays gated).

## Inputs

`pedigree.csv` and `phenotypes.csv` are the same inputs as
`../phase4_multitrait_parity/` (so the diagonal target can be compared directly
against the unstructured target on identical data). Two correlated-error traits,
intercept + covariate `x`, recorded on pedigreed animals.

## Target (fitted by `fit_multivariate_reml(...; genetic_structure = :diagonal)`)

- `expected_genetic_covariance.csv` — the estimated **diagonal** `G0`
  (off-diagonals are exactly `0`).
- `expected_residual_covariance.csv` — the estimated (unstructured) `R0`.
- `expected_beta.csv` — fixed effects (Intercept, x) per trait.
- `expected_ebv.csv` — breeding values per animal per trait.
- `expected_heritability.csv` — per-trait `h²`.
- `expected_metadata.csv` — `genetic_structure = diagonal`, `n_genetic_params`
  (`= t` for diagonal, vs `t(t+1)/2` unstructured), REML `loglik`, `converged`,
  `iterations`, and the per-trait genetic variances.

## Bridge payload contract

`multivariate_result_payload(fit)` returns the "boring" `NamedTuple` the R twin
marshals: `engine`, `target = "multivariate_reml"`, `genetic_structure`,
`n_traits`, `traits`, `genetic_covariance`, `genetic_variances`,
`residual_covariance`, `genetic_correlation`, `residual_correlation`,
`heritability`, `fixed_effects`, `breeding_values`, `loglik`,
`n_genetic_params`, `converged`. It is exposed only for `:unstructured` and
`:diagonal`; `:lowrank`/`:factor_analytic` are rejected (loadings are
rotation-nonidentified).

## LRT usage

`covariance_structure_lrt(diagonal_fit, unstructured_fit)` on the same data gives
the interior-null (`boundary = false`) diagonal-vs-unstructured test:
`statistic = 2(ℓ_full − ℓ_diag)`, `df = t(t+1)/2 − t = t(t−1)/2`, χ²`df` p-value.

## Status

This is an **input + Julia-target bundle**, not external comparator evidence. The
matching `V4-*` rows move toward `covered` only when an R-lane comparator run
records agreement against these targets (Rose audit each). Experimental,
dense/validation-scale.
