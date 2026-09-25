# Univariate fitted animal-model target fixture (#46)

A Julia-native **fitted** univariate Gaussian animal-model target, serialized by
the HSquared.jl engine fitting its **own** model. No textbook EBVs are typed from
memory; the comparator confrontation is run separately by the R lane
(nadiv/pedigreemm/published) and the JWAS.jl comparator (`comparator/`).

## Model

Single trait `y`, fixed effects `Intercept + x` (one numeric covariate `x`), animal
additive random effect with pedigree relationship `A` (`Ainv = pedigree_inverse`),
REML variance components. 20 animals, multi-generation (Mrode/gryphon-shaped),
interior (non-boundary) optimum.

## Files

- `pedigree.csv` — `animal,sire,dam` (unknown parents = `0`).
- `phenotypes.csv` — `animal,x,y` (one record per animal).
- `expected_variance_components.csv` — `name,value` (`sigma_a2`, `sigma_e2`), the
  engine's REML estimate.
- `expected_beta.csv` — `effect,value` (`Intercept`, `x`) at the REML estimate.
- `expected_ebv.csv` — `id,value` animal breeding values, `id` in pedigree order.
- `expected_reliability.csv` — `id,pev,reliability` (`:selinv`).
- `expected_metadata.csv` — `key,value` (`h2`, `loglik` [full REML, package scale],
  `sigma_a2`, `sigma_e2`, `converged`, `n_animals`, `method`).
- `generate.jl` — the reproducible generator (the engine fits + serializes its own
  output; RNG is used once to realize a dataset with genuine additive structure).

## Regenerate

```sh
julia --project=. test/fixtures/animal_model_fitted_target/generate.jl
```

## Comparator protocol (honest)

The committed CI test (`test/runtests.jl`, "Univariate fitted animal-model target
fixture") rebuilds `Ainv` from `pedigree.csv` and checks **self-consistency**: the
Henderson MME at the stored variance components reproduces the serialized fixed
effects, EBVs, PEV/reliability, and the REML loglik. This is NOT external
validation — it pins the serialized target. Note: `beta` (dense GLS vs sparse MME),
`loglik` (dense REML vs the sparse Henderson identity), and PEV/reliability
(`:selinv` vs `:dense`) each agree across two *distinct* numerical routes; the EBV
check re-solves the *same* Henderson MME, so it is a determinism/integrity pin (it
still catches a corrupted serialized EBV) rather than method-independent
corroboration.

External confrontation is separate and opt-in:

- **R lane:** fit the same model with `nadiv`/`pedigreemm` (REML) and compare to the
  stored targets within a documented tolerance; record the package versions.
- **JWAS.jl:** `comparator/run_jwas_animal_model.jl` (opt-in, `HSQUARED_RUN_JWAS=true`,
  separate env). JWAS is MCMC/Bayesian, so its posterior means agree with these REML
  estimates only **approximately** (different estimators) — report agreement
  honestly, never as "parity"/"validation".

No row is promoted to `covered` until an external comparator run records a tolerance
and an evidence chain. This fixture is a serialized target, not external evidence.
