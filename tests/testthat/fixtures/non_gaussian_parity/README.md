# Non-Gaussian Parity Fixture (#44)

Julia-native bridge target for the non-Gaussian `NonGaussianFit` payload
surface. The fixture pins two tiny deterministic fits and their
`nongaussian_result_payload(fit)` outputs so the R lane can write Julia-free
normalizer tests before activating non-Gaussian family syntax.

This is a serialized Julia target, not external comparator evidence and not
public R model-spec activation.

## Cases

- `poisson_laplace`: count response, `family = :poisson`, canonical payload
  `method = "laplace"`, and `n_trials = nothing`.
- `binomial_vector_variational`: success counts with per-record trial totals,
  `family = :binomial`, canonical payload `method = "variational"`, and
  `n_trials` serialized as a semicolon-separated integer vector.

Both cases use the same six-animal pedigree, fixed effects `Intercept + x`,
identity animal incidence, and the low-level dense/validation-scale
`fit_laplace_reml` engine path.

## Files

- `pedigree.csv` - pedigree used to construct `Ainv`.
- `poisson_phenotypes.csv` - IDs, count response, and `x`.
- `binomial_phenotypes.csv` - IDs, successes, per-record trial totals, and `x`.
- `expected_payload_metadata.csv` - payload engine, target, family, method,
  trial denominator encoding, log-likelihood / ELBO value, convergence flag,
  and field lengths.
- `expected_variance_components.csv` - payload variance-component values.
- `expected_fixed_effects.csv` - payload fixed effects.
- `expected_breeding_values.csv` - payload breeding-value IDs and values.
- `generate.jl` - reproducible generator.

## Regenerate

```sh
julia --project=. test/fixtures/non_gaussian_parity/generate.jl
```

## Boundary

The fixture records the bridge payload shape only. It does not activate
non-Gaussian R formula parsing, does not add a public R `family = poisson()` or
`binomial()` path, does not calibrate Bernoulli/binomial intervals, and does not
provide GLLVM/gllvmTMB, ASReml, BLUPF90, MCMCglmm, or other external comparator
parity. The Bernoulli single-trial variance-bias caveat from the validation
ledger still applies outside this fixture.
