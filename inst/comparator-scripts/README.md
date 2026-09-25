# Comparator Scripts

These scripts are manual gates for external multivariate animal-model
comparators. They are not ordinary tests and they do not create validation
claims by themselves.

Current scope:

- `asreml/multivariate-animal.R`: prepares the shared Phase 4 fixture for a
  candidate ASReml-R multivariate animal model. Run with `--run` only on a
  machine with a licensed ASReml-R installation.
- `blupf90/prepare-multivariate-animal.R`: prepares BLUPF90-family flat files
  and parameter templates from the shared Phase 4 fixture.
- `blupf90/multivariate-animal.renf90`: RENUMF90-oriented template.
- `blupf90/multivariate-animal.par`: AIREMLF90/BLUPF90+ application template.

Public-claim rule:

```text
No ASReml/BLUPF90 parity claim exists until a run output is recorded under
docs/dev-log/comparator-runs/ with tool version, input checksum, convergence
status, covariance estimates, and Rose/Fisher/Curie review.
```

These scripts intentionally use the shared tiny Phase 4 fixture first. Larger
recovery and production-scale comparator studies are separate work.
