# Marker Scan Parity Fixture (#45)

R mirror of the Julia-native post-fit marker-scan payload fixture from
`HSquared.jl` PR #142 (`f9fbbb1`). The fixture pins a small
relatedness-corrected `mixed_model_marker_scan(fit, markers)` result and the
`marker_scan_result_payload(scan)` row shape so R can verify payload
normalization without live Julia execution.

This is a serialized Julia target, not calibrated GWAS validation and not
external comparator evidence.

## Files

- `phenotypes.csv` - animal IDs and response values.
- `pedigree.csv` - pedigree used by the Julia generator to construct `Ainv`.
- `markers.csv` - marker dosages, rows aligned to `phenotypes.csv`.
- `expected_marker_scan_payload.csv` - row-aligned payload fields: marker IDs,
  effects, SEs, z-scores, chi-square, p-values, Bonferroni/BH values, LOD,
  denominator, and allele frequency.
- `expected_metadata.csv` - payload target, variance components, marker count,
  and VanRaden scale.

## Boundary

The fixture records nominal Wald scan output plus deterministic Bonferroni/BH
adjustments over the supplied marker set. It does not calibrate genome-wide
thresholds, run LOCO as a public workflow, activate formula-level
`marker_scan()`, join marker maps, draw plots, or promote GWAS/QTL/eQTL support
to covered.
