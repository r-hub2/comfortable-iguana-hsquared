# Hinv-cell R↔engine parity fixture (design-41 §3 #8)

**Status: receipt fixture. NOT a covered flip.**  
`V2-SSHINV` stays **partial**. `public_covered_count` stays **7**.  
Experimental version stays **0.7.0**. Ordinary-route `single_step()` stays held.

```
PLATFORM: cursor | LANE: cursor/08-ss-hinv-parity-20260903
OTHER LANES: R #164 honesty cite-only · Julia #295 cite-only · Codex DRAFT cite-only
```

## What this closes

Rose §3 #8 after Darwin SIGN + `nod SS Boole`: a **committed** element-wise
R↔engine receipt on the **covered cell**, not VanRaden `markers=`.

Live `test-single-step-construct.R` (reorder / all-id / differs-from-pedigree)
does **not** substitute. Those tests build VanRaden `G` from markers
(design-56 §A.2 clause 8 forbids that as the covered-claim `G`).

## Covered cell (design-56 §A.2)

| Knob | Frozen value |
|---|---|
| Path | supplied-`Hinv` (`target = "single_step"`) |
| `G` | `A₂₂ + 0.05 I` (teaching kernel; **not** VanRaden) |
| `tau`, `omega` | 1, 1 |
| `blend_weight`, `ridge` | 0, 0 |
| Pedigree | n=6 half-sib `s1,d1,d2,o1,o2,o3`; genotyped `o1,o2,o3` (rows 4,5,6) |

## Four objects

| Object | How it is pinned | Tolerance |
|---|---|---|
| `Hinv` | this directory `engine_hinv.csv` vs live `single_step_inverse` | `max\|Δ\| ≤ 1e-10` (packet AGREE was `4.24e-12`) |
| `sigma_a2` | live R `single_step(Hinv=)` vs live Julia `fit_single_step_reml` on the **same** `y` | `1e-8` relative |
| `sigma_e2` | same | `1e-8` relative |
| labelled GEBVs | same; ids = pedigree ids | `1e-8` relative, ids exact |

`y` is **not** the n=6 recovery-smoke draw (`engine_recovery_smoke.csv`).
That dump has no `y` / GEBVs, Julia `randn` is not version-stable, and the
smoke is an **unidentified** recovery dump (σ²a ≈ 0.076 vs truth 1.0). This
receipt is **R↔engine agreement**, not recovery-to-truth. Do not quote the
smoke VCs as this receipt.

## Fixture SHA

Julia AGHmatrix construction packet on HSquared.jl [#295](https://github.com/itchyshin/HSquared.jl/pull/295):

- Packet commit: **`0b03d67e`** (`comparator/aghmatrix_hmatrix/`)
- Live PR tip cited by the SS runner: **`e9676014`** (not merged)
- Construction AGREE: `max|Hinv Δ| = 4.24e-12` vs AGHmatrix 3.0.1 Martini τ=ω=1
- Mrode Ch.11 **NO-ANCHOR** (in `metadata.csv`)

Copied files: `A.csv`, `G.csv`, `engine_hinv.csv`, `metadata.csv`.
`ENGINE_PACKET_SHA.txt` repeats `0b03d67e`.

## Comparator evidence (what this is / is not)

- **Is:** R supplied-`Hinv` bridge (`fit_ai_reml` on that `Hinv`) vs engine
  `fit_single_step_reml` that rebuilds the same teaching-kernel `Hinv`.
- **Is:** Hinv construction pin against the #295 AGHmatrix packet.
- **Is not:** AGHmatrix construction AGREE re-badged as REML / fit parity.
- **Is not:** n=240 recovery GATE PASS (`8e6e038b` / `0533e9da`).
- **Is not:** preGSf90 / blupf90+ numbers.
- **Is not:** R-public covered, ordinary-route activation, count 8, or a version bump.

## Test

`tests/testthat/test-single-step-hinv-parity.R`

- Julia-free: fixture shape + SHA + teaching-kernel metadata.
- Live (skip-guarded): the four-object comparison.

## Explicit non-claims

No Darwin SIGN. No `nod SS Boole`. No Rose CLEAN. No field-4 flip.
