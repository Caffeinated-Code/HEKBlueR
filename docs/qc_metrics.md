# QC Metrics

HEKBlueR has four QC layers.

## Design QC

Design QC checks whether the experiment can support interpretation.

| Metric | Pass | Warning | Fail |
|---|---:|---:|---:|
| Technical replicates | >= 3 | 2 | 1 |
| Biological replicates | >= 3 | 2 | missing |
| Dose points | >= 8 | 5 to 7 | < 5 |
| Control wells per class | >= 8 | 4 to 7 | < 4 |
| Metadata completeness | >= 80% | 50 to 79% | < 50% |
| Inter-plate calibrator | every plate | partial | absent |

## Plate QC

| Metric | Purpose |
|---|---|
| Z-prime | control separation and assay quality |
| robust Z-prime | median-based control separation |
| SSMD | standardized separation between controls |
| signal-to-background | positive signal over baseline |
| signal window | positive minus negative control |
| control CV | control stability |
| replicate CV | technical replicate agreement |
| edge effect score | edge vs center difference |
| row bias | row-level drift |
| column bias | column-level drift |
| inter-plate calibrator drift | plate-to-plate comparability |

Starting thresholds:

- Z-prime at least 0.5 is preferred.
- Z-prime from 0.3 to 0.5 needs review.
- Z-prime below 0.3 is usually a failed or weak plate.
- Control CV above 20 percent needs review.
- Replicate CV above 20 percent needs review.

## Dose-Response QC

Flags:

- `GOOD_CURVE`
- `WEAK_RESPONSE`
- `NO_TOP_PLATEAU`
- `NO_BOTTOM_PLATEAU`
- `EC50_OUT_OF_RANGE`
- `HIGH_REPLICATE_NOISE`
- `NON_MONOTONIC`
- `BAD_HILL_SLOPE`
- `CYTOTOXICITY_CONFOUNDED`
- `ASSAY_INTERFERENCE`
- `LOW_CONFIDENCE_HIT`

## Final Actions

| Action | Meaning |
|---|---|
| `PASS_ADVANCE` | result is strong enough to advance |
| `PASS_RETEST` | signal is promising but should be repeated |
| `REVIEW_CONTROL_DRIFT` | controls changed across plate or run |
| `REVIEW_CURVE_NOISE` | curve fit is weak or noisy |
| `FAIL_PLATE_QC` | plate metrics do not support interpretation |
| `FAIL_DESIGN_QC` | missing design elements block interpretation |
| `LIKELY_ARTIFACT` | counter-assay suggests false activity |
| `INSUFFICIENT_REPLICATES` | replicate structure is too weak |

