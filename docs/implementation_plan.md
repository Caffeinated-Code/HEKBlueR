# HEKBlueR Implementation Plan

HEKBlueR is a self-service HEK-Blue screening analysis system. It is designed to remove manual plate review while preserving scientific traceability.

## Main Goal

Upload raw HEK-Blue assay files, plate maps, and metadata. Produce a complete QC-reviewed result package with clean tables, plots, hit calls, curve fits, and reproducibility files.

## Architecture

The project has three layers.

| Layer | Purpose |
|---|---|
| Shiny app | Biologist-facing upload, review, visualization, and export |
| R analysis engine | Reusable analysis functions for QC, normalization, hit calling, and dose-response |
| Nextflow pipeline | Scalable batch execution locally or on AWS Batch |

## App Workflow

1. Upload raw OD data.
2. Upload plate map.
3. Add target, peptide, assay, reagent, and run metadata.
4. Validate experimental design.
5. Clean raw well data.
6. Calculate plate and control QC.
7. Normalize agonist, antagonist, and counter-assay responses.
8. Fit dose-response curves.
9. Flag artifacts and design issues.
10. Export database-ready results.

## Assay Modes

### Agonist

Peptide alone activates the HEK-Blue reporter. The analysis compares peptide wells to vehicle or baseline controls and a positive agonist control.

### Antagonist

Peptide is tested with a fixed agonist challenge. The analysis compares peptide plus agonist wells to agonist-only and baseline controls.

### Counter-Assay

The app identifies likely artifacts from cytotoxicity, no-cell peptide interference, unrelated reporter activity, null-cell activity, or reagent-only signal.

## Experimental Design QC

The first QC layer checks whether the experiment can support the intended interpretation.

| Check | Pass | Warning | Fail |
|---|---:|---:|---:|
| Technical replicates | >= 3 | 2 | 1 |
| Biological replicates | >= 3 | 2 | missing |
| Dose points | >= 8 | 5 to 7 | < 5 |
| Control wells per class | >= 8 | 4 to 7 | < 4 |
| Metadata completeness | >= 80% | 50 to 79% | < 50% |
| Inter-plate calibrator | every plate | partial | absent |

## Plate QC

Core metrics:

- Z-prime
- robust Z-prime
- SSMD
- signal-to-background
- signal window
- positive control CV
- negative control CV
- blank CV
- replicate CV
- edge effect score
- row bias
- column bias
- inter-plate calibrator drift
- saturated OD count
- missing well count

## Dose-Response QC

Curve flags:

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

## AWS Scaling Path

The Shiny app can run locally first. For scale:

- store raw uploads and results in S3
- store searchable summaries in RDS Postgres
- run Nextflow on AWS Batch
- use Spot instances for compute tasks
- use a small on-demand head job
- use CloudWatch for logs
- use Secrets Manager for credentials
- use Cognito or SSO for login

## Output Package

Each run exports:

```text
metadata.json
analysis_config.yml
raw_data.csv
plate_map.csv
cleaned_well_data.csv
normalized_results.csv
design_qc.csv
plate_qc.csv
control_qc.csv
replicate_qc.csv
dose_response_results.csv
dose_response_qc.csv
counter_assay_qc.csv
hit_calls.csv
exclusions.csv
final_qc_table.csv
report.html
plots/
```

