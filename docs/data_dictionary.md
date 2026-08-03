# Data Dictionary

## Raw Data

| Column | Meaning |
|---|---|
| `plate_id` | unique plate identifier |
| `well` | well coordinate, such as A01 |
| `row` | plate row |
| `col` | plate column |
| `assay_stage` | primary, secondary, or counter |
| `assay_mode` | agonist, antagonist, counter, or unknown direction |
| `sample_id` | sample or control name |
| `peptide_id` | peptide identifier |
| `control_type` | well role |
| `concentration_uM` | peptide concentration |
| `technical_replicate` | technical replicate number |
| `biological_replicate` | biological replicate number |
| `raw_od` | raw optical density |
| `expected_activity` | expected compound role, such as agonist, antagonist, or unknown |

## Control Types

Supported values:

- `blank`
- `negative_control`
- `positive_control`
- `agonist_challenge_control`
- `known_antagonist_control`
- `inter_plate_calibrator`
- `test_sample`
- `viability_counter`
- `no_cell_interference`
- `unrelated_reporter`
- `null_cell_reporter`
- `empty`

## Metadata

Metadata uses three columns:

- `field`
- `value`
- `required`

This format is easy to parse into JSON, CSV, SQLite, DuckDB, or Postgres.

## Assay Manifest

| Column | Meaning |
|---|---|
| `assay_identifier` | deterministic run ID based on raw data, plate map, metadata, and QC thresholds |
| `change_code` | six-digit code that changes when inputs, metadata, thresholds, rationale, or analysis strategy change |
| `input_signature` | combined hash for all inputs and thresholds |
| `raw_data_signature` | hash for raw well-level data |
| `plate_map_signature` | hash for plate map or well annotations |
| `metadata_signature` | hash for submitted metadata |
| `threshold_signature` | hash for selected QC thresholds |
| `threshold_changed_from_default` | TRUE when active QC thresholds differ from the defaults |
| `threshold_change_note` | required rationale when default QC thresholds were changed |
| `assay_type` | assay modules found in the raw data |
| `peptide_ids` | peptide or compound IDs found in the raw data |

## Run Documentation

`run_documentation.csv` combines:

- assay identifier
- computed input signatures
- submitted metadata
- selected QC thresholds
- threshold-change rationale
- analysis parameters

This file is intended for audit trails and later database loading.
