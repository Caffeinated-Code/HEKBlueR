# Data Dictionary

## Raw Data

| Column | Meaning |
|---|---|
| `plate_id` | unique plate identifier |
| `well` | well coordinate, such as A01 |
| `row` | plate row |
| `col` | plate column |
| `assay_mode` | agonist, antagonist, or counter |
| `sample_id` | sample or control name |
| `peptide_id` | peptide identifier |
| `control_type` | well role |
| `concentration_uM` | peptide concentration |
| `technical_replicate` | technical replicate number |
| `biological_replicate` | biological replicate number |
| `raw_od` | raw optical density |

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

