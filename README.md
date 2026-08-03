# HEKBlueR

HEKBlueR is a self-service R Shiny app and Nextflow pipeline for HEK-Blue screening analysis.

It helps biologists upload raw plate-reader data, run automated QC, review primary and secondary screening results, detect likely artifacts, and export reproducible tables.

Public app: [https://caffeinated-code.shinyapps.io/HEKBlueR/](https://caffeinated-code.shinyapps.io/HEKBlueR/)

## What It Does

- imports raw HEK-Blue plate-reader data
- lets users choose demo modules: primary screen, secondary confirmation, counter-assay, or any combination
- validates plate maps and experimental design
- previews uploaded raw data, plate maps, and metadata
- filters Uploaded Preview by assay module and physical plate
- runs input EDA with raw data and metadata summaries before interpretation
- keeps analysis logic in reusable R modules separate from the Shiny app
- reuses cached analysis results when the same files and metadata are submitted again in one app session
- creates a cleaned-data review table with cleaning actions
- separates table-heavy pages into tabs for easier review
- scores required and optional metadata completeness separately
- checks controls, replicates, dose series, and metadata completeness
- exports the active QC threshold table used for every automated flag
- requires a rationale when default QC thresholds are changed
- calculates plate QC metrics
- checks reference control stability
- performs inter-plate calibration when a shared calibrator or shared positive control is present
- normalizes agonist and antagonist responses
- fits simple dose-response curves
- reports detailed curve QC, including dose count, curve range, residuals, plateaus, Hill slope, and replicate noise
- reports sample-level PASS, WARN, and FAIL calls based on explicit QC thresholds
- assigns a readable deterministic assay identifier from target, assay type, peptide IDs, assay date, scientist, and a six-digit change code
- provides an interactive dose-response workspace with curve plots, fit tables, QC tables, replicate noise, and QC glossary
- flags noisy or weak curves
- flags likely assay artifacts
- creates major review plots
- exports major plots as named PNG files
- exports database-ready result tables, figures, run documentation, and a QC report
- runs locally through Shiny or in batch through Nextflow

## Demo Scope

The demo dataset contains:

- 1 target
- 6 peptides
- 2 primary screening plates with fixed-dose treatment wells
- 2 secondary dose-response plates
- 1 counter-assay plate
- viability and assay-interference examples
- one artifact-prone peptide

## Run The Shiny App

```r
shiny::runApp()
```

Then click **Load demo data** and **Run analysis**.

## Run The CLI Pipeline

```bash
Rscript scripts/run_pipeline_cli.R \
  --raw data/simulated/raw_plate_reader.csv \
  --plate-map data/simulated/plate_map.csv \
  --metadata data/simulated/run_metadata.csv \
  --out results/demo_cli
```

## Run With Nextflow

```bash
nextflow run main.nf -profile local \
  --samplesheet data/simulated/samplesheet.csv \
  --outdir results/nextflow_demo
```

The samplesheet is a CSV with one row per run:

| Column | Meaning |
|---|---|
| `run_id` | name used for the output folder |
| `raw` | path to the raw plate-reader CSV |
| `plate_map` | path to the plate map CSV |
| `metadata` | path to the run metadata CSV |

You can also run one dataset without a samplesheet:

```bash
nextflow run main.nf -profile local \
  --raw data/simulated/raw_plate_reader.csv \
  --plate_map data/simulated/plate_map.csv \
  --metadata data/simulated/run_metadata.csv \
  --outdir results/nextflow_demo
```

AWS Batch profiles are included in `nextflow.config`.

## App Tabs

| Tab | Purpose |
|---|---|
| Demo | load demo data and run analysis |
| Upload | upload raw data, plate map, and metadata |
| QC Thresholds | adjust review thresholds before running analysis |
| Uploaded Preview | inspect uploaded raw data, plate map, liquid-handler map, and metadata by assay module and plate |
| EDA | review input-level exploratory checks |
| Cleaned Data | inspect cleaning actions, outliers, and cleaned wells |
| Metadata | review run documentation completeness |
| Design QC | check controls, replicates, dose points, and calibrators |
| Plate QC | review Z-prime, CV, SSMD, edge effects, and plate bias |
| Reference & Calibration | review reference control stability and inter-plate calibration |
| Plate Layout | inspect one plate at a time with raw and normalized heatmaps |
| Primary Results | review activation, inhibition, and hit calls |
| Secondary Curves | review EC50 or IC50 and curve QC |
| Counter-Assays | review cytotoxicity and assay-interference flags |
| Plots | create standard and custom plots |
| Final QC | export the final decision table and result package |

## Output Tables

Each run can export individual tables and a structured ZIP package.

The ZIP package contains:

```text
ASSAY_IDENTIFIER_results/
  qc_report.md
  documentation/
  tables/
  figures/
```

Core tables include:

- `metadata.csv`
- `assay_manifest.csv`
- `run_documentation.csv`
- `qc_thresholds.csv`
- `raw_data_summary.csv`
- `metadata_summary.csv`
- `raw_data.csv`
- `plate_map.csv`
- `cleaned_well_data.csv`
- `normalized_results.csv`
- `design_qc.csv`
- `plate_qc.csv`
- `primary_results.csv`
- `sample_qc_table.csv`
- `dose_response_results.csv`
- `dose_response_qc.csv`
- `counter_assay_qc.csv`
- `hit_calls.csv`
- `exclusions.csv`
- `final_qc_table.csv`
- `run_summary.json`
- `analysis_config.yml`

## Why Simulated Data Is Used

Public HEK-Blue assay results are often available as summarized bioactivity values. Raw plate-reader files with full plate maps, control wells, reagent metadata, and counter-assay layouts are less commonly released.

The demo uses simulated raw plate data so the complete QC workflow can be shown without missing design information. Many real primary HTS runs use one fixed concentration. Dose-response curves are usually run later for confirmation. HEKBlueR supports partial uploads, including primary-only, secondary-only, counter-assay-only, or combined runs.

## Documentation

- [Assay primer](docs/assay_primer.md)
- [User manual](docs/user_manual.md)
- [Analysis walkthrough](docs/analysis_walkthrough.md)
- [QC metrics](docs/qc_metrics.md)
- [QC threshold rationale](docs/qc_threshold_rationale.md)
- [Data dictionary](docs/data_dictionary.md)
- [AWS and Nextflow scaling](docs/aws_nextflow_scaling.md)
- [Public data sources](docs/public_data_sources.md)
- [Implementation plan](docs/implementation_plan.md)

## Analysis And App Separation

The analysis engine lives in `R/analysis.R`, `R/qc_metrics.R`, `R/qc_thresholds.R`, `R/plots.R`, and `R/export.R`.

The Shiny app in `app.R` is a front end over those modules. It reads the input files, computes a run signature from the raw data, plate map, and metadata, and reuses the cached result object when the inputs have not changed. UI-only changes do not require rethinking the analysis code.

The same analysis engine is used by:

- the local Shiny app
- `scripts/run_pipeline_cli.R`
- the Nextflow workflow in `main.nf`
- the AWS Batch-ready Nextflow profile in `nextflow.config`

## Key Sources

- [InvivoGen QUANTI-Blue](https://www.invivogen.com/quanti-blue)
- [InvivoGen HEK-Blue Detection](https://www.invivogen.com/hek-blue-detection)
- [HEK-Blue TLR8 antagonist protocol](https://pmc.ncbi.nlm.nih.gov/articles/PMC8715332/)
- [Z-prime factor paper](https://pubmed.ncbi.nlm.nih.gov/10838414/)
- [Nextflow AWS documentation](https://docs.seqera.io/nextflow/aws)
