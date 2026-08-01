# HEKBlueR

HEKBlueR is a self-service R Shiny app and Nextflow pipeline for HEK-Blue screening analysis.

It helps biologists upload raw plate-reader data, run automated QC, review primary and secondary screening results, detect likely artifacts, and export reproducible tables.

## What It Does

- imports raw HEK-Blue plate-reader data
- validates plate maps and experimental design
- checks controls, replicates, dose series, and metadata completeness
- calculates plate QC metrics
- normalizes agonist and antagonist responses
- fits simple dose-response curves
- flags noisy or weak curves
- flags likely assay artifacts
- creates major review plots
- exports database-ready result tables
- runs locally through Shiny or in batch through Nextflow

## Demo Scope

The demo dataset contains:

- 1 target
- 3 peptides
- primary agonist plate
- antagonist dose-response plate
- counter-assay plate
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
nextflow run main.nf -profile local --outdir results/nextflow_demo
```

AWS Batch profiles are included in `nextflow.config`.

## App Tabs

| Tab | Purpose |
|---|---|
| Start | load demo data and run analysis |
| Upload | upload raw data, plate map, and metadata |
| Metadata | review run documentation completeness |
| Design QC | check controls, replicates, dose points, and calibrators |
| Plate QC | review Z-prime, CV, SSMD, edge effects, and plate bias |
| Plate Layout | inspect raw and normalized heatmaps |
| Primary Results | review activation, inhibition, and hit calls |
| Secondary Curves | review EC50 or IC50 and curve QC |
| Counter-Assays | review cytotoxicity and assay-interference flags |
| Plots | create standard and custom plots |
| Final QC | export the final decision table and result package |

## Output Tables

Each run can export:

- `metadata.csv`
- `raw_data.csv`
- `plate_map.csv`
- `cleaned_well_data.csv`
- `normalized_results.csv`
- `design_qc.csv`
- `plate_qc.csv`
- `primary_results.csv`
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

The demo uses simulated raw plate data so the complete QC workflow can be shown without missing design information.

## Documentation

- [Assay primer](docs/assay_primer.md)
- [User manual](docs/user_manual.md)
- [QC metrics](docs/qc_metrics.md)
- [Data dictionary](docs/data_dictionary.md)
- [AWS and Nextflow scaling](docs/aws_nextflow_scaling.md)
- [Public data sources](docs/public_data_sources.md)
- [Implementation plan](docs/implementation_plan.md)

## Key Sources

- [InvivoGen QUANTI-Blue](https://www.invivogen.com/quanti-blue)
- [InvivoGen HEK-Blue Detection](https://www.invivogen.com/hek-blue-detection)
- [HEK-Blue TLR8 antagonist protocol](https://pmc.ncbi.nlm.nih.gov/articles/PMC8715332/)
- [Z-prime factor paper](https://pubmed.ncbi.nlm.nih.gov/10838414/)
- [Nextflow AWS documentation](https://docs.seqera.io/nextflow/aws)

