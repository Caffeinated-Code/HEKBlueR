# Analysis Walkthrough

This walkthrough explains what each app step does and why it matters.

## 1. Demo Or Upload

Use the demo checkboxes to choose a primary screen, secondary confirmation screen, counter-assay, or any combination.

For lab data, upload raw well-level plate-reader data. A plate map is recommended. If the raw file already contains well annotations, the app can use those annotations.

## 2. QC Thresholds

QC thresholds define how the app assigns PASS, WARN, and FAIL calls.

Defaults are suitable starting points for HEK-Blue-style cell assays. Scientists can adjust them when assay validation data supports a different threshold.

The selected thresholds become part of the assay identifier.

## 3. Uploaded Preview

This page answers one question: did the app read the files correctly?

Use the assay module and plate filters first. Then check that plate IDs, wells, assay stages, assay modes, sample IDs, control labels, concentrations, and metadata look right.

## 4. EDA

EDA summarizes the upload before biological interpretation.

Raw data summary checks:

- number of rows
- number of plates
- assay modules detected
- number of tested compounds
- control type diversity
- dose points
- replicate labels
- raw OD range

Metadata summary checks:

- required fields
- optional fields
- blank metadata values
- key fields such as project, scientist, assay date, target, and cell line

## 5. Cleaned Data

Cleaning flags wells that need review.

Common flags include missing raw OD, saturated OD, negative values after blank correction, and robust outliers.

## 6. Design QC

Design QC asks whether the experiment can support interpretation.

It checks controls, replicate counts, dose points, inter-plate calibrators, and metadata completeness.

## 7. Plate QC

Plate QC checks whether the plate behaved like a reliable assay system.

Important metrics include Z-prime, robust Z-prime, SSMD, control CV, edge effect, row bias, column bias, and flagged-well rate.

Intra-plate variability is included because one good plate-level score does not rule out local spatial artifacts.

## 8. Plate Layout

Plate Layout shows each well in its physical plate position.

Hover over a well to see:

- plate ID
- well
- sample
- compound or peptide
- control type
- assay mode
- concentration
- signal value

Use this page to find edge effects, dispense patterns, missing wells, or dose layout issues.

## 9. Primary Results

Primary Results summarizes activation or inhibition at each dose.

The waterfall plot helps identify high-response compounds. The Sample QC table converts the evidence into PASS, WARN, and FAIL review rows.

## 10. Replicate Noise

Replicate noise is measured with technical replicate CV.

Low CV means replicate wells agree. High CV means the average response at that dose may be unreliable.

High CV near the active dose range can weaken hit calls, EC50 estimates, and IC50 estimates.

## 11. Secondary Curves

Secondary Curves fit dose-response models when enough dose points are available.

Review EC50 or IC50 together with:

- dose count
- dynamic range
- replicate CV
- residual error
- monotonicity
- plateau coverage
- Hill slope

## 12. Counter-Assays

Counter-assays help separate true biology from artifacts.

HEKBlueR checks cytotoxicity, no-cell interference, unrelated reporter activity, and null-cell activity.

## 13. Final QC And Exports

Final QC is split into Assay manifest, Run Documentation, Final QC table, and Results Download.

The export package includes:

- assay manifest
- run documentation
- metadata
- raw data
- plate map
- QC thresholds
- cleaned data
- normalized results
- plate QC
- intra-plate variability QC
- primary results
- secondary curve results
- counter-assay QC
- final QC table
- QC report
- generated figures

These files are intended to be searchable and database-ready.

The ZIP download uses this structure:

```text
ASSAY_IDENTIFIER_results/
  qc_report.md
  documentation/
  tables/
  figures/
```
