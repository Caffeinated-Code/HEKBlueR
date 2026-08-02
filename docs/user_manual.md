# User Manual

## 1. Demo

Open the app and click **Load demo data**. Click **Run analysis**.

For your own experiment, upload:

- raw plate-reader CSV
- plate map CSV
- run metadata CSV

## 2. Upload Data

The raw data should include:

- `plate_id`
- `well`
- `assay_mode`
- `sample_id`
- `peptide_id`
- `control_type`
- `concentration_uM`
- `raw_od`

The plate map should describe what each well contains.

After upload, open **Uploaded Preview** to confirm that HEKBlueR parsed the files correctly.

Open **QC Thresholds** before running analysis if the assay needs thresholds different from the defaults. Each editable threshold shows a recommended range. The selected thresholds are saved with the run and affect the assay identifier.

Open **EDA** before interpreting results. This page checks plate count, well count, missing columns, duplicate wells, missing OD values, saturated OD values, and detected control types.

Open **Cleaned Data** to review blank correction, missing wells, saturated wells, negative corrected values, robust outlier flags, and cleaning actions.

Use the search box in the top navigation bar to highlight matching text on the current page. Every table also has column filters and a table search box. These are useful for `target_id`, `assay_mode`, `peptide_id`, `control_type`, status columns, and reviewer notes.

## 3. Add Metadata

Metadata fields are optional unless marked as required. More metadata gives better reproducibility.

The app reports three metadata scores:

- overall completeness
- required-field completeness
- optional-field completeness

Recommended metadata:

- scientist
- project
- assay date
- target
- cell line
- cell passage
- cell lot
- peptide lot
- peptide purity
- reagent lot
- instrument
- protocol version
- incubation time
- readout wavelength
- notes

## 4. Review Design QC

The app checks whether the experiment has enough controls and replicates.

A warning does not always mean the assay failed. It means the interpretation needs caution.

## 5. Review Plate QC

Open **QC Thresholds** before reviewing automated flags. This page lists the rules used for metadata, design, plate quality, primary hit calls, dose-response review, and counter-assay artifacts.

Use Plate QC to check:

- Z-prime
- robust Z-prime
- SSMD
- signal-to-background
- control CV
- edge effects
- row and column bias
- missing wells
- saturated wells

Open **Reference & Calibration** to check reference control stability and inter-plate calibration. HEKBlueR uses an inter-plate calibrator when available. If not available, it can use shared positive controls for plate alignment review.

## 6. Review Primary Results

Primary Results shows peptide response by dose.

Agonist mode reports activation.

Antagonist mode reports inhibition.

Open the **Sample QC** subtab for one row per target, peptide, and assay mode. This table reports PASS, WARN, and FAIL calls based on primary replicate noise, dose-response QC, and counter-assay artifacts.

## 7. Review Secondary Curves

Secondary Curves shows EC50 or IC50 estimates and curve flags.

Use the subtabs in this order:

1. Interactive curves
2. Curve fit table
3. Curve QC table
4. Replicate noise
5. QC glossary

Curve estimates are weaker when:

- there are too few dose points
- the curve has no plateau
- the estimate is outside the tested range
- replicate noise is high
- viability or interference explains the signal

The detailed curve QC table reports dose count, dynamic range, replicate CV, monotonicity, plateau checks, residual error, and whether the EC50 or IC50 is inside the tested dose range.

## 8. Review Counter-Assays

Counter-Assays shows whether a hit may be caused by:

- cytotoxicity
- peptide color interference
- unrelated reporter activation
- null-cell reporter signal

## 9. Export

Use Final QC to export:

- assay manifest
- final QC table
- normalized results
- full run package

The export package is designed for later database loading.

Each run receives a deterministic assay identifier based on raw data, plate map, metadata, and QC thresholds. If any of those inputs change, a new identifier is assigned.

Major plots include download buttons. File names are generated from project metadata, target metadata, plot type, and date.

## 10. Analysis Cache

The app separates analysis from display. The analysis engine lives in the R modules under `R/`. The Shiny app calls that engine and displays the returned tables and plots.

When **Run analysis** is clicked, the app computes a signature from the raw data, plate map, and metadata. If those inputs have not changed during the current app session, HEKBlueR reuses the cached result object. Searches, tab changes, plot downloads, and interface updates do not rerun the analysis.
