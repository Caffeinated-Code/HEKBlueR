# User Manual

## 1. Start

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

## 3. Add Metadata

Metadata fields are optional unless marked as required. More metadata gives better reproducibility.

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

## 6. Review Primary Results

Primary Results shows peptide response by dose.

Agonist mode reports activation.

Antagonist mode reports inhibition.

## 7. Review Secondary Curves

Secondary Curves shows EC50 or IC50 estimates and curve flags.

Curve estimates are weaker when:

- there are too few dose points
- the curve has no plateau
- the estimate is outside the tested range
- replicate noise is high
- viability or interference explains the signal

## 8. Review Counter-Assays

Counter-Assays shows whether a hit may be caused by:

- cytotoxicity
- peptide color interference
- unrelated reporter activation
- null-cell reporter signal

## 9. Export

Use Final QC to export:

- final QC table
- normalized results
- full run package

The export package is designed for later database loading.

