# Public Data Sources For HEK-Blue And Screening Assays

Raw HEK-Blue plate-reader data with full plate maps is not commonly released. Public resources usually provide summarized assay activity, dose-response values, or curated bioactivity records.

HEKBlueR uses simulated data for the demo because full QC requires raw wells, plate maps, controls, reagent metadata, and counter-assay labels.

## PubChem BioAssay

PubChem BioAssay contains public screening assay data and supports downloads through web and FTP routes.

Useful for:

- assay descriptions
- active and inactive calls
- dose-response summaries
- confirmatory assay records

Limitations:

- raw plate maps are often absent
- reagent lots and run metadata are usually absent
- peptide-focused HEK-Blue screens may be difficult to find

Search terms:

- `HEK-Blue`
- `SEAP`
- `QUANTI-Blue`
- `TLR antagonist`
- `NF-kB reporter`
- `AP-1 reporter`

Resource:

- https://pubchem.ncbi.nlm.nih.gov/docs/bioassays

## ChEMBL

ChEMBL is a curated bioactivity database. It is useful for target, compound, assay, and activity summaries.

Useful for:

- curated assay descriptions
- EC50, IC50, potency, and activity records
- target-linked bioactivity
- literature-derived screening results

Limitations:

- usually not raw plate-level data
- control wells and plate layouts are not usually available

Resource:

- https://www.ebi.ac.uk/chembl/

## BindingDB

BindingDB focuses on measured binding and inhibition data.

Useful for:

- target-ligand activity records
- inhibition values
- binding constants
- medicinal chemistry context

Limitations:

- not usually HEK-Blue plate-reader data
- better for binding and potency context than plate QC demos

Resource:

- https://www.bindingdb.org/

## Literature Supplements

Some papers include HEK-Blue assay spreadsheets as supplementary data.

Search terms:

- `HEK-Blue supplementary raw data`
- `QUANTI-Blue plate reader supplementary`
- `HEK-Blue TLR8 antagonist supplementary`
- `SEAP reporter assay raw OD`

Limitations:

- formats vary
- plate maps may be missing
- control definitions may be embedded in methods text

## Best Public Demo Strategy

Use simulated raw HEK-Blue data to teach the full pipeline.

Use PubChem, ChEMBL, BindingDB, and literature supplements to show where real assay summaries can be sourced.

When raw public plate data is found, convert it into the HEKBlueR schema:

- `raw_plate_reader.csv`
- `plate_map.csv`
- `run_metadata.csv`

