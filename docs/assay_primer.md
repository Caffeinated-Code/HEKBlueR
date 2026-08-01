# HEK-Blue Assay Primer

HEK-Blue assays use engineered HEK cells with a reporter system. When the target pathway is activated, the cells secrete SEAP. SEAP is detected with a colorimetric reagent such as QUANTI-Blue.

The readout is optical density. QUANTI-Blue and HEK-Blue Detection workflows commonly use OD at 620 to 655 nm.

## Agonist Assay

An agonist assay asks whether a peptide activates the reporter.

Expected controls:

- blank
- cells plus vehicle
- positive agonist
- test peptide
- optional no-cell peptide control
- optional unrelated reporter line

Output:

- percent activation
- fold over baseline
- robust hit flag

## Antagonist Assay

An antagonist assay asks whether a peptide blocks a fixed agonist challenge.

Expected controls:

- blank
- cells plus vehicle
- agonist-only challenge
- known antagonist if available
- peptide plus agonist
- peptide alone control
- viability counter-assay

Output:

- percent inhibition
- IC50
- Imax
- curve quality

## Counter-Assays

Counter-assays protect against false hits.

Useful counter-assays:

- viability assay
- no-cell peptide plus reagent wells
- unrelated reporter cell line
- parental or null reporter cell line
- reagent-only wells
- vehicle-only controls
- endotoxin-sensitive control when innate immune activation is plausible

## Common Artifacts

- peptide is cytotoxic
- peptide directly changes reagent color
- peptide activates a broad stress pathway
- peptide activates unrelated reporter cells
- edge wells evaporate
- controls drift during long incubations
- positive control loses activity
- plate reader saturation hides the true signal

