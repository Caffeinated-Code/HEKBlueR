# QC Threshold Rationale

HEKBlueR uses conservative defaults for cell-based screening, then lets scientists adjust thresholds when assay validation data supports a different rule.

## Defaults

| Area | Default | Reason |
|---|---:|---|
| Z-prime pass | 0.5 | Common HTS benchmark for strong separation between positive and negative controls. |
| Z-prime warning floor | 0 | Values from 0 to 0.5 are marginal, but may still be reviewable in complex cell-based assays. |
| Control CV warning | 20% | Common assay guidance target for control stability. |
| Intra-plate CV warning | 20% | Extends control CV review to within-plate variability. |
| Edge or spatial bias warning | 0.15 | Flags edge, row, or column effects large enough to distort well-level interpretation. |
| Outlier rate warning | 5% | Flags plates with broad well-level review burden. |
| Primary replicate CV warning | 25% | Keeps primary hit calls stricter than exploratory visualization. |
| Primary replicate CV fail | 50% | Prevents strong conclusions from highly unstable dose-level replicates. |
| Dose points pass | 8 | Practical default for secondary curve fitting. |
| Curve replicate CV warning | 30% | Dose-response curves tolerate some noise, but high replicate noise weakens potency estimates. |
| Hill slope range | 0.2 to 4 | Very shallow or steep slopes often indicate weak biology, artifacts, or unstable fits. |

## Intra-Plate Variability

Intra-plate variability should be considered. Z-prime summarizes separation between control groups, but it does not fully describe where variation occurs on a plate.

HEKBlueR therefore adds `intraplate_variability_qc.csv`, which checks:

- within-plate control CV
- edge effect
- row bias
- column bias
- flagged-well rate

These metrics help identify evaporation, dispense order effects, incubation gradients, reader artifacts, and local outlier burden.

## Sources

- Zhang et al. introduced Z-prime as a screening assay quality metric: https://pubmed.ncbi.nlm.nih.gov/10838414/
- The Assay Guidance Manual is the core reference for assay development and HTS practice: https://www.ncbi.nlm.nih.gov/books/NBK53196/
- NCBI assay guidance notes that Z-prime above 0.5 is common, while complex phenotype assays may still be reviewed in the 0 to 0.5 range: https://www.ncbi.nlm.nih.gov/books/NBK126174/
- Assay operations guidance emphasizes monitoring EC50/IC50 variability and dose-response curve quality: https://www.ncbi.nlm.nih.gov/books/NBK91994/
- CDD plate QC guidance summarizes practical Z-prime interpretation for screening workflows: https://support.collaborativedrug.com/hc/en-us/articles/214359383-Plate-Quality-Control
