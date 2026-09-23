# Calcium dynamics and calcium-relevant bone remodeling models

Status: scope-audited literature map (2026-09-22)

Copyright © 2026 azelys/calcium contributors. Canonical project source and
provenance: https://github.com/azelys/calcium. The cited publications and
source artifacts remain the work of their original authors.

## Inclusion rule

The public atlas is deliberately narrower than a general bone-modeling
bibliography. An entry is retained only when it meets one of these tests:

1. **Direct mineral dynamics:** calcium and/or phosphate is an explicit state,
   flux, balance, or clinical output governed by the model.
2. **Direct calcium-regulating submodel:** the model represents PTH-driven
   bone resorption or a PTH/vitamin-D-controlled bone-cell mechanism that is a
   documented component of the Peterson–Riggs calcium-homeostasis lineage.

Models of bone density, estrogen loss, genetics, fatigue, mechanics, or local
cell dynamics alone are not included merely because bone contains calcium.
They need an explicit mineral-dynamics link to belong in this atlas.

## Audited model inventory

| Model | Why it belongs | Calcium-dynamics status | Julia status |
|---|---|---|---|
| Kroll (2000) | Direct PTH timing mechanism: PTH-mediated resorption releases skeletal calcium. | Supporting submodel; no circulating calcium state. | Implemented |
| Lemaire et al. (2004) | PTH/vitamin-D-regulated osteoblast–osteoclast model used in the integrated-model lineage. | Supporting submodel; no circulating calcium state. | Implemented |
| Raposo et al. (2002) | Gut, kidney, parathyroid, calcitriol, plasma Ca/PO4, and exchangeable bone mineral. | Direct | Implemented |
| Peterson & Riggs (2010) | Integrates Ca, PO4, PTH, calcitriol, gut, kidney, bone, and remodeling. | Direct | Implemented variant |
| Peterson & Riggs (2012) BMD extension | Adds a BMD observation layer to the integrated calcium-and-bone framework. | Direct by inheritance from the parent model | Implemented extension |
| Khurana et al. (2019) | Tests PTH(1-84) schedules against normocalcemia and urinary calcium in hypoparathyroidism. | Direct | Implemented variant |
| Granjon, Bonny & Edwards (2017) | Whole-body Ca/PO4 model with PTH, vitamin D3, FGF23, CaSR, kidney, bone, and chemical complexes. | Direct | Implemented |
| Gaweda et al. (2021) | CKD-MBD QSP model calibrated for Ca, phosphorus, PTH, calcitriol, FGF23, and renal progression. | Direct | Literature target |
| Martínez-Reina et al. (2026) | Couples systemic calcium–phosphorus homeostasis to bone remodeling and mineralization in CKD and hyperparathyroidism. | Direct | Literature target |

The first two entries are retained only as clearly labeled supporting
mechanisms. They do **not** claim to simulate serum calcium; their value is
that PTH-regulated resorption and formation mediate exchange between the
skeleton and the regulated extracellular calcium pool in the integrated
models.

## What each model can answer

| Model | Question answered | Sensible next step |
|---|---|---|
| Kroll | Why can intermittent PTH promote formation while continuous exposure promotes resorption? | Couple its cell dynamics to a calcium mass-balance layer. |
| Lemaire | How do PTH, vitamin D, and osteoblast/osteoclast recruitment shape remodeling? | Link cell states to mineral flux and measured turnover markers. |
| Raposo | How do gut, kidney, parathyroid, calcitriol, plasma, and exchangeable bone buffer Ca/PO4 disturbances? | Repair/document source initialization; add FGF23 and CKD. |
| Peterson & Riggs | How do endocrine mineral regulation and remodeling change together in disease or treatment? | Reconcile the as-coded variant with the publication and quantify identifiability. |
| Peterson BMD | How do turnover changes accumulate into a longitudinal BMD signal? | Add multi-site BMD and fracture-risk observation models. |
| Khurana | Which PTH(1-84) schedules balance normocalcemia against hypercalciuria? | Add PK/PD uncertainty, titration, and patient variability. |
| Granjon | How do PTH, vitamin D3, FGF23, CaSR, kidney handling, and mineral complexes regulate Ca/PO4? | Translate from rat to human and add CKD progression. |
| Gaweda | Which mineral/endocrine mechanisms explain CKD-MBD trajectories? | Reproduce the population calibration and recover source parameters. |
| Martínez-Reina | How do dietary Ca/P, PTH, CKD, and bone-cell mineralization interact? | Reproduce the coupled mineralization algorithm and add FGF23. |

## Implementation and validation policy

The Julia library preserves published or public-source behavior first. A
physiological repair or extension must be separate from the source-faithful
mode and must include:

- equation and parameter provenance;
- dimensional, positivity, and mass-balance tests;
- solver tolerances and reproducible scenarios;
- a stated validation target; and
- sensitivity and practical-identifiability assessment before fitting new data.

The Raposo CellML defaults, for example, are reproducible but are not at a
physiological equilibrium. That is recorded as a source-artifact limitation,
not silently repaired by the solver.

## Primary sources

- Kroll (2000), *Parathyroid hormone temporal effects on bone formation and
  resorption*: https://pubmed.ncbi.nlm.nih.gov/10824426/
- Lemaire et al. (2004), *Modeling the interactions between osteoblast and
  osteoclast activities in bone remodeling*:
  https://pubmed.ncbi.nlm.nih.gov/15234198/
- Raposo et al. (2002), *A minimal mathematical model of calcium homeostasis*:
  https://pubmed.ncbi.nlm.nih.gov/12213894/
- Peterson & Riggs (2010), *A physiologically based mathematical model of
  integrated calcium homeostasis and bone remodeling*:
  https://pubmed.ncbi.nlm.nih.gov/19732857/
- Peterson & Riggs (2012), *Predicting nonlinear changes in bone mineral
  density over time using a multiscale systems pharmacology model*:
  https://pubmed.ncbi.nlm.nih.gov/23835796/
- Khurana et al. (2019), *Use of a systems pharmacology model based approach
  toward dose optimization of parathyroid hormone therapy in
  hypoparathyroidism*: https://pubmed.ncbi.nlm.nih.gov/30350311/
- Granjon, Bonny & Edwards (2017), *Coupling between phosphate and calcium
  homeostasis: a mathematical model*:
  https://pubmed.ncbi.nlm.nih.gov/28747359/
- Gaweda et al. (2021), *Development of a quantitative systems pharmacology
  model of chronic kidney disease: metabolic bone disorder*:
  https://pubmed.ncbi.nlm.nih.gov/33308018/
- Martínez-Reina et al. (2026), *The interrelationship between
  calcium-phosphorus homeostasis and bone remodelling*:
  https://doi.org/10.3389/fbioe.2026.1800350

## Source artifacts

- Raposo CellML:
  https://models.cellml.org/workspace/raposo_sobrinho_ferreira_2002/
- Kroll CellML:
  https://models.cellml.org/exposure/e8cf270f437a0af4f69d1fe49304b570/kroll_2000.cellml/view
- Lemaire CellML:
  https://models.physiomeproject.org/exposure/3094a5c3e0810028a9accc3b773cab36/lemaire_tobin_greller_cho_suva_2004.cellml/view
- CaPO4Sim:
  https://github.com/DivadNojnarg/CaPO4Sim
- OpenBoneMin:
  https://github.com/metrumresearchgroup/OpenBoneMin
