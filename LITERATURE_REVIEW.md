# Calcium homeostasis and bone-remodeling models

Status: initial literature map and replication specification (2026-09-22)

Copyright © 2026 azelys/calcium contributors. Canonical project source and
provenance: https://github.com/azelys/calcium. The publications and model
artifacts reviewed below remain attributed to their original authors.

## Scope

This review focuses on mathematical or computational models that make explicit predictions about calcium/phosphate homeostasis, bone-cell dynamics, bone remodeling, or interventions affecting those systems. A paper is a replication candidate when its equations, parameterization, initial conditions, and enough simulation detail are recoverable from the paper, supplement, or an associated model repository.

The field is not one homogeneous model family. It spans organ-level endocrine feedback, well-mixed cell-population ODEs, spatial BMU models, and integrated quantitative-systems-pharmacology (QSP) models. We should therefore replicate canonical model families rather than claim that one implementation represents every publication.

## Working taxonomy

| Layer | Representative models | Main question | Julia representation |
|---|---|---|---|
| Systemic calcium/phosphate | Raposo 2002; FGF23 extension by Kuro-o/Hypotheses-style models | How do gut, kidney, parathyroid, calcitriol, phosphate, and exchangeable bone interact? | Nonlinear ODEs; optional events and interventions |
| Coupled mineral/FGF23 homeostasis | Gaweda 2010; Granjon et al. 2017; Gaweda et al. 2021 | How does the bone–kidney–parathyroid FGF23 axis stabilize phosphate and calcium, and fail in CKD? | Nonlinear ODE/QSP models with hormone and renal-function modules |
| Well-mixed bone-cell dynamics | Kroll 2000; Komarova 2003; Lemaire 2004 | How do osteoblast/osteoclast recruitment and feedback create remodeling balance, oscillation, or pathology? | ODEs and delay differential equations |
| Signaling/control models | Pivonka 2008/2010; Bellido/Peterson intracellular pathway | How do RANKL/OPG, TGF-beta, PTH, Runx2, CREB, and Bcl-2 regulate cell populations? | ODEs, Hill/mass-action terms, sensitivity analysis |
| Spatial BMU/mechanobiology | Ryser 2009/2010; Buenzli/Pivonka; V-Bone 2020 | How do space, diffusion, mechanics, osteocytes, and BMU geometry shape remodeling? | PDEs, method-of-lines, reaction-diffusion, or agent-based models |
| Integrated calcium/bone QSP | Peterson & Riggs 2010; Peterson et al. 2012; Khurana et al. 2019 | What are the coupled clinical consequences of disease and PTH-based interventions? | Large ODE system with dosing events and observables |
| CKD-MBD QSP | Riggs et al. 2012; Gaweda et al. 2021; later FGF23/vascular models | How does declining renal function produce mineral-bone disease, and which therapy combinations are plausible? | ODEs with disease progression, FGF23, soft tissue, and drug modules |

## Priority replication queue

### Tier 1: canonical and tractable

1. **Kroll (2000), PTH temporal effects.** A delayed cell-differentiation model that directly tests continuous versus intermittent PTH. It is the shortest path to a rigorous pulse-scheduling experiment.
2. **Komarova et al. (2003), osteoclast autocrine regulation.** A compact power-law model with osteoblasts, osteoclasts, and bone mass. It is useful for equilibrium, stability, parameter-sensitivity, and Paget-like behavior.
3. **Lemaire et al. (2004), osteoblast–osteoclast coupling.** A mass-action model with precursor and mature cell populations. It is a strong baseline for hormone, cytokine, estrogen, vitamin-D, and combination-treatment perturbations.
4. **Raposo et al. (2002), minimal calcium homeostasis.** An 11-state organ-level ODE model linking serum calcium/phosphate, PTH, calcitriol, gut, kidney, parathyroid, and exchangeable bone mineral.
5. **Pivonka et al. (2008/2010), RANKL/OPG control.** A mechanistic extension that tests alternative ligand-expression structures and virtual interventions in a BMU-style cell system.

### Tier 2: integrated and clinically oriented

6. **Peterson & Riggs (2010), integrated calcium homeostasis and bone remodeling.** A 28-equation model assembled from Raposo, Lemaire, and an intracellular PTH-response hypothesis. This is the closest published predecessor of the R code in this folder. The public BioModels entry and Metrum's `OpenBoneMin` R/mrgsolve repository are the best source artifacts for a canonical version; the local R code is a later 32-state Khurana-oriented variant.
7. **Peterson et al. (2012), BMD extension.** Adds a nonlinear link from formation/resorption biomarkers to lumbar-spine BMD. It is essential if the project is intended to make bone-density predictions rather than only biomarker predictions.
8. **Khurana et al. (2019), hypoparathyroidism/PTH(1-84) QSP.** An adaptation of the Peterson model with PTH pharmacokinetics, oral calcitriol, and cumulative 24-hour urinary calcium. It is directly relevant to the local code, but its assumptions and dosing labels need auditing.
9. **Lemaire & Cox (2019), drug-combination bone-remodeling model.** Extends the Lemaire lineage with teriparatide, denosumab, and romosozumab PK/PD and calibrates against 57 literature checks. This is a clinically ambitious treatment model, but the full parameter/code artifact is less accessible than the CellML predecessors.
10. **Jörg et al. (2022), multi-drug osteoporosis model.** Open full text with dynamic modules for denosumab, romosozumab, bisphosphonates, and teriparatide. It is a strong replication target for combination therapy and a useful comparator for the 2019 model.
11. **Granjon, Bonny & Edwards (2017), coupled calcium/phosphate model.** Adds FGF23, vitamin D3, CaSR regulation, and calcium–phosphate–fetuin-A complexes to a rat mineral-homeostasis model. It is the cleanest next systemic extension beyond Raposo.
12. **Gaweda et al. (2021), CKD-MBD QSP.** A modern benchmark that builds on an open-source calcium/phosphorus model, estimates against 5,496 CRIC participants, and shows that FGF23 and a soft-tissue compartment are needed to reproduce CKD trajectories.
13. **Ruiz-Lozano et al. (2024) and the 2026 integrated Ca/P–bone-remodeling framework.** This is the most important new bridge for the project: it replaces the simplistic Peterson/Riggs bone compartment with a mechanistic five-cell BCPM, variable PTH, plasma–marrow–matrix mineral exchange, mechanobiology, and a mineralization algorithm. It is the highest-value next replication target after the current ODE suite because it directly addresses the gap between systemic Ca/P homeostasis and cell-level remodeling.

### Tier 3: spatial and mechanobiological extensions

14. **Ryser, Nigam & Komarova (2009) and Ryser, Komarova & Nigam (2010).** Spatial BMU models with diffusion, cell movement, time delays, and RANKL/OPG. These answer questions unavailable to well-mixed ODEs, but require a more explicit numerical validation plan.
15. **Buenzli/Pivonka BMU models.** Spatial refilling and osteocyte-generation models that connect cell distributions, geometry, and formation rates.
16. **V-Bone (2020).** A mechano-biochemical, spatial simulation platform using reaction-diffusion signaling and cell-genesis/apoptosis rules. It is a high-value extension target rather than the first replication target.
17. **Osteocyte-centered models.** The 2013 osteocyte/sclerostin/RANKL model, and newer mechanobiological models, provide a route to mechanical loading and targeted remodeling that is absent from the local R implementation.
18. **Nelson et al. (2026), surgical menopause.** A public MATLAB ODE model with estrogen forcing and BMD output. It is a useful modern benchmark for endpoint validation, although its endocrine and mineral-homeostasis modules are much simpler than Peterson/Khurana. A seven-state Julia translation is now included in this project.
19. **Garzón-Alvarado et al. (2026), remodeling plus fatigue damage.** Couples a Komarova-type BMU model to Miner/S–N damage accumulation under prescribed mechanical loads. This is a useful mechanobiology benchmark, but the authors explicitly frame it as a local, hypothesis-generating model rather than a validated fracture-risk model.
20. **Rattsev et al. (2026 preprint), genetic postmenopausal bone loss.** A 27-ODE model with estrogen, RANK/RANKL/OPG, Wnt, TGF-beta, BMD, drug response, and 22 genetic variants. It is high-value for personalized modeling, but it should remain a watch-list target until the preprint's code/data availability and peer-review status are stable.

### Current landscape update

A recent systems-biology review organizes roughly 80 bone-remodeling models by
power-law, mass-action, spatial/PDE, and agent-based approaches. Its most
useful conclusion for this project is that the field has many local remodeling
models but few validated bridges from cell/signaling dynamics to organ-level
calcium, phosphate, BMD, or fracture outcomes. It also highlights persistent
underrepresentation of reversal cells, bone-lining cells, immune cells, and
mechanical regulation. The review is a ChemRxiv preprint, so it is useful as a
map rather than as a definitive clinical consensus.

There is also a newly useful open-source benchmark: Nelson et al. (2026),
“Mathematical Modeling of Bone Remodeling in Surgical Menopause.” Its GitHub
repository includes MATLAB equations, steady-state initialization, an
estrogen intervention, BMD output, and human spine-BMD data. It is not a
replacement for the mechanistic calcium models above, but it is a valuable
modern validation target because it connects an endocrine perturbation to a
clinical density endpoint with public code.

Two additional 2026 directions are important but should not be conflated with
validated clinical QSP. Garzón-Alvarado et al. add fatigue damage to a
Komarova-type BMU model and explicitly identify the framework as local and
hypothesis-generating. Rattsev et al. extend mechanistic postmenopausal bone
modeling into polygenic virtual subjects and drug-response calibration, but the
work is currently a preprint. Both are good extension targets after the core
Julia models have equation-level regression tests.

## Reproducibility status

| Model | Equation access | Parameter access | Code/accessibility | Replication confidence |
|---|---|---|---|---|
| Kroll 2000 | Published equations and public CellML pulse artifact | Published/CellML tables | Public CellML translation | High for the CellML ODE; separate DDE reconstruction still needed |
| Komarova 2003 | Published equations | Published/secondary tables | Straightforward hand port | High for the core ODE |
| Lemaire 2004 | Published equations and Physiome CellML | CellML parameters | Public CellML translation | High, with extension variants tracked separately |
| Pivonka 2008 | Published equations and Physiome CellML | CellML/Matlab translation | Public CellML translation | High once ligand-expression variant is fixed |
| Raposo 2002 | Published equations plus public CellML artifact | Published/CellML tables | Public Physiome CellML translation | High for the as-coded artifact; low for physiological validity until initial conditions are repaired |
| Peterson & Riggs 2010 | Full equations in paper and local R source | Local R parameter file | Local R implementation plus Julia port | High for “as coded”; medium for “as published” until discrepancies are resolved |
| OpenBoneMin canonical Peterson/Riggs implementation | Public C++/mrgsolve source and documentation | Public model parameters and initial conditions | Public GitHub repository; Julia BMD layer now cross-checked against the published equation | High for the source artifact; exact figure reproduction still requires matching regimen conventions |
| Peterson et al. 2012 | Model extension and BMD equation available in open paper | BMD parameters published; denosumab data digitized in original study | Julia 33-state extension now added; original R site is no longer reliable | High for the equation-level extension; medium for exact figure reproduction |
| Khurana 2019 | Adaptation described in paper and local R source | Local R parameter file plus paper | Local R implementation plus Julia port | High for “as coded”; lower for independent parameter recovery |
| Lemaire & Cox 2019 | Published model description and 57 calibration checks | Article-level values; full artifact not located | No verified public code found | Medium |
| Jörg et al. 2022 | Full open article with appendices | Article/supplement | Open full text; code artifact not yet located | Medium to high after appendix transcription |
| Granjon et al. 2017 | Full equations in open journal article | Article tables/supplement plus CaPO4Sim defaults | Public CaPO4Sim R/C repository plus Julia port | High for the public artifact; strong biological scope, but rat-specific and not a bone-cell remodeling model |
| Gaweda 2021 CKD-MBD | Paper-level QSP description | Article/supplement data need recovery | Open-source predecessor referenced | Medium until source code is obtained |
| Ruiz-Lozano 2024 / integrated Ca/P–BCPM 2026 | Open equations and supplement | Tables plus calibrated disease scenarios | No verified public Julia/C++ artifact found | Medium to high for equation-level replication; high-value next target |
| Ryser/V-Bone spatial models | Equations and numerical method described | Varies by supplement | Public papers/repositories, but higher implementation burden | Medium to high depending on artifact |
| Nelson 2026 | Open paper and public MATLAB repository | Model parameters and BMD data available in repository | Public GitHub code plus Julia seven-state port | High for as-coded; moderate for clinical generalization |
| Garzón-Alvarado 2026 | Open paper and equations | Article-level parameters/scenarios | No stable code artifact located | Medium for a qualitative extension |
| Rattsev et al. 2026 preprint | Open preprint and model description | 50 fixed + 14 calibrated parameters described; data by request | No stable public code artifact located | Medium for conceptual replication; low for independent calibration |

The Physiome Model Repository is especially useful here because it provides machine-readable CellML artifacts for the Raposo 2002, Lemaire 2004, Pivonka 2008, and Kroll 2000 models. Those files will be treated as source artifacts, not as substitutes for checking the original publication. In particular, the Raposo artifact is reproducible but its default state initialization is not physiologically self-consistent: the initial bone-to-plasma calcium flux is 100 while the plasma-to-bone flux is 3.3, so the plasma calcium pool crosses below zero. That is an important source-level finding, not a Julia solver failure.

Granjon et al. is now a particularly strong reproducibility case because the authors' model is available in the public [CaPO4Sim repository](https://github.com/DivadNojnarg/CaPO4Sim). Its compiled C core, R core, fixed-parameter file, and model-engine initial states make it a better replication source than manually transcribing the paper alone. The Julia port follows the compiled core and is derivative-tested against the R core at the public hypoparathyroidism state.

## What each core model can answer

| Model | Scientific question it answers | High-value extension |
|---|---|---|
| Kroll | Why can the same PTH molecule be anabolic when intermittent and catabolic when continuous? | Replace fixed delays with measured transit-time distributions; fit pulse frequency, amplitude, and exposure jointly. |
| Komarova | Which autocrine/paracrine feedbacks can stabilize bone mass or create oscillatory/pathological dynamics? | Add explicit RANKL/OPG, osteocyte signals, stochasticity, or patient-specific turnover markers. |
| Lemaire | How do precursor/mature OB and OC interactions respond to PTH, estrogen loss, vitamin-D deficiency, senescence, and glucocorticoids? | Add pharmacokinetics, uncertainty, and a calibrated mapping to BSAP/CTX/BMD. |
| Raposo | How do endocrine feedbacks distribute a calcium/phosphate perturbation across gut, kidney, parathyroid, and exchangeable bone? | Add FGF23, renal disease progression, soft tissue, acid-base effects, and explicit mass-balance diagnostics. |
| Pivonka | Which RANKL/OPG expression architecture best controls BMU behavior and treatment response? | Add osteocyte-derived RANKL, sclerostin/Wnt, denosumab pharmacology, and spatial localization. |
| Peterson & Riggs | What is the integrated response to hypoparathyroidism, renal insufficiency, continuous/intermittent PTH, and denosumab? | Correct/document uncertain equations, add BMD and site specificity, modernize vitamin-D/FGF23 biology, and quantify identifiability. |
| Peterson BMD | How do short-term formation/resorption changes accumulate into nonlinear BMD trajectories? | Use multi-site BMD, fracture-risk outputs, and longitudinal patient data rather than one lumbar-spine mapping. |
| Khurana | How do PTH(1-84) schedules affect serum calcium and urinary calcium in hypoparathyroidism? | Add titration, adherence, inter-individual variability, alternative vitamin-D analogs, and explicit PK/PD calibration. |
| Lemaire & Cox 2019 | Which combinations and sequences of teriparatide, denosumab, and romosozumab can improve turnover and BMD across osteoporosis phenotypes? | Add patient-level variability, fracture-risk endpoints, treatment safety constraints, and prospective validation. |
| Jörg et al. 2022 | How do multiple osteoporosis drug mechanisms interact over long time scales? | Add calcium/phosphate/PTH coupling, drug resistance or rebound, site-specific BMD, and probabilistic calibration. |
| Granjon et al. 2017 | How does FGF23 coordinate renal, intestinal, skeletal, and intracellular Ca/PO4 fluxes? | Add human/CKD calibration, explicit bone-cell remodeling, soft-tissue calcification, and intervention PK/PD. |
| CKD-MBD QSP | Which combinations of phosphate handling, FGF23, calcitriol, PTH, dialysis, binders, and calcimimetics can control CKD-MBD? | Couple mineral markers to bone turnover, vascular calcification, and patient-level Bayesian calibration. |
| Ruiz-Lozano / integrated Ca/P–BCPM | How do variable PTH, Ca/P availability, RANK/RANKL/OPG, Wnt–sclerostin, osteocytes, mechanics, and matrix mineralization jointly determine bone density in health, CKD, HPT, and vitamin-D deficiency? | Replicate the queue/mineralization algorithm, add FGF23 and explicit clinical PK/PD, and separate local RVE mechanics from whole-skeleton observables. |
| Spatial BMU/V-Bone | How do geometry, diffusion, mechanical stimulus, and osteocyte signaling create a moving remodeling unit? | Couple organ-level endocrine inputs to spatial bone mechanics and a clinically observable BMD/fracture endpoint. |
| Nelson 2026 | How do estrogen trajectories after surgical menopause propagate through bone-cell dynamics into spine-BMD loss? | Add calcium/phosphate/PTH/FGF23 coupling, bone-site heterogeneity, uncertainty, and fracture-risk outputs. |

## Assessment of the local R model

The local code is not an independent new model. It is a modified implementation of the Peterson/Riggs integrated model, adapted toward the Khurana hypoparathyroidism/PTH(1-84) use case. It is valuable as the first end-to-end Julia target because it already contains dosing, calcium/phosphate, bone-cell, and urinary-calcium components. The Julia port now reproduces the local R run's full 210-day headline ranges and final values, including the 10,081-row event-aware output convention.

Before calling it a replication, the Julia port must preserve a frozen “as-coded” version and a separately documented “corrected/physiologic” version. The current R run is finite but produces a negative immediate-exchangeable bone-calcium state (`Q`), which is physically impossible. The code also contains an apparent calcitriol-balance transcription issue, parameter/comment mismatches around PTH clearance and gland loss, and a drug-label mismatch between “teriparatide” comments and full-length PTH(1-84)-like dosing. These are not to be silently fixed: each change should become a regression test and be traceable to the source paper or an explicit modeling decision.

## Replication standard

Every Julia model should ship with:

- equations and state/parameter names in a human-readable file;
- source-specific parameter and initial-condition tables;
- a baseline reproduction test for every published figure that can be digitized;
- dimensional, positivity, mass-balance, and event-order tests;
- solver tolerances and a documented numerical method;
- an “as published” mode and a “repaired/extended” mode;
- global and local sensitivity analysis, with practical identifiability assessed before fitting new parameters;
- intervention scenarios that answer a question, not just a plot;
- a machine-readable provenance record for every parameter and equation.

## Initial interpretation

The strongest scientific opportunity is not merely porting R to Julia. It is creating a reproducible hierarchy in which the short models explain mechanisms, the integrated model connects mechanisms to systemic observables, and the modern CKD-MBD/mechanobiology models test whether those mechanisms survive added biology. The main unresolved gap is the bridge between endocrine mineral homeostasis, spatial/mechanical remodeling, and clinically measurable outcomes. That is where the extensions should converge.

## Primary sources and local source papers

- Raposo et al. 2002: https://pubmed.ncbi.nlm.nih.gov/12213894/
- Raposo 2002 public CellML artifact: https://models.cellml.org/workspace/raposo_sobrinho_ferreira_2002/@@rawfile/3dcc0383f02438af3b98d3e74b025c870ff95bf1/raposo_2002.cellml
- Kroll 2000: https://pubmed.ncbi.nlm.nih.gov/10824426/
- Kroll 2000 public CellML artifact: https://models.physiomeproject.org/workspace/kroll_2000/rawfile/5686223dbf499c1738659fde4c45c9eddbfd0dc2/kroll_2000.cellml
- Komarova et al. 2003: https://pubmed.ncbi.nlm.nih.gov/14499354/
- Lemaire et al. 2004: https://pubmed.ncbi.nlm.nih.gov/15234198/
- Lemaire 2004 public CellML artifact: https://models.physiomeproject.org/exposure/3094a5c3e0810028a9accc3b773cab36/lemaire_tobin_greller_cho_suva_2004.cellml/view
- Pivonka et al. 2008: https://doi.org/10.1016/j.bone.2008.03.025
- Pivonka 2008 public CellML artifact: https://models.physiomeproject.org/workspace/pivonka_zimak_smith_gardiner_dunstan_sims_martin_mundy_2008/rawfile/9848d85152c7ff7bae01d8ed69a9464992d8acf5/pivonka_2008_matlab.cellml
- Pivonka et al. 2010: https://pubmed.ncbi.nlm.nih.gov/19782692/
- Ryser et al. 2009: https://pubmed.ncbi.nlm.nih.gov/19063683/
- Peterson & Riggs 2010: https://pubmed.ncbi.nlm.nih.gov/19732857/
- Peterson et al. 2012: https://pubmed.ncbi.nlm.nih.gov/23835796/
- Khurana et al. 2019: https://pubmed.ncbi.nlm.nih.gov/30350311/
- Lemaire & Cox 2019: https://pubmed.ncbi.nlm.nih.gov/30460589/
- Jörg et al. 2022 osteoporosis multi-drug model: https://pubmed.ncbi.nlm.nih.gov/35942681/
- Gaweda et al. 2010 FGF23/CKD model: https://pmc.ncbi.nlm.nih.gov/articles/PMC2901635/
- Granjon et al. 2017 coupled calcium/phosphate model: https://doi.org/10.1152/ajprenal.00271.2017
- Granjon/CaPO4Sim public implementation: https://github.com/DivadNojnarg/CaPO4Sim
- Gaweda et al. 2021: https://pubmed.ncbi.nlm.nih.gov/33308018/
- V-Bone 2020: https://pubmed.ncbi.nlm.nih.gov/32181336/
- Osteocyte/sclerostin/RANKL model 2013: https://pubmed.ncbi.nlm.nih.gov/23717504/
- Cook and Lighty et al. review preprint (2024): https://doi.org/10.26434/chemrxiv-2024-5vrcc
- Nelson et al. surgical-menopause model and public code (2026): https://github.com/ashleefv/SurgicalMenopauseBone
- Garzón-Alvarado et al. remodeling plus damage (2026): https://link.springer.com/article/10.1007/s10441-026-09522-x
- Rattsev et al. genetic postmenopausal bone-loss preprint (2026): https://doi.org/10.64898/2026.06.04.26354968
- Ruiz-Lozano et al. age/postmenopausal BCPM (2024): https://pubmed.ncbi.nlm.nih.gov/38700787/
- Ruiz-Lozano et al. integrated Ca/P–bone-remodeling framework (2026): https://doi.org/10.3389/fbioe.2026.1800350
- Metrum OpenBoneMin public model repository: https://github.com/metrumresearchgroup/OpenBoneMin
- Peterson–Riggs BioModels entry: https://www.omicsdi.org/dataset/biomodels/BIOMD0000000613

Local PDFs reviewed: `literature/raposo2002.pdf`, `literature/peterson2010.pdf`, and `literature/khurana2018.pdf`.
