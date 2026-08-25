# SSc Notch-HES1 fibroblast programme and nirogacestat repurposing — reproducibility code

**Analysis code repository for the manuscript:** *Single-cell regulatory analysis and exploratory spatial mapping motivate nirogacestat testing of Notch-HES1 signalling in systemic sclerosis* (Manuscript in preparation; repository remains private until publication). Analysis code only; all input data are public (accessions listed below).

---

## Overview

This repository contains the reproducible analysis code for a patient-tissue study that connects the Notch transcriptional effector **HES1** to the **myofibroblast** state in systemic sclerosis (SSc) skin, and that nominates the FDA-approved γ-secretase inhibitor **nirogacestat** for direct pathway testing.

The study integrates:

1. **A large patient-tissue SSc skin atlas** (325,636 cells from 230 donors, five non-redundant cohorts), with donor-level composition and fibroblast-state resolution.
2. **Complete virtual-knockout screens** (CellOracle, skin-scATAC and promoter priors) that supply **descriptive perturbation rankings** — HES1 ranks 28/107 and 30/109 and is not a top-ranked driver; its translational priority arises from convergent Notch biology, donor-level fate association, spatial context and a druggable upstream mechanism.
3. **Donor-level regulon activity and myofibroblast-fate trajectories** (CollecTRI/ULM + CellRank2), including within-donor association (median Spearman ρ = 0.282, 98.1% positive) and a non-significant adjusted SSc–HC contrast (BH q = 0.061).
4. **Cross-cohort machine-learning transfer** across eight prespecified adult held-out cohorts.
5. **SCENIC+ enhancer-driven eGRN** (266 TFs, 502 direct / 285 extended eRegulons, 126 HES1 targets), with donor-pseudobulk accessibility.
6. **Multi-platform spatial analyses** (Visium and Xenium), with composition, depth and domain adjustment and block-restricted null models.
7. **Exploratory cell–cell communication and contact analyses** (NicheNet, LIANA+, COMMOT, segmentation-mask contact), reported as non-significant and hypothesis-generating.
8. **Structural layer**: the deposited 2.60 Å γ-secretase–nirogacestat complex (PDB 8KCT) and an explicit-membrane molecular dynamics simulation, reported as pocket-compatibility evidence without affinity, potency or efficacy claims.

**Evidence boundary:** the study is computational, hypothesis-generating; it makes no causal or therapeutic-efficacy claims. The spatial communication/contact layers are exploratory and non-significant; no MM-GBSA/MM-PBSA or docking-score-as-affinity estimate was used.

---

## Pipeline / stages

| # | Stage | Folder | Key tools |
|---|-------|--------|-----------|
| 1 | scRNA integration + annotation | `01_scrna_atlas/` | Scanpy, Harmony, emptyDrops, DecontX |
| 2 | Virtual-KO engine + skin base-GRN | `02_celloracle_ko/` | CellOracle, Cicero |
| 3 | Myofibroblast-fate trajectory | `03_trajectory/` | CellRank2 |
| 4 | Bulk regulon activity + ML transfer | `04_bulk_regulon_ml/` | decoupler (CollecTRI/ULM), scikit-learn (LOCO), IOBR |
| 5 | Multi-platform spatial co-localisation | `05_spatial/` | squidpy/Scanpy, RCTD, BANKSY |
| 6 | Niche cell–cell communication | `06_ccc_nichenet/` | NicheNet, LIANA+, COMMOT |
| 7 | Complementary eGRN model | `07_scenicplus_egrn/` | SCENIC+, pycisTopic, Snakemake |
| 8 | Structure curation, docking QA + explicit-membrane MD | `08_docking_md/` | CHARMM-GUI, GROMACS 2025.2, MDAnalysis |
| 9 | Publication figures | `09_figures/` | matplotlib, seaborn |

Run order within a stage is encoded in file-number prefixes. The MD subproject lives under `08_docking_md/gamma_secretase_8kct_nirogacestat_rebuild_20260808/` with its own `scripts/`, `config/`, `mdp/` and protocol documents.

---

## Data availability (all public)

- **scRNA/scATAC**: GEO GSE249279, GSE292979, GSE195452, GSE138669, GSE236111, GSE312129, GSE320020; cross-cohort validation and sensitivity sets: GSE130955, GSE58095, GSE249550, GSE181549, GSE9285, GSE32413, GSE76807, GSE231692, GSE95065, GSE125362, GSE76885, GSE264508, GSE288490.
- **Spatial**: Visium and Xenium datasets identified in the manuscript Methods (GSE312932 for Xenium).
- **Structure**: PDB 8KCT (γ-secretase–nirogacestat, ligand O6U).
- Raw data are downloaded separately and are not included in this repository.

---

## Reproducibility notes

- Frozen analysis configuration and thresholds are documented in `08_docking_md/.../ANALYSIS_CONFIG_FREEZE_20260821.md`, `STUDY_DESIGN.md`, `LITERATURE_AND_FAILURE_GATES.md` and the MD `config/` manifests.
- MD production protocol and failed-closed gates: `08_docking_md/.../config/production_protocol_hmr4fs_303K_v1.json`, `POSTPROCESSING_FAIL_CLOSED_GATES.md`.
- O6U ligand parameterization evidence chain and its bounds are summarised in `08_docking_md/.../O6U_PARAMETERIZATION_FINAL_AUDIT_20260823.md`.

---

## License

See `LICENSE`.
