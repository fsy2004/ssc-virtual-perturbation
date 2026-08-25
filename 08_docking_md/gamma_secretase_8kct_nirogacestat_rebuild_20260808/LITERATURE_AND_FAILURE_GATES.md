# Literature basis and pre-production failure gates

Evidence freeze: 2026-08-09. Downloaded full texts, supplements, official records, and hashes are indexed under `references/`.

## Why this single model

- 8KCT is the 2.60 Angstrom cryo-EM structure of human gamma-secretase bound to nirogacestat. It supplies an experimental heavy-atom pose, so docking is unnecessary for the production start: [RCSB 8KCT](https://www.rcsb.org/structure/8KCT) and [Guo et al. 2025](https://doi.org/10.1038/s41594-024-01439-8).
- The deposited model contains O6U, 3 CLR, 2 PC1/DSPC, 18 NAG, and 3 BMA. Retaining resolved components is more evidence-led than stripping them to imitate older lower-resolution models.
- Dehury et al. used one full-complex pure-POPC gamma-secretase system and three 500 ns trajectories with different initial velocity seeds, with the last 300 ns used for mature summaries: [Dehury et al. 2019](https://doi.org/10.1039/C9RA02623A).
- A second triplicate full-complex study also used 500 ns runs and excluded 0-200 ns from primary trajectory analyses because the early segment of a large system can retain starting-configuration artifacts: [Dehury et al. 2020](https://doi.org/10.1039/D0RA04683C). The present design therefore adopts the same 200-500 ns primary window while retaining a narrow local native-pose compatibility claim.
- Other full-complex studies used pure POPC with 0.15 M NaCl and confirm that POPC is an established simplified membrane model: [Aguayo-Ortiz et al. 2017](https://doi.org/10.1039/C7SC00980A).
- POPC/cholesterol mixtures are relevant to explicit lipid-composition questions, but 60:40 was one studied arm rather than a uniquely validated physiological default: [Aguayo-Ortiz et al. 2018](https://doi.org/10.1039/C8CP04138E). Because the present question is ligand-pose compatibility and only one build is allowed, bulk pure POPC avoids an unnecessary cholesterol-placement variable while preserving the three experimentally resolved cholesterol molecules.
- The Communications Biology checklist asks for at least three independent simulations and convergence/statistical analysis. It also emphasizes that a well-sampled simpler model is preferable to a more complex poorly sampled model: [MD reliability checklist](https://doi.org/10.1038/s42003-023-04653-0).

Pure POPC, symmetric 0.15 M NaCl, fixed dyad protonation, and one membrane construction remain model choices. Literature precedent does not make them uniquely physiological.

## Chemistry gates

- O6U must match the official neutral `(2S,2S)` component, InChIKey `VFCRKLWBYMDAED-REWPJTCUSA-N`, formal charge 0, and the deposited 8KCT heavy-atom pose: [RCSB O6U](https://www.rcsb.org/ligand/O6U).
- The previous +2 topology and unvalidated CGenFF penalties of 35.5/32.824 cannot enter the rebuild. Scores in the 10-50 range indicate analogy requiring validation; they are not automatically proof of a wrong parameter. Every high-penalty charge/bonded term must have retained QM target data and an approved disposition: [CGenFF methodology](https://doi.org/10.1002/jcc.21367) and [FFParam workflow](https://doi.org/10.1021/acs.jpcb.4c01314).
- CHARMM and GROMACS representations of the approved ligand must pass charge, atom-order, stereochemistry, and single-point energy regression before membrane construction.
- The fixed PSEN1 Asp257-deprotonated/Asp385-protonated dyad follows the structure-specific 8KCT PROPKA assignment and published monoprotonated gamma-secretase models. It remains a predeclared approximation: cryo-EM does not establish hydrogen positions, and context-dependent simulations do not support calling any one state universally physiological.

## Docking protocol QA

- AutoDock Vina 1.2 introduced the current extensible docking implementation; the official Meeko self-docking tutorial defines a native-ligand envelope plus 5 Angstrom padding workflow: [Vina 1.2](https://doi.org/10.1021/acs.jcim.1c00203) and [Meeko tutorial](https://meeko.readthedocs.io/en/develop/tutorial1.html).
- The frozen three-seed 8KCT/O6U self-redocking experiment found native-like poses in all three searches, but only one of three highest-ranked poses was within the preregistered 2 Angstrom threshold. The protocol therefore failed its ranking gate and is excluded from the manuscript evidence chain. `DOCKING_QA_REPORT.md` and the all-pose CSV retain the complete result.
- No redocked pose enters MD. The experimentally deposited O6U pose is used directly. Docking scores are not free energies, and a lower-ranked favorable pose is not selected after viewing the result.

## Build and replication gates

- Use one manually approved PDB Reader job and one Quick Bilayer job through the [official API](https://www.charmm-gui.org/?doc=api&module=quickb).
- The API must retain heteroatoms, preserve all signed structure decisions, and return a complete, visually inspected, hash-verified GROMACS package.
- Replication begins at first NVT velocity generation after a common deterministic minimization. Three seeds must generate three complete independent dynamic equilibration histories. Reseeding only at production after one shared equilibration is prohibited.
- `grompp -maxwarn 0` is mandatory. No warning is ignored to make a TPR.

## PBC, membrane, and numerical gates

- Current GROMACS guidance requires ordered, separate coordinate transformations and warns against applying PBC operations after fitting: [`gmx trjconv`](https://manual.gromacs.org/current/onlinehelp/gmx-trjconv.html) and [trajectory terminology](https://manual.gromacs.org/current/user-guide/terminology.html).
- `whole`, conditional `cluster`, processed first-frame `nojump`, centering/reboxing, and final fitting are logged and inspected independently.
- A large plotted peak is a review trigger, not an automatic outlier. It is retained unless it is demonstrably a visualization/PBC artifact; even then the raw frame remains immutable.
- Fatal errors, nonfinite values, unexplained constraint warnings, corrupt restarts, TPR/trajectory incompatibility, nonphysical cell collapse, membrane rupture, or unreconstructable topology are hard failures.
- Membrane assessment uses time-dependent thickness, local area, order, leaflet integrity, water defects, pressure tensor, and cell dimensions. Whole-system RMSD alone cannot establish equilibration: [lipid membrane best practices](https://doi.org/10.33011/livecoms.1.1.5966).
- Protein-aware area per lipid must use a validated Voronoi/protein-footprint route such as [APL@Voro](https://doi.org/10.1021/ci400172g) or FATSLiM; POPC deuterium order parameters must use a validated exact CHARMM36 mapping such as [`gorder`](https://doi.org/10.1016/j.softx.2025.102254). Both are hard pre-production gates, not analyses to improvise after the trajectories exist.

## Analysis and claim gates

- Fixed primary window: 200-500 ns; full 0-500 ns traces remain visible. Matching a published duration does not by itself establish global convergence, so the claim remains limited to native-pose compatibility under the modeled conditions.
- There is no outcome-dependent replacement window. Failure of any mandatory gate in any realization makes the planned result inconclusive; the cutoff and realization set remain fixed.
- All three realization-level results must be shown. Autocorrelated frames cannot be converted into independent sample size.
- MM/GBSA and MM/PBSA are omitted. With no comparator or calibration they do not answer the narrow question, and frame-level intervals would be pseudoreplicated. The available membrane-PB settings are system-specific rather than a transferable validation of absolute affinity.
- PCA, occupancy-derived FEL, and 3D population/free-energy surfaces are hard-prohibited in this 3x500 ns protocol; the retained legacy executables are rejection stubs only.
- PLIP supplies rule-based putative contacts, not interaction forces or energies. Structural figures must identify the native 8KCT pose and distinguish literature-reported contacts from computed geometric annotations: [PLIP](https://doi.org/10.1093/nar/gkv315).

## Fail-closed matrix

| Gate | Passing rule | Failure response |
|---|---|---|
| Input identity | RCSB/version/component/stereochemistry/charge/hashes agree | Correct the source and restart curation |
| Protein topology | PSEN1 NTF/CTF remain separate; no silent gap filling or terminal patch | Correct PDB Reader decisions and rebuild |
| Native components | O6U=1, CLR=3, PC1=2, NAG=18, BMA=3 with all expected links | Reject the build |
| Ligand parameters | Neutral non-all-zero charges; every high-penalty term validated; CHARMM/GROMACS regression passes | Block all membrane/MD work |
| Membrane packing | Pure bulk POPC, adequate image separation, no ring penetration/void/duplicate structural lipids | Rebuild the single common system under a new manifest version |
| Replica branching | Three unique seeds at first NVT; complete dynamic equilibration independently | Regenerate all branches from the common minimized state |
| Preprocessing | Exact inputs and zero ignored `grompp` warnings | Resolve cause; never use `-maxwarn` |
| Raw integrity | Complete 500 ns, compatible TPR/XTC/EDR, monotonic times, finite values, valid restarts | Resume a clean hardware interruption or correct the common cause and regenerate affected work |
| Numerical stability | No unexplained LINCS/SETTLE/SHAKE warning, fatal error, or cell/membrane failure | Diagnose; no selective frame deletion |
| PBC reconstruction | Continuous minimum-image pocket geometry and inspected whole complex | Fix processing; if impossible, do not analyze |
| Stationarity | Frozen QC/block rules pass in all three realizations over the fixed window | Classify the planned result as inconclusive; do not move the cutoff or select realizations |
| Ligand behavior | Every realization reported without rescue selection | Treat exit/contact loss/heterogeneity as the result |
| Endpoint energies | Not executed | Remove any output or manuscript claim |
| PCA/FEL/3D population or free-energy surface | Prohibited; no execution gate exists | Do not run or report the map |
| Figure provenance | Native source, contact source, camera, selections, logs, and hashes complete | Regenerate automatically |

No protocol can make a model immune to criticism. These gates prevent the avoidable failures that invalidated the retired analysis and keep the conclusion no broader than the calculation.
