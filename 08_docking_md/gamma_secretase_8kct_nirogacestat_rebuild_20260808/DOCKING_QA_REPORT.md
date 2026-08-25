# Native O6U self-redocking QA report

Status: **FAILED AS A RANKING PROTOCOL; RETIRED FROM THE MANUSCRIPT EVIDENCE CHAIN**

This calculation is a prospective self-redocking quality-control experiment. It is not the source of the MD start pose and is not an affinity, potency, stability, or mechanism calculation. MD starts from the deposited 2.60 Angstrom cryo-EM O6U heavy-atom pose in 8KCT.

## Frozen protocol

- Receptor: hydrogenated 8KCT protein with the then-preregistered double-deprotonated PSEN1 Asp257/Asp385 docking-QA approximation and no O6U in the receptor. This is retained as exact provenance of the failed docking experiment and is not the subsequently frozen production-MD dyad state.
- Ligand: exact neutral O6U stereoisomer reconstructed from the RCSB CCD and deposited heavy-atom geometry; formal charge 0.
- PDBQT preparation: Meeko 0.7.1. The ligand has 38 PDBQT atoms, no zero partial charge, and a rounded partial-charge sum of -0.003 e. The small residual is PDBQT decimal rounding, not the formal chemical charge.
- Docking: AutoDock Vina 1.2.7, `vina` scoring, exhaustiveness 32, 20 modes, energy range 5 kcal/mol, CPU 4, seeds 11111/22222/33333.
- Search box: the native O6U heavy-atom envelope plus 5 Angstrom padding on every side, following the official Meeko self-docking tutorial. Center `(164.105, 174.6835, 142.9475)` Angstrom; size `(19.594, 18.381, 18.711)` Angstrom.
- Pose metric: RDKit symmetry-corrected heavy-atom `GetBestRMS` against deposited O6U. The success threshold was fixed at 2.0 Angstrom before inspecting results.
- Scoring success: the highest-ranked pose is at most 2.0 Angstrom in every seed.
- Sampling success: at least one of 20 poses is at most 2.0 Angstrom in every seed.

## Complete seed-level result

| Seed | Top-ranked Vina score (kcal/mol) | Top-ranked RMSD (Angstrom) | Best sampled rank | Best sampled RMSD (Angstrom) | Top-rank pass | Sampling pass |
|---:|---:|---:|---:|---:|:---:|:---:|
| 11111 | -8.746 | 4.366846 | 5 | 1.090288 | No | Yes |
| 22222 | -8.540 | 4.421699 | 4 | 1.172755 | No | Yes |
| 33333 | -7.923 | 1.081226 | 1 | 1.081226 | Yes | Yes |

All-seed sampling passed, but all-seed top-rank scoring failed because only one of three highest-ranked poses reproduced the native geometry. The protocol can find native-like poses but does not rank them reproducibly. Selecting ranks 4 or 5 from favorable runs would be outcome-dependent cherry-picking.

## Decision

1. Do not report a docking score as binding energy or use docking to confirm affinity.
2. Do not select a lower-ranked pose for MD. Use the experimental 8KCT/O6U coordinates.
3. Do not place the colorful overlay in the biological-results figure. It is retained only as a method-failure QA artifact.
4. The manuscript structural panel will show native 8KCT/O6U. PLIP may annotate only rule-based putative geometric contacts on the prepared native structure; it will not display forces, energies, strengths, or occupancies.

## Reproducibility artifacts

- `docking_native_redock/runs/formal_meeko/redocking_run_manifest.json`: exact commands, versions, hashes, charge audit, box, seeds, and run outputs.
- `docking_native_redock/runs/formal_meeko/redocking_all_poses.csv`: all 60 poses without favorable-pose filtering.
- `docking_native_redock/runs/formal_meeko/redocking_analysis_report.json`: preregistered top-rank and sampling decisions.
- `docking_native_redock/figures/redocking_qa/`: deterministic 600 dpi PyMOL QA overlay, session, camera, generated PML, logs, and hashes.

Method sources: [AutoDock Vina 1.2](https://doi.org/10.1021/acs.jcim.1c00203), [original AutoDock Vina validation](https://doi.org/10.1002/jcc.21334), [official Meeko docking tutorial](https://meeko.readthedocs.io/en/develop/tutorial1.html), [PLIP](https://doi.org/10.1093/nar/gkaf361), and [8KCT nirogacestat structure](https://doi.org/10.1038/s41594-024-01439-8).
