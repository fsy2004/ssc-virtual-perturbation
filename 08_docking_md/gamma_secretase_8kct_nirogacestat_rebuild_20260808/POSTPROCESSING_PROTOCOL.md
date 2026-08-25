# Prespecified postprocessing overview

The executable pass/fail specification is `POSTPROCESSING_FAIL_CLOSED_GATES.md`. This overview fixes the scientific output set and interpretation before production begins.

## Immutable source data

For `rep01`, `rep02`, and `rep03`, retain the exact topology/includes, MDP and `mdout.mdp`, TPR, GRO, NDX, CPT, EDR, LOG, XTC, stdout/stderr, command transcript, GROMACS/container/GPU metadata, and SHA-256 manifest. Raw trajectories are never overwritten by analysis transformations.

Before coordinate processing, run `gmx check`, confirm TPR/trajectory atom and bond compatibility, verify monotonic times and box vectors, and scan logs for fatal errors, NaN/Inf, constraint warnings, ignored preprocessing warnings, and broken restarts.

The raw-output validator records integrity only and deliberately leaves QC/stationarity as `not_evaluated`. A separate signed report must record `qc_status=pass` and `stationarity_status=pass` for `rep01`, `rep02`, and `rep03`; its SHA-256 must match the frozen analysis plan and it must be supplied through `make_analysis_trajectories.py --production-qc-report`. No analysis trajectory is released from file-integrity checks alone.

## Fixed time policy

- Show complete unsmoothed 0-500 ns traces for all three realizations.
- Calculate primary structural distributions and occupancies over 200-500 ns.
- Do not choose burn-in separately by realization or by outcome.
- If any mandatory QC or stationarity gate fails in any realization, classify the planned result as inconclusive; do not move the cutoff, drop a realization, or add an outcome-dependent extension.
- A ligand exit, loss of native contacts, or unfavorable heterogeneity remains a scientific result and does not authorize a rescue rerun.

## PBC and spike policy

Follow the version-pinned GROMACS order in the fail-closed protocol: make molecules whole; cluster the protein-ligand complex when topology inspection requires it; create a processed first-frame reference; remove jumps; center/rebox; and fit the structural-analysis branch last. No PBC operation is allowed after fitting.

Use the centered system trajectory for membrane/cell measurements and the fitted complex trajectory for pocket-aligned protein/ligand measurements. Preserve every intermediate command, group-selection input, log, hash, and first/middle/last inspection frame.

Every large peak is traced through the raw and intermediate trajectories and matched to coordinates, box vectors, energies, and logs. A peak may be classified as a PBC artifact only when the validated transformation removes it while minimum-image physical geometry remains continuous. Every structural, membrane, and energy flag must appear exactly once in an approved schema-v3 adjudication record bound to all three component `COMPLETE.json` hashes, a prospectively frozen policy hash, and one or more evidence-file hashes. Missing, duplicated, unresolved, or blocking dispositions make the result ineligible. Persistent finite events remain in the data. No smoothing, clipping, deletion, interpolation, winsorization, or replacement is permitted.

## Required analyses

Report every realization separately, using shared axes and unsmoothed traces:

- temperature, pressure tensor, density, potential/total energy, volume, and box area;
- bilayer thickness, protein-aware local/Voronoi POPC area, POPC acyl-chain order, leaflet integrity, water defects, and protein orientation;
- protein and transmembrane-core RMSD plus protein C-alpha RMSF;
- pocket-aligned O6U heavy-atom RMSD to 8KCT;
- O6U center-of-mass displacement;
- frozen native-contact fraction/occupancy and prespecified pocket/hydrogen-bond distances.

Frames are autocorrelated within a realization. Estimate correlation scales and use sufficiently long blocks for within-realization summaries. Independently, every mandatory scalar in every realization must pass the prespecified fixed five-block 200-500 ns stationarity audit: block-median linear change, first-to-last shift, maximum adjacent shift, and maximum eligible change-point shift, each normalized by the larger of `1.4826*MAD` and its frozen metric-specific scale floor. These effect-size gates are frozen before production and are not p-values. The three realization-level estimates may be summarized by median and range; they do not justify normal-theory inference or a frame-level p value.

Protein-aware area per lipid and POPC deuterium order parameters are external hard pre-production gates. Validate APL@Voro v3.3 or FATSLiM against the exact post-build lipid identities and protein footprint, and validate `gorder` against the exact CHARMM36 POPC carbon-hydrogen and unsaturated-chain mapping. The APL mapping JSON must bind the topology-identity hash, bulk-POPC atom-index hash, protein atom-index hash, and protein-footprint method. The `gorder` mapping JSON must bind the topology-identity hash and explicitly enumerate each chain/carbon atom and its hydrogens, including unsaturated-chain geometry. The executable route ingests hash-bound per-saved-frame APL and fixed five-by-60 ns `gorder` profiles for rep01-rep03, verifies exact profile coverage plus source code/version/command/mapping/validation records, and independently rechecks stationarity. Until the tool route itself passes on pre-production test data, production is not authorized; until all production outputs pass, the publication claim remains blocked.

The executable primary bundle is `analyze_primary_structure_mdanalysis.py`, `analyze_membrane_qc_mdanalysis.py`, `gmx_energy_qc.py`, and the independent `validate_primary_postprocessing.py`. Their manifest and mapping templates are under `templates/`. Production energy analysis must extract directly from each manifest-hashed EDR with the exact approved GROMACS version; the existing-XVG parser is restricted to synthetic self-tests.

## Excluded analyses

MM/GBSA, MM/PBSA, per-residue endpoint-energy decomposition, interaction entropy, normal-mode entropy, docking-score comparisons, and any absolute binding-free-energy claim are excluded.

PCA, occupancy-derived `-kBT ln P` maps, FEL-labelled plots, and 3D population/free-energy surfaces are prohibited in this 3x500 ns protocol. `analyze_fel.py` and `validate_analysis_outputs.py` are deliberate rejection stubs with no production path.

## Decision rule

The native-pose compatibility statement is allowed only if raw integrity, membrane behavior, PBC reconstruction, structural stationarity, and realization-level pose/contact reproducibility all pass. Any prespecified no-gap continuous RMSD-egress, COM-egress, or native-contact-loss event over 0-500 ns is a scientific failure even if the primary-window frame fraction would otherwise pass.
