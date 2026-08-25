# Prespecified one-system study design

## Scientific question

Does the experimentally observed neutral nirogacestat pose in 8KCT remain geometrically compatible with its PSEN1 pocket during three reproducible stochastic realizations of one fully specified membrane model?

This is a native-pose retention study. It does not estimate affinity, inhibition potency, efficacy, target engagement, spontaneous association, or a causal Notch/SSc mechanism.

## Single biochemical model

| Item | Frozen choice |
|---|---|
| Experimental scaffold | RCSB 8KCT, biological assembly 1, 2.60 Angstrom cryo-EM |
| Ligand | Native O6U heavy-atom coordinates, exact neutral CCD stereoisomer, formal charge 0 |
| Resolved cofactors | Retain 3 CLR and 2 PC1/DSPC molecules |
| Glycans | Retain all resolved 18 NAG and 3 BMA residues at 12 nicastrin N-glycosylation sites |
| Bulk membrane | Symmetric pure POPC; resolved CLR/PC1 are counted separately as structural lipids |
| Solvent and ions | CHARMM-modified TIP3P, 0.15 M NaCl |
| Force-field family | CHARMM36m protein, matching CHARMM36 lipid/carbohydrate, validated matching CGenFF ligand parameters |
| Fixed protonation | Neutral O6U; PSEN1 Asp257 deprotonated and Asp385 protonated |
| Thermodynamic state | 310.15 K, 1 bar, semi-isotropic production coupling |
| Integration | 2 fs, PME, bonds to hydrogen constrained, no HMR |

Pure POPC and a symmetric NaCl bath are controlled modeling choices with direct full-complex gamma-secretase precedent; they are not asserted to be the unique physiological membrane or ionic environment.

## Realization hierarchy

There is one membrane construction and one deterministic energy minimization. Replication begins at the first stochastic thermalization step:

1. freeze and hash the accepted minimized coordinates, topology, parameters, index, and MDP inputs;
2. create `rep01`, `rep02`, and `rep03` from those identical minimized coordinates;
3. generate Maxwell velocities with three distinct recorded seeds in the first NVT step;
4. run every subsequent NVT/NPT equilibration stage independently in each branch;
5. start a 500 ns production trajectory from each branch's own equilibrated endpoint.

Production release is staged. Each branch's original 500-ns TPR is created once, then used for an exact 5.0-ns (2,500,000-step) technical canary. All three canaries must pass raw file, checkpoint/restart, energy/trajectory/log, PBC-tool, output-rate, and storage gates before checkpoint continuation is released. Canary frames remain the first 5 ns of the corresponding 500-ns trajectory. The canary is never used to select the most favorable ligand behavior or replace a realization.

The three realizations are independent conditional on one fixed preparation. They quantify stochastic trajectory reproducibility, not membrane-build uncertainty or biological variability. They must never be described as three independently constructed systems.

## Prespecified time windows

- Raw/QC display: 0-500 ns for every realization.
- Primary structural summaries: 200-500 ns for every realization.
- No realization-specific burn-in, frame exclusion, or cutoff change is permitted.
- If any mandatory integrity, membrane, PBC, stationarity, or correlation-aware block gate fails in any realization, the planned analysis is inconclusive and the MD claim is omitted. The cutoff is not moved and no realization is dropped or selectively extended.
- A continuous ligand exit, contact loss, or between-realization difference is a scientific observation and cannot trigger a rescue rerun.

Hardware interruption may be resumed from an exact compatible checkpoint. A setup error that changes chemistry, topology, protonation, parameters, or MDP physics invalidates every downstream branch and requires regeneration from the last approved common checkpoint.

## Outcomes

Primary descriptive outcomes are calculated after fitting a frozen, experimentally resolved PSEN1 transmembrane/pocket atom selection:

- O6U heavy-atom RMSD to the deposited 8KCT pose;
- O6U center-of-mass displacement;
- fraction and occupancy of a frozen native heavy-atom contact set;
- prespecified K380/L432 hydrogen-bond geometry and other literature-reported pocket distances;
- whole-protein C-alpha and transmembrane-core C-alpha RMSD plus whole-protein C-alpha RMSF as structural context.

Required physical-integrity outcomes include energy/temperature/pressure/cell behavior, membrane thickness, protein-aware local area per lipid, POPC order parameters, leaflet integrity, water defects, protein orientation, and topology/trajectory compatibility.

Every realization is shown separately. Frames are autocorrelated samples within a realization, not independent replicates. Report realization-level estimates and their median/range; do not create frame-level p values or naive bootstrap confidence intervals.

## Prohibited analyses and comparisons

- apo, ligand-deleted, semagacestat, or second-GSI production arms;
- redocking as an MD start or docking score as an affinity estimate;
- MM/GBSA, MM/PBSA, per-residue endpoint-energy decomposition, interaction entropy, or normal-mode entropy;
- selecting the most stable realization or deleting a ligand-exit realization;
- smoothing, interpolation, winsorization, peak clipping, or undisclosed frame removal;
- treating an occupancy-derived landscape as an absolute binding free-energy surface;
- running PCA, an occupancy-derived `-kBT ln P` map, an FEL-labelled plot, or a 3D population/free-energy surface in this 3x500 ns protocol;
- changing protonation, contacts, fitting atoms, membrane composition, or the analysis window after viewing outcomes.

## Reporting rule

Use the native-pose compatibility statement only if all three realizations pass the raw-integrity, membrane, PBC, stationarity, and native-pose reproducibility gates. Otherwise report the observed instability or heterogeneity, or classify the planned result as inconclusive when a technical/QC gate fails, without changing the window or realization set.
