# Fail-closed postprocessing gates for the one-system design

Status: frozen method specification; no production result exists yet.

System: experimental 8KCT gamma-secretase with native O6U nirogacestat.

Sampling design: one audited membrane construction and deterministic minimization, followed by three independent thermalization/equilibration branches and three 500 ns production realizations with distinct Maxwell velocity seeds. The raw trajectory, energy, log, checkpoint, run input, topology, index, and final coordinate file from every realization are immutable source data.

The primary analysis interval is fixed at 200-500 ns for every realization; complete 0-500 ns traces remain visible for QC. No outcome-dependent burn-in or extension is permitted.

## Scope boundary

The MD question is limited to physical stability of the constructed membrane-protein system and retention of the experimental nirogacestat pose. The production workflow does not run MM/GBSA, MM/PBSA, energy decomposition, interaction entropy, normal-mode entropy, or a frame-level hypothesis test. It does not estimate an absolute binding free energy.

PCA, occupancy-derived `-kBT ln P` maps, FEL-labelled plots, and 3D population/free-energy surfaces are prohibited for this 3x500 ns protocol. The retained PCA/FEL command names are rejection stubs and cannot produce an analysis.

## Required frozen index groups

Create the following named groups once from the accepted `prod.tpr` and `index.ndx`, then record the group number, atom count, atom-index SHA-256, and selection expression in the analysis manifest:

- `System`: every atom in the run.
- `Complex`: all gamma-secretase protein atoms plus O6U; no lipid, water, or ion.
- `Protein`: all protein atoms.
- `TMCore_CA`: pre-specified C-alpha atoms in experimentally resolved transmembrane helices used only for alignment. The residue list must be frozen before viewing production results.
- `O6U_Heavy`: all non-hydrogen O6U atoms. Its atom count must equal the accepted ligand structure record.
- `Bulk_POPC`: all bulk POPC atoms.
- `Structural_Lipids`: the three retained CLR and two retained PC1/DSPC molecules, with atom-index identities fixed from the accepted build.
- `Membrane`: `Bulk_POPC` plus `Structural_Lipids`.
- `Phosphate`: the named POPC and PC1 phosphate selections used for thickness and leaflet tracking, with species kept distinguishable.

Interactive GROMACS group numbers are not hard-coded. Instead, the preflight script resolves names to group numbers and writes small selection-input files. Every command log must include those files and their hashes.

## Gate 1: raw-run integrity

Run these checks before any coordinate transformation for each realization:

```bash
gmx check -f prod.xtc -s1 prod.tpr > qc/gmx_check_xtc_tpr.txt 2>&1
gmx check -e prod.edr > qc/gmx_check_edr.txt 2>&1
gmx check -c prod.gro > qc/gmx_check_final_gro.txt 2>&1
```

Hard failures:

- a command exits nonzero;
- the trajectory or energy file is unreadable, truncated, or has non-monotonic/duplicate time after concatenation;
- the final time differs from 500 ns by more than one saved-frame interval without a documented clean checkpoint restart after hardware interruption;
- the trajectory and `tpr` are incompatible, including a bond/virtual-site mismatch reported by `gmx check`;
- any coordinate, velocity where stored, box vector, or energy value is NaN or infinite;
- any production log contains `Fatal error`, `Segmentation fault`, `LINCS WARNING`, `constraint error`, `SETTLE`, or an unrecovered domain-decomposition failure;
- restart segments overlap or leave a time gap;
- the recorded `tpr`, topology, MDP, coordinate, or ligand-parameter hash differs from the preflight manifest.

No `-maxwarn` override is permitted in a production `grompp` command.

## Gate 2: energy, cell, and spike classification

Extract at least the following terms with the exact GROMACS build that ran production: temperature, pressure, pressure tensor components, potential energy, kinetic energy, total energy, density, volume, box-X, box-Y, and box-Z. Preserve the unsmoothed values.

Every scalar series is evaluated in the following order:

1. verify finite values and strictly increasing times;
2. plot the raw series at the saved resolution;
3. calculate the first difference and flag a point when its absolute robust z-score exceeds 12, using the median and `1.4826 * MAD` of the first-difference series;
4. inspect the matching log steps, box vectors, and coordinate frames for every flag;
5. compute fixed-duration block summaries only after the raw flags are adjudicated.

The robust-z threshold is an engineering review trigger, not a universal physical cutoff. A flagged point is never automatically removed.

Every structural, membrane, and energy flag must be covered exactly once by the approved schema-v3 adjudication record. That record is bound to the structural, membrane, and energy `COMPLETE.json` hashes, a prospectively frozen policy file and hash, the exact flag identity, and at least one matched evidence file and hash. Missing or duplicate coverage, changed evidence, an unresolved disposition, or any allowed `analysis_blocked` disposition blocks publication eligibility.

Spike decisions are fail-closed:

- A positional RMSD/distance spike that disappears after the validated PBC pipeline is a postprocessing artifact. Record it and use only the corrected trajectory for positional analysis.
- A raw energy spike associated with a constraint warning, NaN, impossible box change, or corrupted frame fails that realization.
- A finite raw energy spike without a log or coordinate anomaly remains in the data. Report its component terms; do not delete, winsorize, interpolate, or replace it.
- An instantaneous pressure peak alone is not a failure. Pressure is judged with its tensor, sustained block behavior, and cell dimensions; no frame is deleted because pressure looks large.
- A ligand-RMSD spike that persists after correct PBC handling and has continuous minimum-image coordinates is a sampled event. It is not repaired or censored. If it is ligand exit, the realization remains part of the result and a global pose-retention claim fails.

## Gate 3: deterministic PBC processing

Use separate GROMACS calls in the order recommended by the version-pinned GROMACS documentation. `sel_*.txt` files below contain the resolved numeric group choices and are archived with the logs.

```bash
# 1. Make all molecules whole.
gmx trjconv \
  -s prod.tpr -f prod.xtc -n index.ndx \
  -o analysis/01_whole.xtc -pbc whole \
  < selections/sel_system.txt \
  > logs/01_whole.log 2>&1

# 2. Cluster protein and O6U if they are separate topology molecules.
gmx trjconv \
  -s prod.tpr -f analysis/01_whole.xtc -n index.ndx \
  -o analysis/02_cluster.xtc -pbc cluster \
  < selections/sel_complex_then_system.txt \
  > logs/02_cluster.log 2>&1

# 3. Extract the processed first frame; it becomes the nojump reference.
gmx trjconv \
  -s prod.tpr -f analysis/02_cluster.xtc -n index.ndx \
  -dump 0 -o analysis/02_first.gro \
  < selections/sel_system.txt \
  > logs/02_first.log 2>&1

# 4. Remove jumps using that processed first frame as the reference.
gmx trjconv \
  -s analysis/02_first.gro -f analysis/02_cluster.xtc -n index.ndx \
  -o analysis/03_nojump.xtc -pbc nojump \
  < selections/sel_system.txt \
  > logs/03_nojump.log 2>&1

# 5. Center on the protein-ligand complex and rebox whole molecules.
gmx trjconv \
  -s prod.tpr -f analysis/03_nojump.xtc -n index.ndx \
  -o analysis/04_centered_system.xtc -center -pbc mol -ur compact \
  < selections/sel_complex_then_system.txt \
  > logs/04_center.log 2>&1

# 6. Create the fitted protein-ligand analysis branch.
# No PBC operation is permitted after this call.
gmx trjconv \
  -s prod.tpr -f analysis/04_centered_system.xtc -n index.ndx \
  -o analysis/05_fitted_complex.xtc -fit rot+trans \
  < selections/sel_tmcore_then_complex.txt \
  > logs/05_fit.log 2>&1
```

If `Complex` is already one whole topology molecule, step 2 may be omitted only after topology inspection and beginning/middle/end visual confirmation. That decision is recorded; it is not guessed from the animation.

PBC validation failures:

- any command exits nonzero or changes the frame count/time stamps;
- O6U is separated from the protein because the wrong cluster group was selected;
- protein bonds are broken in first, middle, or last processed snapshots;
- minimum-image protein-O6U contact distances differ between raw and processed trajectories by more than `0.01 nm` at any matching frame over the complete 0-500 ns trace; this tolerance accommodates repeated compressed-coordinate quantization and is not a biological cutoff;
- any PBC or reboxing operation is applied after fitting;
- only a smoothed/movie trajectory is retained while the raw source is discarded.

Use `04_centered_system.xtc` for membrane and box analyses. Use `05_fitted_complex.xtc` only for protein/ligand structural analyses.

## Gate 4: membrane quality control

All membrane metrics are reported by realization and over the complete trajectory. The frozen 200-500 ns analysis window exactly matches the mature 300 ns interval used in the direct triplicate gamma-secretase precedent and was fixed before any production outcome exists. It supports the narrow pose-compatibility question but is not proof of global conformational convergence. Membrane/core QC determines whether that common window is admissible, not where a realization-specific cutoff should be moved.

Required checks:

- POPC phosphate-to-phosphate bilayer thickness from the two leaflet density maxima;
- protein-aware local/Voronoi area per lipid by species and leaflet; the naive `Lx*Ly/Nleaflet` value is cell QC only and is not reported as the protein-containing membrane's area per lipid;
- POPC acyl-chain deuterium order profiles using the exact CHARMM36 atom mapping;
- stable leaflet membership for phospholipid headgroups and no unexplained leaflet-count change;
- time series of lateral area, box-Z, volume, density, and pressure tensor;
- beginning, middle, and final views showing no vacuum gap, box collapse, persistent lipid void outside the protein boundary, or broken membrane molecule.

Protein-aware area per lipid and POPC deuterium order parameters are hard pre-production NO-GO gates. The approved route must use APL@Voro v3.3 or FATSLiM for the exact post-build lipid identities and protein footprint, and `gorder` for the exact CHARMM36 POPC carbon-hydrogen and unsaturated-chain mapping. The topology identity, protein/POPC atom-index sets, every `gorder` chain/carbon/hydrogen identity, and unsaturated-chain handling are explicit and hash-bound; output profile keys must match that mapping exactly. Commands, versions, mappings, outputs, and hashes must be retained. Periodic lateral cell area is not accepted as area per lipid, and no generic lipid atom-name guess is accepted for order parameters.

Water in the gamma-secretase internal cavity is not automatically a defect because the catalytic region is hydrated. Only an unintended bulk-water defect through the lipid phase outside the protein or a persistent vacuum/lipid void is a failure candidate. Core-water exclusion uses minimum-image three-dimensional water-to-protein distance within the membrane slab and three-dimensional periodic water clustering. A protein atom in the nicastrin ectodomain cannot mask a membrane-core water merely because its XY projection overlaps.

No universal numeric thickness or area-per-lipid value is imposed before construction because the protein footprint, retained structural lipids, temperature, and force field determine the distribution. Failure is based on nonstationary block behavior, structural disruption, or inconsistency among realizations, not on tuning the system to a desired number.

## Gate 5: primary descriptive structural outputs

The following definitions are frozen before production analysis:

- `protein-core RMSD`: C-alpha RMSD of `TMCore_CA` after fitting on `TMCore_CA`;
- `pocket-aligned ligand RMSD`: O6U heavy-atom RMSD to the experimental 8KCT pose after fitting the pre-specified pocket/TM-core atoms; atom mapping must be one-to-one and chirality preserving;
- `ligand COM displacement`: minimum-image displacement of O6U's center of mass from its experimental-pose reference after protein alignment;
- `native contact`: a protein-O6U heavy-atom pair at or below `0.45 nm` in the accepted experimental 8KCT model;
- `native-contact fraction`: fraction of those frozen contact pairs at or below `0.45 nm` in each frame;
- `hydrogen bond`: donor-acceptor distance at or below `0.35 nm` and hydrogen-donor-acceptor angular deviation at or below 30 degrees, with donor/acceptor identities frozen in the contact manifest.

For each realization, provide the raw time series, distribution, autocorrelation estimate, and block summaries. Show all three realization-level summaries. Do not pool frames and report the resulting standard error as if the number of frames were the sample size.

An unfavorable but continuous event is a result, not a technical failure. Before production, freeze a positive duration no longer than 5 ns. Search the full 0-500 ns trace separately and without gap bridging for (i) O6U heavy-atom RMSD above its egress cutoff, (ii) O6U COM displacement above its egress cutoff, and (iii) native-contact fraction below its contact-loss cutoff. Any no-gap event meeting that short duration blocks a general pose-retention conclusion even when the overall 200-500 ns passing-frame fraction remains above threshold.

## Gate 6: stationarity and inconclusive-result rule

For every retained scalar observable, estimate the integrated autocorrelation time after the pre-defined equilibration exclusion. A valid uncertainty block is at least ten times the longest relevant integrated autocorrelation time, and at least five complete blocks must remain per realization. If those conditions cannot be met, report the trace descriptively and do not report a confidence interval for that observable.

Autocorrelation/block counts do not establish stationarity. For every mandatory structural, membrane, external APL, and exact EDR scalar in each realization, also split the fixed 200-500 ns interval into five calendar-time blocks `[200,260)`, `[260,320)`, `[320,380)`, `[380,440)`, and `[440,500]`. Compute block medians, the linear change across 300 ns, first-to-last median shift, maximum adjacent-block shift, and maximum eligible two-sided change-point shift. Normalize absolute effects by `max(1.4826*MAD, frozen metric-specific scale floor)`. Every limit and scale floor must be positive, source-rationalized, and frozen before production review. Failure of any one gate makes that realization nonstationary; no frame is removed and no window is moved.

The 500 ns endpoint and 200-500 ns primary window are fixed. There is no QC-triggered extension or replacement window in this protocol. The planned pose-compatibility result is inconclusive if any realization has fewer than five valid blocks for a quantitative summary, continues to drift in a mandatory membrane/core observable, undergoes an unresolved late non-ligand transition, or otherwise fails the beginning/middle/end and block audit. The affected quantitative analysis and, when relevant, the complete MD claim are omitted; the cutoff is not moved and the best-looking realizations are not selected.

Ligand exit, loss of native contacts, or materially discordant pose-retention behavior is a scientific result. It blocks a general pose-retention conclusion and does not authorize a rescue rerun. A clean checkpoint continuation after a hardware interruption may finish the originally frozen 500 ns trajectory but may not change its physics, endpoint, or analysis interval.

## Prohibited repairs

- deleting an outlier realization;
- deleting high-energy, high-pressure, high-RMSD, or ligand-exit frames solely because they are inconvenient;
- applying smoothing before a failure decision;
- treating a PBC artifact as physical or treating a persistent physical transition as PBC;
- changing the atom selection, alignment group, contact cutoff, equilibration window, or replica inclusion after seeing the result without recording a protocol deviation;
- reporting a single pooled trajectory or pooled-frame confidence interval while hiding realization-specific results;
- labeling an occupancy histogram as an absolute free-energy calculation.
