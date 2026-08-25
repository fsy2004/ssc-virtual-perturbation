# Server-rental go/no-go checklist

Current status: **NO-GO**. Do not authorize full production until every item below is complete and `python scripts/validate_preflight.py --manifest config/study_manifest.json --stage production --strict` exits zero. Earlier strict stages use the same explicit `--manifest` argument and their staged status; they do not require unresolved future-stage artifacts.

## Structure and chemistry

- [ ] The downloaded 8KCT/O6U records, wwPDB validation report, versions, dates, and SHA-256 values match the manifest.
- [ ] `O6U_PREPARAMETERIZATION_AUDIT.json` still passes and is bound to the official neutral 76-atom SDF, O6U CCD, 8KCT mmCIF, and the immutable 76-row preparameterization mapping.
- [ ] `STRUCTURE_PRECURATION_AUDIT.json` still passes and the signed PDB Reader checklist is bound to it.
- [ ] The signed structure record contains all segments, artificial/authentic termini, four disulfides, 12 protein-to-NAG links, nine glycan-internal links, 18 NAG, 3 BMA, 3 CLR, 2 PC1, O6U, and every PDB Reader decision; its curated coordinates have not already been PPM-transformed.
- [ ] PSEN1 residues 292-376 were not built, and no bond crosses the NTF/CTF break.
- [ ] O6U is the exact neutral CCD stereoisomer; all partial charges are finite/nonzero and sum to the accepted total charge.
- [ ] Every initial CGenFF penalty above 10 has retained QM/MM target data and a signed disposition; every acceptance criterion in `LIGAND_PARAMETERIZATION_PROTOCOL.md` passes.
- [ ] CHARMM and GROMACS ligand atom order, stereochemistry, charge, and single-point energy regression pass.
- [ ] The fixed O6U/Asp257/Asp385 protonation record is signed and no alternate state is present in an input directory.

## CHARMM-GUI common build

- [ ] One manual PDB Reader job and its report/screenshots/hashes are approved.
- [ ] One Quick Bilayer response/archive exists for symmetric pure POPC, 0.15 M NaCl, the accepted margins, and retained heteroatoms.
- [ ] Quick Bilayer applied PPM exactly once, and the build-level orientation record retains the transformation matrix, pre/post hashes, membrane boundaries, screenshots, job ID, and its own SHA-256.
- [ ] Token, password, cookie, email, SSH secret, and historical plaintext credential material are absent from the package and logs.
- [ ] Archive type/path safety, generator/force-field versions, all topology includes, ligand parameter hashes, component counts, total charge, ions, water, box, and staged equilibration files pass validation.
- [ ] A complete sorted SHA-256 manifest covers every extracted GROMACS input file; topology, both indexes, starting coordinates, production MDP, all six unchanged CHARMM-GUI equilibration MDPs, archive, jobs, and orientation record are cross-bound to the strict build report.
- [ ] Visual audit confirms no ring penetration, unresolved steric clash, native-lipid duplication, vacuum gap, lipid void, broken glycan, pocket disruption, or protein-periodic-image contact; at least 2-3 lipid shells separate periodic protein images.

## Three stochastic branches

The MD-runtime checks in the next section are upstream gates: they must already be sealed and bound before the common minimization or any equilibration command is released.

- [ ] One deterministic minimization endpoint and all upstream hashes are frozen.
- [ ] `rep01`, `rep02`, and `rep03` start from that endpoint, use three distinct recorded seeds at first NVT velocity generation, and run every dynamic equilibration stage independently.
- [ ] No branch inherits velocities or a fully thermalized/equilibrated coordinate set from another branch.
- [ ] Every equilibration MDP, TPR, log, final coordinate, checkpoint, energy file, and hash is present; temperature, density, cell, membrane, and constraint checks pass before production.

## MD runtime

- [ ] Exact OS/container digest, GROMACS build, compiler, CUDA, GPU, driver, FFT library, and command line are recorded; a bundled non-scientific benchmark passes.
- [ ] Production physics are frozen: 310.15 K, 1 bar, 2 fs, no HMR, PME, matching CHARMM cutoffs/dispersion treatment, constraints to hydrogen, semi-isotropic Parrinello-Rahman coupling, COM-removal policy, output cadence, and 500 ns endpoint.
- [ ] `grompp` succeeds with zero warnings and no `-maxwarn` override for all three branches.
- [ ] Checkpoint/restart is tested without overlapping/duplicate time, and the runner cannot silently skip a realization.
- [ ] Each original 500-ns TPR is created exactly once after all-three equilibration PASS; no resume path reruns `grompp` or replaces `production.tpr`.
- [ ] Each realization runs exactly 5.0 ns (2,500,000 steps) as a mandatory same-TPR canary; XTC, EDR, log, and CPT must independently report the exact endpoint, while the final GRO hash must match the latest successful mdrun record. The measured storage projection plus 30 percent headroom must fit both the frozen budget and live filesystem free space. The 5 ns frames remain the start of production and ligand behavior is not a realization-selection criterion.
- [ ] Continuation uses the exact per-realization TPR hash plus `production.cpt -append`; hash-chained segment command/stdout/stderr records are append-only and preserve every interruption/resume.

## Analysis readiness

- [ ] Frozen atom-index hashes exist for System, Complex, Protein, TM-core alignment, pocket alignment, O6U heavy atoms, POPC, structural CLR/PC1, and leaflet/phosphate groups.
- [ ] Native contact pairs, hydrogen-bond geometry, membrane definitions, fixed 200-500 ns window, block rules, and the fail-closed no-extension policy are hashed.
- [ ] Native-pose and build-specific thermodynamic acceptance cutoffs, their units, sources, rationale, approval, and SHA-256 records are frozen before production; no field remains `TODO`.
- [ ] Protein-aware POPC area per lipid has a validated APL@Voro v3.3 or FATSLiM route for the exact post-build lipid identities and protein footprint; lateral cell area is not substituted.
- [ ] POPC deuterium order parameters have a validated `gorder` route using the exact post-build CHARMM36 carbon-hydrogen and unsaturated-chain mapping.
- [ ] Output validation checks raw logs, TPR/XTC/EDR compatibility, restarts, energy/cell behavior, membrane QC, PBC continuity, final time, and all input hashes.
- [ ] A separate sealed QC/stationarity report covers rep01--rep03 and binds the manifest, normalized analysis-plan contract, raw-output report, build report, archive, criterion artifacts, evaluator, approver, and validator hash. Its final file SHA-256 is bound in the analysis plan before `make_analysis_trajectories.py` can run.
- [ ] Matching-frame raw and centered/reboxed trajectories independently reproduce the protein-O6U heavy-atom minimum-image distance within 0.01 nm over the complete 0-500 ns trace in every realization; the XVG series, commands, input/output hashes, and PASS/FAIL report are retained.
- [ ] Each CHARMM-GUI equilibration stage has exactly one immutable grompp/TPR record. Interrupted mdrun stages resume only from their own retained CPT with append; a partial stage without a valid checkpoint is quarantined instead of overwritten.
- [ ] Raw, centered-system, and fitted-complex trajectories use separate immutable paths; no cleanup or plotting code can delete or replace raw frames.
- [ ] Plotting shows all three realizations on shared scales and has no smoothing, peak omission, interpolation, winsorization, or hidden inclusion/exclusion path.
- [ ] MM/GBSA, MM/PBSA, decomposition, docking-score analysis, and primary PCA/FEL commands are absent from the production run list.
- [ ] PCA, occupancy-derived FEL, and 3D population/free-energy surfaces are hard-prohibited with no supplementary exception; the legacy executables are rejection stubs.

## Structural figure readiness

- [ ] The figure starts from native 8KCT/O6U, never a redocked pose.
- [ ] Literature-reported contacts and PLIP-derived geometric contacts have separate provenance fields.
- [ ] PyMOL input, contact record, selections, camera, colors, versions, PML/PSE, 600 dpi master PNG, logs, and SHA-256 manifest regenerate deterministically with no hand-added interaction.
- [ ] The legend calls PLIP annotations putative geometric contacts and makes no force, energy, affinity, or interaction-strength claim.

## Operations and storage

- [ ] Historical credentials have been rotated and are not reused.
- [ ] A measured short-segment output rate is used to estimate 3 x 500 ns raw/derived storage, checkpoints, and backups; at least 30 percent working headroom remains after one complete off-server copy.
- [ ] Periodic checkpoint sync, immutable SHA-256 manifests, job-failure notification, and off-server backup restoration are tested.
- [ ] Production commands and environment can be reproduced from a clean server directory without files outside the signed package.

Any unchecked item is `NO-GO`. Successful launch is not scientific readiness.
