# Local pre-server audit

## Frozen study design

- One biochemical system: the native 2.60 Angstrom cryo-EM 8KCT gamma-secretase-O6U complex.
- No apo or active-compound comparator is included.
- One CHARMM-GUI membrane construction is planned.
- Three independent dynamic branches start at first NVT velocity generation after one common deterministic minimization.
- Each branch is 500 ns; aggregate production is 1,500 ns.
- The fixed primary analysis window is 200-500 ns in each branch.
- The claim is limited to short-timescale compatibility of the experimentally observed pose under the modeled membrane conditions.

The MD does not use a docking-derived starting pose. Native-pose redocking is retained only as a method QA exercise, and its failure to rank the native-like pose first in all seeds excludes docking from any binding-mode validation claim.

## Local checks completed

- Reference 8KCT PDB/mmCIF, O6U CCD/SDF, and the wwPDB validation report are present with checksums.
- The local O6U preparameterization audit passed exact identity, 76-atom/35-heavy-atom/41-hydrogen composition, neutral charge, 78-bond connectivity, two S stereocentres, and one-to-one correspondence with the single deposited 8KCT O6U heavy-atom pose. It remains explicitly blocked for MD parameter use.
- The deposited 8KCT precuration audit passed the resolved protein segments, component counts, four disulfides, 12 protein-to-NAG links, and nine glycan-internal links. It remains explicitly blocked pending manual PDB Reader review and two-person sign-off.
- Five complete, checksum-pinned full-text articles covering CGenFF, penalty interpretation, CHARMM-GUI Ligand Reader, FFParam-v2, and the expanded CGenFF training set are retained locally with an evidence-to-action matrix.
- The fixed neutral O6U and Asp257(-1)/Asp385(0) model is documented as a predeclared model choice.
- The CHARMM-GUI API token exists outside the project in a user-restricted local secret file. No account credential or bearer token was found in the project tree.
- All package Python and JSON parse checks pass; the new O6U and structure audits also pass their valid-input plus negative-tamper self-tests.
- Core release-gate, staged runner, raw-output, tree-manifest, sealed QC-report, PBC-invariance, and primary postprocessing synthetic/adversarial self-tests passed.
- Native-structure PLIP 3.0.1 was run twice in an isolated local Open Babel environment. The raw text reports were byte-identical, and all 11 canonically sorted O6U interaction records were identical across runs; all displayed endpoint coordinates matched atoms in the prepared structure exactly.
- The native 8KCT-O6U PyMOL contact panel passed deterministic-script reconstruction, input/output hash, log, camera, PSE and 600-dpi PNG validation. The validated native-structure panel is included in the English manuscript as Figure 7; no docking-derived pose or score is shown.
- The strict validator rejects the draft template, as intended.
- The build and run directories are empty; no production manifest or analysis plan exists. Accidental production launch is therefore blocked.

## Hard blockers before CHARMM-GUI construction

1. **Ligand parameters are not accepted.** Neutral O6U still needs the frozen CGenFF release, full penalty inventory, targeted FFParam/QM validation, a CHARMM-to-GROMACS energy regression, and independent sign-off. The old +2, high-penalty ligand model is retired and must not be reused.
2. **PDB Reader curation is not approved.** The segment breaks, termini, four disulfides, 12 protein-to-NAG links, nine glycan-internal links, 18 NAG, 3 BMA, 3 resolved CLR, 2 resolved PC1/DSPC, native O6U, dyad protonation, and PSEN1 NTF/CTF separation require a signed structure record and visual audit.
3. **The membrane is not built.** No PDB Reader job ID, Quick Bilayer job ID, PPM orientation record, CHARMM-GUI archive, GROMACS topology, coordinate file, or build-validation report exists.
4. **The real execution environment is not frozen.** Exact GROMACS/container/CUDA/GPU/driver versions and measured storage throughput remain unknown until the server is available.
5. **Build-specific analysis definitions do not yet exist.** Atom indices, native contacts, membrane groups, APL/gorder mappings, and numeric physical acceptance records can only be finalized against the accepted build.

## Cost-controlled execution sequence

### Stage 1: chemistry and structure only

Use a short CPU-capable server allocation for FFParam/QM work and Linux PLIP. Do not rent the long GPU allocation at this stage. Stop if any high-penalty term remains unresolved, stereochemistry/charge changes, the native pose is altered, or the CHARMM/GROMACS energy regression fails.

### Stage 2: one CHARMM-GUI build

Approve the PDB Reader record manually, then submit exactly one Quick Bilayer job through the pinned API client. Validate and visually inspect the downloaded archive locally before any MD server is started. The planned bulk membrane is symmetric POPC with 0.15 M NaCl while retaining the experimentally resolved 3 CLR and 2 PC1/DSPC lipids.

### Stage 3: environment, minimization, and equilibration

Rent the GPU server only after chemistry and build validation pass. Freeze the exact runtime, run the common minimization, and create all three independently seeded full equilibration branches. Any LINCS, topology, component-count, membrane, density, temperature, pressure, PBC, or checkpoint defect is a hard stop.

### Stage 4: mandatory 3 x 5 ns canary

Generate each 500-ns production TPR exactly once, run 5 ns from each TPR, and retain those frames as the beginning of production. Continue only when all three technical canary reports pass and the measured storage projection leaves at least 30% headroom. The runner resumes the same TPR/checkpoint; it may not regenerate a TPR, discard a realization, or choose a seed based on attractive ligand behavior.

If all three canaries show immediate gross native-pose loss, the study should stop for scientific adjudication rather than spending the remaining budget. A single unfavorable realization cannot be dropped or replaced. A later loss of pose is a scientific result, not a technical excuse to rerun.

### Stage 5: checkpoint continuation

Only after the all-three canary release gate passes may the same three trajectories continue to 500 ns. Progress is checkpointed and hash-chained; every restart must append without overlapping time.

## What can and cannot be guaranteed before 1,500 ns

The local and canary gates can eliminate identifiable chemistry, topology, build, runtime, PBC, output, storage, and restart failures before full expenditure. No scientifically honest workflow can guarantee in advance that an unbiased 500-ns trajectory will preserve a ligand pose or show the desired result. The study is protected from wasted cost by staged technical release and from biased inference by prohibiting replica replacement, window movement, smoothing, spike deletion, and post hoc model changes.

## Analysis boundary

Primary outputs are full 0-500 ns QC traces and prespecified 200-500 ns replica-level summaries of protein, pocket, ligand, native-contact, hydrogen-bond, membrane, and water-defect metrics. MM/GBSA, MM/PBSA, per-residue decomposition, PCA, FEL, and 3D free-energy/population surfaces are prohibited. Frames are not independent replicates; the independent stochastic unit is the velocity-seeded trajectory (`n = 3`).

## Current decision

**NO-GO for membrane submission and GPU rental.** The code and local failure gates are ready, but ligand-parameter acceptance and signed PDB Reader curation must be completed first.
