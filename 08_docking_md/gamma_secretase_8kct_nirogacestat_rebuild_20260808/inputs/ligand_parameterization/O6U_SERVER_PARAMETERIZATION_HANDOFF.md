# O6U server parameterization handoff

## Release state

Local identity validation has passed, but MD parameterization remains blocked. The server must start from the audited explicitly 3D derivative of the official hydrogen-complete neutral O6U SDF and the 76-row correspondence table in this directory. The derivative changes only the molfile 3D marker and provenance fields; its atom order, coordinates, graph, charge, and stereochemistry were revalidated against the official source. No old `LIG`, +2-charge, docking-PDBQT, MOL2 charge, RTF, PRM, ITP, or topology file may be reused.

## Immutable inputs

| Input | SHA-256 |
|---|---|
| `../reference/O6U_ideal.sdf` (immutable official source) | `985ac8899eb97efff73a57682da4794fbe55acf06be25717f9d885d64cb5962a` |
| `O6U_neutral_hydrogen_complete_3D.sdf` (audited server input) | `2cb9d769cde4157181a6199b83294cad56cade14ab34a5e86a6deb6790fc28d5` |
| `O6U_preparameterization_atom_correspondence.tsv` | `a773708cf030ba03cdee924a51359b84810742fbce63b99cc95b8d6309e2eca3` |
| `../reference/O6U.cif` | `c5c3c1ac73c9cb512612a855e39caddd9d840eaf823822f677c2675256659acd` |
| `../reference/8KCT.cif` | `2ed75442ca2c503a014b4e5e8bac67e107201c31776238d0ae94069b39013da9` |

Expected identity: `C27H41F2N5O`, formal charge 0, 76 total atoms, 35 heavy atoms, 41 explicit hydrogens, 78 bonds, one connected component, no radicals, and RDKit-assigned stereocentres `(10,S)` and `(13,S)` in the official SDF zero-based ordering. The native MD heavy-atom coordinates are the single deposited O6U at 8KCT chain B residue 502.

## Server stages

### 1. Freeze the toolchain before generating parameters

- Record operating system/container digest, scheduler, CPU/GPU model, CHARMM executable/version, GROMACS exact version, QM engine/version, FFParam version/commit, and CGenFF program/rule/topology/parameter releases.
- Obtain the initial O6U assignment through the official CHARMM/CGenFF-compatible route (for example, CHARMM-GUI Ligand Reader or the licensed current CGenFF/ParamChem service).
- Keep the exact uploaded SDF, unedited raw response/download, commands, standard output/error, and SHA-256 inventory.
- Do not combine atom typing or charges from one CGenFF release with topology/parameter files from another release.

### 2. Fail immediately on identity or initial-assignment errors

- Map all 76 CCD atoms one-to-one into the CGenFF and final GROMACS names; no missing, duplicate, reordered without an explicit mapping, or extra atom is allowed.
- Confirm total topology charge `0.0000 +/- 0.0001 e`, a nonzero partial-charge vector, unchanged connectivity, and the two S stereocentres.
- Export every atom charge penalty and every bond, angle, Urey-Bradley, improper, and dihedral penalty. Inventory all values above 10; values above 50 require extensive validation/reoptimization.
- A generated topology is only an initial hypothesis. Stop here if the raw output, version, penalty table, or atom mapping is incomplete.

### 3. Generate retained QM target data with the available official FFParam release plus independent engines

- Use official FFParam 1.2.0 only for the additive CHARMM operations it implements and the exact model chemistry documented in `LIGAND_PARAMETERIZATION_PROTOCOL.md`. Generate missing Hessian/PED and headless automation artifacts with separately versioned Psi4/OpenMM scripts. Never report the combined workflow as FFParam-v2 unless an official v2 package is later obtained and retained.
- Retain the complete validated 771-conformer neutral CREST ensemble with its continuous energy annotations. Use the frozen deterministic TFD 0.20 diversity selection for expensive MP2 calculations: the global minimum (frame 1) plus centroid frames 342, 641, 679, and 768. Keep the native experimental conformer as a separate mapping/reference structure. Do not introduce a post hoc energy-window exclusion or expand expensive QM to every conformer. After MP2 optimization, collapse two representatives only when both heavy-atom RMSD is below 0.1 Angstrom and absolute energy difference is below 0.1 kcal/mol. Retain full frequency/Hessian/PED evidence for every representative used as a fitting target.
- Calculate applicable polar-site water interactions and the molecular dipole using the published CHARMM neutral-molecule target convention. Add weakly weighted MP2 electrostatic-potential data to constrain sites that water probes cannot determine. Freeze the published expanded-training relative weights before fitting: water interaction 10.0 kcal^-1 mol, dipole 3.0 Debye^-1, and ESP 1.0 kcal^-1 mol e^-1. Apply the published initial-charge restraint atom weight `1/(0.1 + charge_penalty)`. This is an auxiliary CHARMM charge-fitting target, not a RESP/GAFF charge model.
- Classify every high-penalty and pocket-relevant candidate torsion before target generation. Use full-range adiabatic scans only for genuine free rotors, stiff local PES targets for double/aromatic/conjugated/amide-like bonds, and low-energy conformer/puckering plus local PES targets for ring-internal terms. Treat multiple Fourier components for one quartet as one coupled target.
- Use the whole O6U molecule unless a chemically justified fragment is predeclared; if a fragment is necessary, retain its atom mapping and verify the final whole-molecule conformer ordering, electrostatics, and geometry.
- Keep input/output files for every converged and failed QM job. Failed jobs may be corrected prospectively and rerun, but may not be silently omitted.

### 4. Fit only the unsupported terms and cross-validate the whole molecule

- Preserve initial and final RTF/STR/PRM files and a machine-readable change table for charges and each changed bonded term. Fit the high-penalty charge region as a constrained molecular vector; neighbouring lower-penalty atoms may change only when supported by the joint water/dipole/weak-ESP targets, penalty-weighted restraint toward the initial charges, exact total charge, and chemically justified equivalence constraints.
- Recheck water interaction energies/distances, dipole magnitude/direction, weak-ESP agreement, optimized geometry, Hessian/normal-mode balance, torsional minima/barriers/profile shape, and conformer ordering after every fit cycle.
- Do not optimize a single attractive contact or one docked/native pose. The fit targets are transferable molecular properties, not the 8KCT protein environment.
- Do not use RESP/GAFF-style charges inside the CHARMM additive model.

### 5. Cross-engine energy regression before any membrane job

- Generate CHARMM and GROMACS representations from the same final parameter set and the same coordinates.
- Compare bonded, electrostatic, van der Waals, and total single-point energy components on multiple retained conformers. Record unit conversions and exclusions explicitly.
- Require the prospective tolerance in `LIGAND_PARAMETERIZATION_PROTOCOL.md`; any unexplained component mismatch blocks PDB Reader/Quick Bilayer.

### 6. Cheap pre-membrane canaries

- Run vacuum/isolated-ligand minimization and a short explicit-water stability rehearsal only as gross implementation checks. Inspect energy, geometry, chirality, constraints, and every warning; these runs do not validate binding or replace QM targets.
- Re-run `audit_o6u_preparameterization.py` against the immutable source and bind the accepted final 76-row CGenFF/GROMACS mapping to the signed ligand parameter record.
- Two reviewers must sign: one for QM/parameter fitting and one for topology/atom mapping.

## Required return bundle

1. Environment and licensed-tool version record.
2. Raw initial CGenFF/Ligand Reader inputs, outputs, and full penalty tables.
3. Completed 76-row CCD/SDF/CGenFF/GROMACS atom correspondence.
4. All QM inputs/outputs plus convergence manifest.
5. FFParam configuration, commands, target tables, plots, initial/final parameter files, and machine-readable fit metrics.
6. CHARMM/GROMACS energy-component regression inputs and results.
7. Cheap-canary commands, logs, trajectories, and audit report.
8. Completed `ligand_parameter_record` with hashes and two signatures.

## Hard stop

Do not create or submit the PDB Reader/Quick Bilayer production model until `O6U_PREPARAMETERIZATION_AUDIT.json` still matches the immutable inputs and the final ligand parameter record passes every gate. A high-penalty term without retained target data, an unexplained energy mismatch, or an unreviewed atom-name conversion is a NO-GO, not a reason to continue and hope that the 5-ns membrane canary will detect the problem.
