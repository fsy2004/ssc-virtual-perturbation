# Secondary endpoint-energy execution supplement (2026-08-20)

## Status and scope

This supplement is frozen before any MM/PBSA or decomposition result exists.
It operationalizes the fixed scientific scope in
`MMGBSA_PBSA_SECONDARY_ANALYSIS_AMENDMENT_20260818.md`, as superseded for
sensitivity withdrawal by
`MMGBSA_PBSA_SECONDARY_ANALYSIS_SENSITIVITY_WITHDRAWAL_20260822.md`, and
`config/secondary_endpoint_energy_plan_v1.json`. It does not alter production
MD, seeds, sampling, replica inclusion, or the interpretation ceiling.

## Version-pinned implementation defaults

All calculations use gmx_MMPBSA 1.6.5 at source commit
`64e994c71aaff315f3c82dd0852919aecb1ab62e`, Python 3.11.8,
AmberTools 23.3, GROMACS 2023.4, and the captured isolated-environment record.
OpenMPI 4.1.6 and ParmEd 4.3.0 are fixed. CHARMM radii are selected with
`PBRadii=7`.

The single membrane PB model retains the already frozen `memopt=1`, `emem=7`,
`indi=4`, `exdi=80`, `istrng=0.15 M`, `poretype=1`, the measured common
`mthick`, and `mctrdz=0`. The remaining
numerical fields are fixed to the version-pinned CHARMM membrane example:
`radiopt=0`, `fillratio=1.25`, `inp=2`, `sasopt=0`, `solvopt=2`, `ipb=1`,
`bcopt=10`, `nfocus=1`, `linit=1000`, `eneopt=1`, `cutfd=7.0`,
`cutnb=99.0`, `maxarcdot=15000`, and `npbverb=1`.

The previously listed PB `indi=1`, GB-OBC2, and GB-Neck2 sensitivity branches
are withdrawn before formal endpoint-energy execution and must not be run,
summarized, or used in the manuscript evidence chain.

## Fixed descriptive decomposition

Per-residue decomposition uses `idecomp=2`, `dec_verbose=0`, and
`csv_format=1`. The printed set is restricted to the nine unique protein
residues in the approved native PLIP record plus O6U itself:

`B/261,268,272,282,287,380-381,431-432,502`

The source normalized contact record has SHA-256
`69c57604a13cb454ee55bd98b4614eefaac0527fc0a2636351393aefb874969d`.
The ligand entry `B/502` is retained because gmx_MMPBSA requires the printed
decomposition mask to include receptor and ligand. No distance-derived,
trajectory-derived, or result-ranked residue can be added.

## Execution and reporting

The same single model input, decomposition mask, 200–500 ns midpoint frames,
and fixed blocks apply to all replicas. The three-frame canary must produce
finite total and decomposition outputs before the formal release. Raw outputs
are immutable. Decomposition is descriptive; no residue ranking, residue-level
test, affinity claim, or absolute binding-free-energy claim is permitted.
