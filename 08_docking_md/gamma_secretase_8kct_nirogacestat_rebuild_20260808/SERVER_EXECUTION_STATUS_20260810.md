# Server execution status, 2026-08-10

> Historical checkpoint. The live authority is the server run state and the
> dated gate records. As of 2026-08-11, nirogacestat is the unique frozen
> HES1-Notch-aligned candidate; official FFParam 1.2.0 is the authenticated,
> hash-pinned downloadable release. The earlier request below to obtain an
> FFParam-v2 archive is superseded and must not be used as an execution step.

## Release decision

Current state: **NO-GO for PDB Reader, membrane construction, equilibration, canary production, and 3 x 500 ns production**. The server is being used only for retained environment, chemistry, and runtime qualification until the O6U parameter model passes the prospective ligand gates.

## Completed and retained

- The local package was transferred to the AutoDL data volume and the transfer archive passed SHA-256 verification.
- GROMACS 2025.2 reports CUDA support and detects one compatible NVIDIA GeForce RTX 5090.
- The isolated analysis environment imports MDAnalysis 2.10.0, Gemmi 0.7.3, h5py 3.16.0 against HDF5 2.1.0, NumPy 2.4.6, SciPy 1.17.1, and pandas 3.0.5 without warnings. Exact conda package and explicit environment records are retained.
- The nine core release and postprocessing rehearsals pass on the server. This is an automation/runtime result, not scientific validation.
- The hardened release suite was rerun after synchronization and passed on both Windows and the server analysis environment. It now rejects altered production physics, regenerated or changed equilibration TPRs, broken checkpoint chains, incomplete XTC/EDR/log/CPT endpoints, a final GRO not bound to the latest successful mdrun, stale storage budgets, insufficient live filesystem space, full-trace PBC changes, and unreviewed structural/membrane/energy spikes.
- A disposable 884-water, 5-ps NVT CUDA canary completed with short-range nonbonded work, PME, coordinate update, and constraints assigned to the RTX 5090. No fatal error, LINCS warning, segmentation fault, or token-level `NaN` was present. The canary is not a study system and contributes no scientific evidence.
- The neutral O6U full-molecule HF/6-31G(d) calculation completed as a Psi4 runtime canary. Its result is not FFParam target data and is not accepted as force-field validation.
- Official CHARMM-GUI Ligand Reader output was retained from CGenFF program 4.0 with topology/parameter release 5.0. O6U has 76 mapped atoms, 35 heavy atoms, 41 hydrogens, 78 bonds, formal/topology charge zero, and nonzero partial charges.
- The initial assignment contains seven charge terms above penalty 10 and seven bonded parameter lines above penalty 10, including a maximum parameter penalty of 35.500 and a maximum charge penalty of 32.824. Every term is mapped to immutable CCD atom identities. These terms remain unapproved pending targeted QM fitting and validation.
- The PDB Reader derivative preserves all 10,839 deposited coordinate records and splits the unresolved PSEN1 NTF/CTF break without modeling residues 292-376. The CTF is placed on source-unused author chain Z. A primary audit and an independent Gemmi review both pass the exact derivative hash.

## Runtime-canary development record

All failed harness iterations remain retained because they exposed conditions that must not be rediscovered in the study system:

1. A zero-count molecule line exposed a GROMACS 2025.2 `grompp` failure and was removed from the generated disposable topology.
2. The image-level `OMP_NUM_THREADS=25` setting conflicted with an explicit 24-thread command. The runner now sets and records the intended thread count.
3. GROMACS correctly rejected GPU PME for the non-dynamical steep integrator. Disposable minimization is therefore explicitly CPU-only; the dynamical NVT stage is the CUDA test.
4. Pure rigid water contains no GPU-implemented bonded term. The disposable water canary explicitly keeps bonded work on CPU; project-specific bonded offload remains a later accepted-build canary requirement.
5. A substring search incorrectly read `lincs-warnangle` as containing `NaN`. The retained parser now uses a token-boundary expression and still rejects a true numerical `NaN`.
6. The corrected version passed and is the only accepted generic CUDA canary.

## Bound evidence

- `server_records/gromacs_gpu_canary_v6/GROMACS_GPU_CANARY_REPORT.json`, SHA-256 `ccd8d6b843c7aad7ed3f64503d1114ed50ecaa66c930520a5010461ce5300cc5c`.
- `server_records/toolchain/analysis_environment_v3/SERVER_TOOLCHAIN_RECORD.json`, SHA-256 `11cb5c5ead5621732b1c80eee72aa4bfb088f5c37d838f1f95b5db3d63ecc614`.
- `server_records/toolchain/analysis_environment_v3/gromacs_executable_identity.json`, SHA-256 `eafdcbd3cf28b7b0495d2e84d02d2cc6c389e144ad20264555485ccc31676ee4`.
- `server_records/server_rehearsal/SERVER_REHEARSAL_REPORT.json`, SHA-256 `114b58bbf1ff2800ae68ba7925194831f6eff8bdf613e159343c6ce0a68f41d28`.
- `inputs/structure_precuration/8KCT_pdbreader_input_chain_split.pdb`, SHA-256 `d48f06eea12d000d49ab1d7c539d7d842895dd9b1b0457e218e7897cf0d08610`.
- `inputs/ligand_parameterization/O6U_CCD_CGENFF_ATOM_CORRESPONDENCE.tsv`, SHA-256 `62b5a9500a0c5e0c2d85eb3fa51fb4e4cb82881dafd3d0b2b38f02231d6935f5`.

## Mandatory next release sequence

1. Use the retained official FFParam 1.2.0 release only for operations it implements; use separately versioned Psi4/OpenMM tooling for missing target-generation and validation operations, and do not report the combined workflow as FFParam-v2.
2. Complete and independently validate the frozen five-member MP2/6-31G(d) representative ensemble. Generate the prescribed formal water-interaction inputs from the frozen charge/water geometry, independently reconstruct their geometry audit, and freeze a chemical/visual disposition for every registered direction before any water-interaction QM is authorized.
3. Generate and retain the prescribed whole-molecule QM targets, fit only unsupported terms, cross-validate the complete target set, and independently review the final 76-atom mapping and parameter changes.
4. Pass the CHARMM-versus-GROMACS single-point component regression and the isolated-ligand implementation canaries.
5. Submit the already audited chain-split derivative to one PDB Reader job, review every generated segment, terminus, disulfide, glycan link, heteroatom, missing atom, and ligand mapping, then sign the structure record.
6. Create one Quick Bilayer system only after chemistry and PDB Reader release. Validate orientation, components, topology tree, membrane packing, and project-specific GROMACS preprocessing before minimization.
7. Run one common minimization, three first-NVT-seeded complete equilibration branches, and the mandatory same-TPR 3 x 5 ns project canary. Each equilibration stage keeps one immutable grompp/TPR and resumes only from its own checkpoint. Full continuation remains locked until all three XTC/EDR/log/GRO/CPT endpoints and the measured plus live-storage gates pass.
8. Continue each accepted realization to 500 ns without regenerating its TPR. Analyze only the fixed 200-500 ns interval while retaining complete 0-500 ns traces.

No project production trajectory has been started.
