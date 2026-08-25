# 8KCT structure-curation protocol

Status: fail-closed manual checkpoint before CHARMM-GUI submission.

## Source identity

- Use RCSB entry 8KCT, biological assembly 1, version recorded in the structure manifest.
- Preserve the deposited 2.60 Angstrom cryo-EM heavy-atom coordinates and the author-assigned four-subunit human gamma-secretase assembly.
- Verify the local `8KCT.cif`, `8KCT.pdb`, wwPDB validation report, and O6U CCD files against `inputs/reference/SHA256SUMS.txt`.
- Do not substitute an AlphaFold model, another gamma-secretase entry, a docked complex, or a coordinate file modified by a visualization program.

## Resolved protein model

The accepted coordinate record must reproduce these resolved regions before any topology is generated:

| Subunit | Author chain | Resolved residues | Required treatment |
|---|---|---|---|
| Nicastrin | A | 34-700 | retain resolved coordinates; do not build unresolved terminal residues automatically |
| PSEN1 | B | 76-291 and 377-467 | treat the mature NTF and CTF as separate topology segments; never create a peptide bond across 291-377 and never fill residues 292-376 |
| APH-1A | C | 2-244 | retain resolved coordinates; do not build unresolved terminal residues automatically |
| PEN-2 | D | 6-101 | retain resolved coordinates; do not build unresolved terminal residues automatically |

There are no accepted alternate-location choices or partial-occupancy conformers in the deposited model. If a later RCSB version changes that fact, the current manifest is invalidated and curation restarts.

Artificial coordinate truncations must not create spurious charged termini in the membrane or pocket. The curator must map every terminus to UniProt and distinguish a biologically authentic free terminus from an unresolved internal continuation. Exact CHARMM terminal patches are recorded and reviewed in PDB Reader; artificial truncations are neutralized, while authentic mature termini retain the chemically appropriate state. No patch is inherited silently from a default.

## Covalent and glycan topology

Preserve and verify the four nicastrin disulfides:

- Cys50-Cys62;
- Cys140-Cys159;
- Cys230-Cys248;
- Cys586-Cys620.

Retain the deposited 18 NAG and 3 BMA residues at the 12 nicastrin N-glycosylation sites: Asn45, Asn55, Asn187, Asn264, Asn387, Asn435, Asn464, Asn506, Asn530, Asn562, Asn573, and Asn580. Build only the glycan atoms resolved in 8KCT. Do not complete a consensus glycan tree beyond the density and do not strip glycans to reproduce older lower-resolution simulations.

## Native small molecules and structural lipids

Retain exactly:

- one O6U nirogacestat molecule in its deposited heavy-atom pose;
- three CLR cholesterol molecules;
- two PC1 molecules, chemically 18:0/18:0 phosphatidylcholine (DSPC), not POPC.

O6U hydrogens are added according to the neutral CCD microstate and relaxed with deposited heavy atoms restrained. The ligand must map one-to-one to the signed CGenFF/FFParam parameter record.

The resolved CLR and PC1 molecules are structural seed lipids. They are not converted into bulk POPC and are counted separately when the final leaflets are audited. CHARMM-GUI must pack the bulk bilayer around them without duplication, deletion, ring penetration, or steric overlap.

8KCT contains no resolved water molecules. Do not hand-place favorable pocket waters. Water enters through the recorded solvation procedure and is allowed to equilibrate dynamically.

## Fixed protonation and orientation handoff

- Apply the single state in `PROTONATION_MODEL.md`: neutral O6U, deprotonated PSEN1 Asp257, and protonated PSEN1 Asp385.
- Record all other titratable-residue assignments produced by PDB Reader and manually review buried, interfacial, and pocket residues. A software default is not evidence.
- Submit un-oriented curated PDB Reader coordinates to the single Quick Bilayer job and request `ppm=true`. PPM is applied exactly once during that build. Afterward, save and inspect the transformation matrix, oriented coordinates, membrane boundaries, job ID, screenshots, and hashes in the build-level orientation record. Do not pre-transform the input or reorient individual realizations.

## Mandatory PDB Reader evidence

The signed structure record must contain:

- RCSB input hashes and PDB Reader job ID;
- exact sequence/residue and segment mapping;
- every terminal and internal patch;
- all four disulfides, all 12 protein-to-NAG links, and all nine glycan-internal links;
- O6U, CLR, PC1, NAG, and BMA counts before and after PDB Reader;
- complete ligand atom-name mapping and total topology charge;
- all protonation decisions;
- visual screenshots of the pocket, termini, glycans, and structural lipids before membrane construction;
- curator and independent reviewer sign-off.

## Hard stops

Production is blocked if PDB Reader or any later step fills PSEN1 residues 292-376, bonds the NTF to the CTF, strips or duplicates a native component, changes O6U stereochemistry/charge/pose, loses a disulfide/glycan link, adds an undocumented terminal charge, or leaves a serious wwPDB clash/geometry issue unresolved.

The subsequent build is also blocked if PPM is applied zero or more than one time, or if its build-level transformation/orientation record is absent.
