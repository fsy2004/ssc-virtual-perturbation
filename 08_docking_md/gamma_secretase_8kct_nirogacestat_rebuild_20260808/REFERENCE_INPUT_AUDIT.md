# Frozen reference-input audit

The reference files were downloaded from RCSB on 2026-08-08. Exact SHA-256 values are stored in `inputs/reference/SHA256SUMS.txt`; the acquisition script can repeat the download but must never overwrite an approved release without a new manifest version.

## 8KCT

- Entry 8KCT is a 2.60 Angstrom single-particle cryo-EM model of human gamma-secretase bound to nirogacestat.
- The deposited biological assembly contains nicastrin, APH-1A, PEN-2, and the mature PSEN1 fragments.
- O6U is the deposited nirogacestat component. The production ligand pose is the deposited heavy-atom pose, not a docking prediction.
- The deposited model contains one O6U, three CLR, two PC1/DSPC, 18 NAG, and 3 BMA residues at 12 nicastrin N-glycosylation sites. Its `struct_conn` records specify four disulfides, 12 protein-to-NAG links, and nine glycan-internal links. All are retained and explicitly counted in the structure/build records.
- Nicastrin residues 34-700, APH-1A residues 2-244, PEN-2 residues 6-101, and PSEN1 residues 76-291 plus 377-467 are resolved. The PSEN1 292-376 gap is not modelled and no topology bond may cross it.
- The coordinate model has no alternate locations, partial-occupancy conformers, or resolved water molecules. A future RCSB-version change invalidates that statement and triggers recuration.

## O6U

- Formula: C27 H41 F2 N5 O.
- RCSB formal charge: 0.
- Atom count including hydrogens: 76.
- Chiral centres: 2.
- The exact RCSB atom names, bonds, stereochemistry, and neutral hydrogen assignment define the only ligand microstate in this study.

The prescribing information reports pKa values of 5.77 and 7.13. This means the neutral free base is the majority free-solution macrostate at pH 7.4, but it does not prove a unique bound-state population. The present study deliberately adopts the neutral deposited chemical-component model as a single fixed approximation and does not claim to characterize protonation equilibria.

## Sources

- RCSB entry: https://www.rcsb.org/structure/8KCT
- RCSB component: https://www.rcsb.org/ligand/O6U
- Primary structure report: https://doi.org/10.1038/s41594-024-01439-8
- Prescribing information: https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=f172e6ff-3190-41b0-b95a-58a7ef9e9e1e
