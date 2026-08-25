# Nirogacestat parameterization protocol

## Frozen chemical model

Only neutral O6U is parameterized. Its identity is fixed by the RCSB O6U component definition, including both stereocentres, bond orders, atom names, 76 atoms in the hydrogen-complete component, and formal charge zero. The native 8KCT heavy-atom coordinates are preserved. No alternate protomer, tautomer, salt form, or docking-derived geometry enters production.

## Initial assignment

1. Run the frozen CGenFF release on a hydrogen-complete O6U structure whose heavy atoms map one-to-one to the RCSB component.
2. Retain the exact input, program and database versions, full output, atom mapping, topology, parameter file, charge sum, and every per-atom/per-term penalty.
3. Reject an all-zero partial-charge vector, a non-integer total charge, missing atom/bond, changed stereocentre, or any topology charge differing from zero by more than 0.001 e.
4. Inventory every charge or bonded term with penalty above 10. A penalty in the 10--50 range requires targeted validation; a value above 50 requires extensive optimization. A successful topology build is not validation.

## CHARMM-compatible QM targets

Use the officially distributed FFParam 1.2.0 package for the CHARMM additive target-generation and comparison operations that it implements. Use separately retained Psi4/OpenMM calculations for the Hessian, normal-mode/PED, command-line automation, and other validation operations that are not implemented in that release. Record the exact tool responsible for every artifact and preserve all raw inputs and outputs; do not label the workflow as FFParam-v2 unless an official v2 distribution is actually obtained and hash-pinned. At minimum:

- optimize the relevant neutral conformer ensemble and calculate molecular dipoles at MP2/6-31G(d), or use a documented higher-level method only if the mapping back to the CHARMM target convention is justified. After QM optimization, treat two conformers as duplicates only when heavy-atom RMSD is below 0.1 Angstrom and the absolute energy difference is below 0.1 kcal/mol, following the expanded CGenFF training workflow;
- generate donor/acceptor water-interaction targets at HF/6-31G(d), apply the standard neutral-compound CHARMM scaling and distance correction, and compare every chemically applicable polar interaction with the MM result;
- fit the large, coupled heteroaromatic-amide charge region against the standard CHARMM additive water-interaction and molecular-dipole targets, augmented by weakly weighted MP2 electrostatic-potential data where water probes do not constrain buried carbon sites. Use the published expanded-training relative target weights: water interaction 10.0 kcal^-1 mol, dipole 3.0 Debye^-1, and ESP 1.0 kcal^-1 mol e^-1. Fit charges as one constrained molecular vector using exact total charge, chemically justified equivalence constraints, and a restraint toward the initial CGenFF vector with atom weight `1/(0.1 + charge_penalty)`; do not optimize the seven high-penalty atoms independently and do not import RESP/GAFF charges;
- calculate the Hessian, normal modes, and potential-energy distribution used to assess high-penalty bonds and angles, including low-frequency mode balance;
- classify each candidate torsion by central-bond chemistry before generating a target. Use a full-range adiabatic PES only for a genuinely free rotor. Use a stiff PES near the accessible minimum for double, aromatic, conjugated, and amide-like bonds. Treat ring-internal terms through low-energy ring conformers/puckering states and local PES data. Multiple Fourier terms on one quartet are one coupled target, not separate scans;
- compare QM and MM minima, barriers, conformer ordering, molecular geometry, water-interaction energies/distances, dipole magnitude/direction, weakly weighted ESP agreement, Hessian/PED balance, and category-correct torsional profiles.

The original CGenFF protocol identifies approximately 0.2 kcal/mol as an ideal water-interaction energy agreement. A 0.5 kcal/mol line shown in published validation plots is used here only as a prespecified project review cap for an applicable individual interaction, not as a universal CGenFF acceptance rule. The signed report must explain every excluded, unfavourable, or weak water orientation and every larger residual, and must demonstrate that improvements to high-penalty terms do not degrade the remaining target set.

## Fail-closed acceptance

Production remains blocked unless all of the following are true. The numerical limits below are prespecified project acceptance gates applied to the retained QM/MM target set; they are not universal claims about every CGenFF molecule:

- the atom mapping and stereochemistry are independently reviewed;
- topology and listed partial charges sum to 0.0000 +/- 0.0001 e, and the vector is not all zero;
- every initial penalty above 10 has a raw target-data artifact, fit comparison, and signed disposition;
- no unresolved high-penalty term remains;
- all QM jobs used as targets converged and their outputs are retained, not summarized only in images;
- initial and final parameters, fitting configuration, software versions, commands, numerical metrics, and SHA-256 hashes are recorded;
- targeted bond/angle fitting aims for at most 0.02 Angstrom and 2 degrees, with hard review limits of 0.03 Angstrom and 3 degrees for the targeted internal-coordinate set;
- targeted vibrational-frequency mean absolute relative deviation is at most 5 percent and the signed PED review confirms acceptable low-frequency mode ordering and internal-coordinate balance; the scalar average alone cannot release the model;
- every retained applicable water-interaction energy differs by at most the project review cap of 0.5 kcal/mol, the mean absolute difference is at most the ideal 0.2 kcal/mol target, and corrected target distances differ by at most 0.2 Angstrom. Exclusions must be chemically justified and retained before fitting decisions are made;
- the MM/QM dipole-magnitude ratio lies from 1.20 to 1.50 under the published CHARMM charge-fitting convention. Dipole direction is reported and reviewed together with the water targets, but no unsupported universal angular cutoff is asserted;
- torsion validation uses prespecified, category-specific angular domains, energy-zero alignment, periodic matching, and low-energy weighting. The report must reproduce accessible minima, relevant barriers, conformer ordering, and profile shape without degrading other targets; published case-study RMSE, angle, or barrier values are not treated as universal CGenFF cutoffs;
- CHARMM and GROMACS single-point energy components agree within the project implementation-regression tolerance of 0.1 kcal/mol or 0.1 percent, whichever is larger, on identical coordinates and explicitly identical exclusions, 1-4 rules, nonbonded settings, units, and numerical precision. A failure is corrected in the conversion/settings, not by changing chemical parameters to force agreement;
- a second reviewer signs the final parameter record before the single CHARMM-GUI job is submitted.

Primary references: https://doi.org/10.1002/jcc.21367, https://doi.org/10.1021/acs.jpcb.4c01314, and https://doi.org/10.1021/acs.jctc.5c00046.
