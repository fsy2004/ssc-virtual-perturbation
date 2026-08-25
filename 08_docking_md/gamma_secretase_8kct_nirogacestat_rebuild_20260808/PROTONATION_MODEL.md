# Frozen protonation model

## Single production state

The following fixed state is used in the single production model at the modelled pH of 7.4:

- nirogacestat: neutral O6U, total formal charge 0;
- PSEN1 Asp257: deprotonated carboxylate, total side-chain charge -1;
- PSEN1 Asp385: protonated neutral carboxylic acid, total side-chain charge 0, with the proton assigned by the accepted CHARMM residue patch and recorded atom mapping.

No alternative ligand microstate or catalytic-dyad state is simulated. These choices are frozen before membrane construction and may not be changed after trajectories are inspected. They define one controlled approximation, not a sensitivity analysis.

## O6U identity gate

The official PDB chemical-component record defines O6U as the neutral `(2S,2S)` stereoisomer with InChIKey `VFCRKLWBYMDAED-REWPJTCUSA-N`. The hydrogen-complete component has 76 atoms. Its CCD bonds include N04--H39, N06--H41, and N07--H43; N05 and N08 have no attached hydrogen. The 8KCT deposited heavy-atom coordinates are retained, only CCD-consistent hydrogens are added, and hydrogen relaxation is performed with heavy atoms restrained.

The prescribing information reports macroscopic pKa values of 5.77 and 7.13. Under a simple sequential free-solution model at pH 7.4, the neutral macrostate is the largest population. The binding pocket can alter this equilibrium, so the neutral state is described as a predeclared modelling choice anchored to the deposited chemical component, not an experimentally measured bound-state protonation.

## Catalytic dyad rationale

Cryo-EM does not resolve the relevant hydrogens, so this state is a model choice. Structure-specific PDB2PQR/PROPKA preparation of the 8KCT protein at pH 7.4 predicted pKa values of 5.97 for Asp257 and 7.77 for Asp385, assigning Asp257 as deprotonated and Asp385 as protonated. Independent gamma-secretase simulations provide convergent, although context-dependent, support: multiscale work estimated Asp257 and Asp385 pKa values of approximately 5.12 and 9.91, respectively, and found that protonated Asp385 favoured the active-state ensemble; pH-replica-exchange work found physiological-pH apo dyad deprotonation but shifted the substrate-bound dyad toward monoprotonation and assigned Asp257 as general base and Asp385 as general acid. Recent QM/MM work found that protonated Asp257 can be favoured for the APP chemical step, underscoring that the preferred site depends on biochemical state and question.

The present system is inhibitor-bound rather than an apo enzyme or a reacting substrate complex. We therefore fix the structure-specific Asp257(-1)/Asp385(0) assignment before membrane construction. No post hoc state change or second protonation arm is allowed.

This evidence is model-dependent. The manuscript must therefore call this a `predeclared fixed-protonation model`; it must not call the dyad state experimentally established or universally physiological. The MD claim is limited to pose compatibility conditional on this fixed state.

## Sources

- O6U record: https://www.rcsb.org/ligand/O6U
- 8KCT structure: https://www.rcsb.org/structure/8KCT
- Nirogacestat label: https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=f172e6ff-3190-41b0-b95a-58a7ef9e9e1e
- Physiological-pH dyad modelling: https://doi.org/10.1021/acschemneuro.2c00563
- Multiscale dyad-state modelling: https://doi.org/10.1039/C7SC00980A
- Inhibitor-bound gamma-secretase protonation analysis: https://pmc.ncbi.nlm.nih.gov/articles/PMC9282858/
- Substrate-hydrolysis protonation study: https://doi.org/10.1021/acs.jpcb.4c04085
