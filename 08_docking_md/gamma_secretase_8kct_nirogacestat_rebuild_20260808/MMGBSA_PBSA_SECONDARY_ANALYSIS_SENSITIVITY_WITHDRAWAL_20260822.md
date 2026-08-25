# Versioned endpoint-energy sensitivity-withdrawal amendment

Date frozen: 2026-08-22, before any formal endpoint-energy result was generated.

This amendment supersedes only the sensitivity-model portion of
`MMGBSA_PBSA_SECONDARY_ANALYSIS_AMENDMENT_20260818.md`. It does not change any
production coordinate, topology, TPR, trajectory endpoint, primary
structure/membrane/energy gate, sampling window, frame-selection rule, residue
decomposition set, entropy policy, or claim ceiling.

## Decision

The secondary exploratory endpoint-energy analysis now contains one executable
model only:

- membrane-aware MM-PBSA: `memopt = 1`, `emem = 7.0`, `indi = 4.0`,
  `exdi = 80.0`, `istrng = 0.150`, `poretype = 1`, `mctrdz = 0.0`, with
  `mthick` frozen from the fixed-window membrane geometry before endpoint
  energies are inspected.

The following sensitivity branches are withdrawn and must not be run or used in
the manuscript evidence chain:

- PB `indi = 1.0`
- GB-OBC2 (`igb = 5`)
- GB-Neck2 (`igb = 8`)

## Rationale and reporting boundary

The withdrawn branches were planned only as model-sensitivity diagnostics. They
do not strengthen the biological or clinical interpretation of this one-system
application study, increase compute and audit burden, and could distract from
the already narrow endpoint-energy claim ceiling. Removing them before formal
endpoint-energy execution reduces nonessential multiplicity without selecting
on results.

The remaining endpoint-energy output is still secondary and exploratory. It is
not an experimental affinity, an absolute binding free energy, potency,
target-engagement, efficacy, target-occupancy, or causal evidence. It must be
reported as one membrane-aware PB endpoint-energy score under a specified
implicit-solvent approximation, with the same three-realization/five-block
descriptive and bootstrap limitations defined in the 2026-08-18 amendment.

## Audit rule

Any stale canary, resource plan, formal-run script, summary, automation prompt,
or report that still expects `PB_membrane_indi1`, `GB_OBC2`, or `GB_Neck2` is
out of date for formal execution and must fail closed or be regenerated under
the v2 plan before secondary endpoint-energy work starts.
