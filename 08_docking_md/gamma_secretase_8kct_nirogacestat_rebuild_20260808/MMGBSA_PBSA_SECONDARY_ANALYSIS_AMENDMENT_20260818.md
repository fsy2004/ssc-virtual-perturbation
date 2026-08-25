# Versioned MM-GBSA/MM-PBSA secondary-analysis amendment

Date frozen: 2026-08-18, before inspection of any endpoint-energy result.

## Reason for amendment

The rebuilt native-8KCT/O6U system now has a hash-bound, nonzero-charge ligand
model, a reproducible topology and production release, and three prospectively
seeded 500 ns trajectories. The earlier prohibition on endpoint-energy
calculations is retained for the primary pose-compatibility analysis but is
superseded for the separate secondary analysis defined here.

This amendment does not change any production coordinate, topology, TPR,
random seed, thermostat, barostat, timestep, checkpoint, trajectory endpoint,
or primary structural/membrane gate. It cannot authorize a rerun based on an
energy result.

## Scientific role and claim ceiling

MM-PBSA and MM-GBSA are secondary, exploratory endpoint-energy summaries of
one protein-ligand system under alternative implicit-solvent models. They are
not experimental binding affinities, absolute binding free energies, potency
estimates, target-engagement evidence, or efficacy evidence. There is no
between-ligand ranking and no apo/unbound simulation arm.

The reported quantity is an endpoint-energy score, with and without the
nonpolar term, under the stated model. Entropic terms are not added to the
primary score. Interaction entropy, C2 entropy, quasi-harmonic entropy, and
normal-mode entropy are disabled. No frame is treated as an independent
biological or inferential replicate.

## Eligibility and anti-selection rules

- All calculations use the frozen release archive SHA256
  `5a421f28afee664b5a8919db5f415f1205f35200950117bb3a67fceaba544a98`
  and production-TPR release SHA256
  `a6e41f920f5af4860b7452c4cbdb2afeed8243bf65fb23b4fd6730e3ebbca4aa`.
- A realization must pass endpoint integrity, PBC-distance invariance,
  finite-energy, membrane, and trajectory-continuity checks before its energy
  score is interpretable.
- The formal analysis includes `rep01`, `rep02`, and `rep03`. No realization,
  frame, or time window may be selected because of its energy or ligand
  behavior.
- If any realization is technically ineligible, the all-three endpoint-energy
  summary is No-Go. A diagnostic calculation may be retained with the failure
  label but cannot replace the missing realization.
- No interim single-replica endpoint-energy calculation is permitted. Batch
  execution begins only after all three 500 ns realizations pass their required
  integrity and scientific QC gates.

## Receptor, ligand, and trajectory definition

- Ligand: the neutral 76-atom O6U residue used in the frozen production TPR.
- Receptor: all gamma-secretase protein chains, covalently linked NAG/BMA
  glycans, three retained CLR molecules, and the two retained PC1/DSPC
  structural components.
- Bulk POPC, water, sodium, and chloride are removed only in the derived
  endpoint-energy input; raw trajectories remain immutable.
- Use the validated whole/cluster/nojump/center/rebox pipeline. The membrane
  midplane is placed at `z = 0` without rotating the membrane normal away from
  z. The protein-O6U minimum-image distance must agree with the matching raw
  frame within 0.01 nm.
- Use a single-trajectory protocol: complex, receptor, and ligand coordinates
  are derived from the same frame.

## Frozen sampling

- Window: 200.0 to 500.0 ns for every realization.
- Divide the window into 300 half-open 1 ns strata. Select the frame nearest
  each stratum midpoint: 200.5, 201.5, ..., 499.5 ns. Ties select the earlier
  frame. Thus each eligible realization contributes exactly 300 frames.
- Preserve five fixed 60 ns blocks: [200,260), [260,320), [320,380),
  [380,440), and [440,500] ns.

## Frozen calculation models

The target implementation is `gmx_MMPBSA` 1.6.5 with a hash-pinned isolated
environment and compatible AmberTools. The exact executable versions,
environment export, input files, generated AMBER topologies, commands, and
logs must be retained. A three-frame technical canary must pass before study
frames are processed. A compatibility failure is a technical No-Go and cannot
be repaired by changing scientific parameters after viewing scores.

### Primary endpoint model: membrane-aware MM-PBSA

- `memopt = 1`
- `emem = 7.0`
- `indi = 4.0`
- `exdi = 80.0`
- `istrng = 0.150`
- `poretype = 1`
- `mctrdz = 0.0` after the validated membrane-centering operation
- `mthick` is the across-realization median of the median upper-to-lower
  leaflet phosphate-plane separation over the fixed window, converted to
  Angstrom. The value is computed without reference to endpoint energies and
  is then frozen for all realizations.

Prespecified PB sensitivity: repeat with `indi = 1.0`; all other fields remain
unchanged. A homogeneous-water PB calculation is not a substitute for the
membrane-aware primary model.

### MM-GBSA sensitivity models

GB is reported as a solvent-model sensitivity because it does not represent
the explicit bilayer as the primary PB slab does.

- GB-OBC2: `igb = 5`, `intdiel = 4.0`, `extdiel = 80.0`,
  `saltcon = 0.150`.
- GB-Neck2: `igb = 8`, with the same dielectric and salt settings.
- The nonpolar surface term uses the tool's documented version-pinned default;
  the resolved numerical value is captured in the generated input and cannot
  vary between realizations.

## Validation and reporting

- Verify residue/atom counts, total charge, nonzero O6U charge vector, CHARMM
  topology conversion, receptor/ligand masks, and single-point component
  consistency before batch execution.
- Report each energetic component, total endpoint score, and failure flags for
  every realization and fixed block.
- Report per-realization means, the three-realization mean and SD, and the
  complete five-block values. No frame-level p value is permitted.
- A hierarchical bootstrap over realizations and fixed blocks may be reported
  only as an exploratory precision interval; its resampling seed is 20260818
  and the limitations of three realizations must be explicit.
- Residue decomposition is restricted to the hash-pinned native PLIP contact
  set and is descriptive. No data-driven top-residue list or residue-level
  significance test is permitted.
- Agreement in sign or magnitude across PB and GB models is a robustness
  observation, not validation against experiment. Discordance is retained and
  reported.

## Resource and immutability policy

Endpoint-energy work does not run concurrently with production. It begins only
after all three trajectories and their required gates are complete. Derived
files use a new versioned directory. Raw XTC, EDR, CPT, GRO, log, TPR,
topology, and parameter files are never overwritten.
