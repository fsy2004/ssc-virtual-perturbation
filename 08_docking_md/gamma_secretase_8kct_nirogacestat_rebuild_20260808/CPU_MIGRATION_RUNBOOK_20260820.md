# O6U secondary endpoint-energy GPU-to-CPU runbook

## Purpose

This runbook stages all trajectory-dependent work on the current MD node and
moves only a sealed, hash-bound package to a CPU/RAM node. Formal endpoint
energies remain forbidden until rep01–rep03 are each complete at 500 ns, all
integrity/PBC/membrane/energy gates pass, and every production runner exits.

No step in this document modifies raw XTC, EDR, CPT, GRO, log, TPR, topology,
or parameter files. Every output directory must be new.

## Phase A — current node, only after all-three release

1. Run `scripts/status_new_md_run.py` from the local project and confirm the
   current endpoint, no active production runner, fixed archive/manifest/canary
   hashes, production-release hash, healthy disk, and zero OOM events.
2. Use `scripts/seal_secondary_endpoint_all_three_gate.py` to seal the
   all-three eligibility record from each replica's production-completion,
   trajectory-integrity, PBC, membrane, and energy reports plus a runner-status
   snapshot. It must contain `status=pass`,
   `eligible_replicas=[rep01,rep02,rep03]`,
   `all_three_500ns_complete=true`, `all_required_gates_passed=true`, and
   `production_runners_active=false`, with hashes of every completion and
   analysis gate used to reach that decision.
3. Install the isolated preprocessing environment into a new prefix:

   ```bash
   bash scripts/install_endpoint_preprocess_env.sh /root/autodl-tmp/o6u_secondary_endpoint_energy_v1_20260818/preprocess_env
   ```

   This environment contains Python 3.11.8, NumPy 1.26.4, and MDAnalysis
   2.10.0. It is separate from production GROMACS.
4. For each replica, run the validated whole→cluster→nojump→center/rebox→fit
   workflow first. Then call `prepare_secondary_endpoint_energy_inputs.py`
   against the hash-bound centered/reboxed trajectory, its provenance, its PBC
   report, its membrane atom mapping, the production TPR hash, and the common
   all-three eligibility record.
5. The preparation step writes exactly 300 full-system midpoint frames at
   200.5–499.5 ns, a 3-frame canary at 200.5/350.5/499.5 ns, a full-system GRO,
   a complex-only PDB reference, endpoint groups, the frame map, and a SHA-256
   manifest. Protein–O6U minimum-image distances must agree with the validated
   source to at most 0.01 nm.
6. Run `freeze_endpoint_energy_membrane_geometry.py` using all three passing
   membrane-QC tables and all three preparation manifests. This produces one
   common phosphate-plane `mthick` and fixes `mctrdz=0` for every model and
   replica.
7. Run `build_endpoint_energy_cpu_migration_package.py` with the three
   preparation directories, three self-contained topology bundles, the
   relative `topol.top` path in each bundle, the all-three gate, common geometry,
   frozen plan/amendment, execution defaults/supplement, and normalized native
   PLIP contact record. Retain the uncompressed TAR SHA-256 outside the archive.

The current GPU is not used by these Python operations. Keeping this phase on
the existing node avoids transferring the three 31 GB raw production
trajectories; only the approximately 6% midpoint subsets and small metadata are
migrated.

## Phase B — CPU/RAM node

1. Verify the external archive SHA-256, extract once into a new directory, and
   do not edit the extracted package.
2. Install the CPU toolchain into a new isolated prefix:

   ```bash
   bash code/install_gmx_mmpbsa_1_6_5_cpu.sh /data/o6u_endpoint/toolchain/gmx_mmpbsa_1_6_5
   ```

   The installer fixes Python 3.11.8, AmberTools 23.3, GROMACS 2023.4,
   OpenMPI 4.1.6, ParmEd 4.3.0, mpi4py 4.0.1, tqdm 4.67.1, and
   gmx_MMPBSA 1.6.5 source
   commit `64e994c71aaff315f3c82dd0852919aecb1ab62e`. It installs no CUDA stack.
3. Run `capture_gmx_mmpbsa_toolchain.py` with that prefix. The report records
   package metadata, exact executable hashes, safe version output, and the
   one-thread BLAS/OpenMP environment.
4. Run `run_gmx_mmpbsa_canary.py --execute` on rep01's fixed 3-frame trajectory.
   All four models, generated topology counts and charges, finite energy terms,
   the fixed PLIP-only decomposition, and peak RSS must pass. Do not start the
   formal batch if the canary fails. The command explicitly binds both
   `FINAL_DECOMP_MMPBSA.dat` (`-do`) and `FINAL_DECOMP_MMPBSA.csv` (`-deo`);
   the canary rejects any residue mapping outside the frozen 10-entry set.
5. Run `collect_endpoint_cpu_inventory.py`, then
   `plan_secondary_endpoint_resources.py`. The latter reserves two logical
   CPUs, at least 5%/16 GiB RAM, and at least 10%/50 GiB disk. It uses canary
   peak RSS to prevent PB oversubscription.
6. Run `run_secondary_endpoint_energy_cpu.py`. For each model it starts
   rep01–rep03 concurrently in disjoint `taskset` CPU partitions; the four
   models run sequentially. OMP, OpenBLAS, and MKL threads remain one, so MPI
   ranks consume the node without nested oversubscription. A failing replica
   stops later models after all already-started jobs in that frozen batch exit;
   it never triggers a selective rerun.
7. Run `summarize_secondary_endpoint_energy.py`. It retains every 300-frame
   series, five 60 ns block means, three replica means, and a seeded
   realization/block bootstrap interval. Fixed-residue decomposition is written
   in the prespecified order as window and five-block descriptive means only.
   It produces no frame-level p value, residue-level test or ranking, absolute
   binding-free-energy claim, affinity/potency claim, target-occupancy claim,
   efficacy claim, or ligand ranking.

## CPU node sizing and time target

The RTX 5090 is unnecessary for PB/GB. CPU count, RAM, and scratch throughput
are the limiting resources. A 128–208 logical-CPU node with 512–768 GiB RAM is
the preferred one-to-two-day target for this large receptor; a smaller node can
be cheaper but may miss that target. The three-frame canary, not a guess,
determines the final rank count and whether the selected CPU node has enough
RAM. GB is expected to finish before membrane PB; the formal runner still
preserves the fixed model order and complete three-replica inclusion.
