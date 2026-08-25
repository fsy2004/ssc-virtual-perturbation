# Analysis-configuration freeze record (pre-registration)

Date: 2026-08-21 (frozen before any production trajectory result is inspected by
this analysis pipeline; the only prior read of production data was the Codex
exploratory quick-look on rep01, which is explicitly NOT used as a basis for any
value below).

System: 8KCT–O6U nirogacestat membrane MD, three 500 ns realizations,
fixed primary window 200–500 ns, PBC tolerance 0.01 nm.

## Source of every value

All values below are fixed from (i) the fail-closed gate specification
(`POSTPROCESSING_FAIL_CLOSED_GATES.md`), (ii) the accepted experimental 8KCT
structure and its Guo 2025 Supplementary Table 3 contact record, (iii) the
constructed membrane geometry (step5_input.pdb), (iv) the frozen production
protocol (`config/production_protocol_hmr4fs_303K_v1.json`), and (v) standard
practice in the cited literature. No trajectory-derived number is used.

## 1. Structural analysis definitions (from POSTPROCESSING_FAIL_CLOSED_GATES)

| Item | Value | Source |
|---|---|---|
| Native contact cutoff | 0.45 nm | Gate 5 definition (heavy-atom pair at/below 0.45 nm in accepted 8KCT model) |
| Reference distance tolerance | 0.0001 nm | numerical exactness |
| H-bond distance cutoff | 0.35 nm | Gate 5 definition |
| H-bond angular deviation cutoff | 30 deg | Gate 5 definition |
| PBC distance-invariance tolerance | 0.01 nm | Gate 3; repeated compressed-coordinate quantization |
| Continuous-event search window | 0–500 ns | Gate 5; no gap bridging |
| Minimum continuous event duration | 5.0 ns | Gate 5: "positive duration no longer than 5 ns" (frozen at the maximum allowed) |

## 2. TM-core / pocket alignment

The alignment set is the 13 PSEN1 residues reported for nirogacestat in Guo et
al. 2025, Supplementary Table 3 (hydrogen bonds: K380, L432; van der Waals:
Y77 (TM1), V261 (TM6), L268/L271/V272/L282/I287 (TM6a), V379/L381 (TM7),
L425 (TM8), A431 (TM9)). These are experimentally resolved 8KCT residues;
no trajectory-derived residue may be added. The O6U heavy-atom set is the
35-atom non-hydrogen atom list mapped CCD→CGenFF via
`O6U_CCD_CGENFF_ATOM_CORRESPONDENCE.tsv`; atom mapping is one-to-one and
chirality preserving (identity audit: `O6U_PREPARAMETERIZATION_AUDIT.json`).

## 3. Stationarity scale floors (normalization denominator floors)

A scale floor is the minimum denominator in
`max(1.4826*MAD, floor)` for the fixed five-block stationarity audit. Floors
are set from the 2.6 Å experimental resolution and the physical size of the
measured quantity; they prevent a near-zero MAD from making a negligible
change "significant".

| Metric | Floor | Rationale |
|---|---|---|
| pocket_aligned_o6u_heavy_rmsd_nm | 0.05 | ~half the 2.6 Å cryo-EM resolution; ligand heavy-atom RMSD resolution floor |
| pocket_aligned_o6u_com_displacement_nm | 0.05 | same resolution argument |
| tm_core_ca_rmsd_nm | 0.05 | C-alpha precision at 2.6 Å |
| protein_ca_rmsd_nm | 0.05 | C-alpha precision at 2.6 Å |
| native_contact_fraction | 0.05 | 5% contact-fraction resolution |
| phosphate_peak_thickness_nm | 0.05 | phosphate peak position resolution (0.05 nm bandwidth) |
| protein_aware_area_per_lipid_nm2 | 0.05 | per-lipid area resolution (0.05 nm2) |
| cell_lateral_area_nm2_not_apl | 0.10 | box area resolution at 0.05 nm edge quantization |
| box_z_vector_length_nm | 0.05 | box-vector resolution |
| cell_volume_nm3 | 0.50 | volume resolution at 0.05 nm edges |
| protein_tilt_deg | 1.0 | tilt-angle resolution (1 deg) |
| temperature_k | 1.0 | thermostat fluctuation scale |
| pressure_bar / pressure components | 5.0 | instantaneous pressure fluctuation scale |
| potential/kinetic/total energy kj_mol | 1000.0 | total-energy fluctuation scale (system ~330k atoms) |
| density_kg_m3 | 1.0 | density resolution |
| volume_nm3 | 5.0 | volume resolution |

## 4. Fixed five-block trend/change limits

| Gate | Normalized limit | Basis |
|---|---|---|
| maximum_abs_normalized_linear_change | 3.0 | effect-size gate; a 300-ns drift beyond 3×(MAD/floor) is non-stationary |
| maximum_abs_normalized_first_last_shift | 3.0 | same |
| maximum_abs_normalized_adjacent_shift | 3.0 | same |
| maximum_abs_normalized_change_point_shift | 3.0 | same |

The 3.0 normalized limit follows the fail-closed gate specification's
"block-median linear change, first-to-last median shift, maximum adjacent-block
shift, maximum eligible two-sided change-point shift" effect-size design; the
absolute value 3.0 (robust z-scale) is chosen as the prespecified change
threshold consistent with the spike-adjudication robust-z criterion (12.0) but
much tighter because it is a mean/median drift, not a single-point spike. No
post-hoc adjustment is permitted.

## 5. Native-pose event gates (frozen before production review)

| Item | Value | Basis |
|---|---|---|
| maximum_pocket_aligned_o6u_heavy_rmsd_nm | 0.50 | Ligand RMSD above 0.5 nm from the experimental pose = substantial pocket exit; standard practice for ligand pose retention in membrane-protein MD (see Dehury 2019; Lee 2018) |
| maximum_o6u_com_displacement_nm | 0.50 | same argument for COM displacement |
| minimum_native_contact_fraction | 0.50 | fewer than 50% of the 11 frozen native contacts = contact loss |
| minimum_fraction_of_primary_frames_meeting_all_pose_gates | 0.80 | primary-window passing-frame fraction |
| ligand_egress_or_contact_loss_is_scientific_failure | true | protocol |
| failure_triggers_rerun_or_extension | false | protocol |

## 6. Thermodynamic/cell QC gates

| Item | Value | Basis |
|---|---|---|
| target_temperature_k | 303.15 | frozen production protocol |
| maximum_absolute_primary_mean_temperature_deviation_k | 5.0 | V-rescale 200–500 ns mean deviation tolerance |
| approved_primary_mean_pressure_range_bar | [-100, 100] | semi-isotropic membrane simulation without dispersion correction; 200–500 ns mean pressure window |
| approved_primary_mean_density_range_kg_m3 | [950, 1100] | CHARMM TIP3P + POPC mixture density window |
| maximum_relative_total_energy_closure_error | 1.0e-3 | numerical closure tolerance |
| maximum_absolute_pressure_trace_closure_bar | 10.0 | pressure-tensor trace closure |
| maximum_relative_orthorhombic_volume_closure_error | 1.0e-3 | box-volume closure |

## 7. Membrane QC gates (from build geometry)

| Item | Value | Basis |
|---|---|---|
| phosphate_density_bandwidth_nm | 0.05 | standard kernel density bandwidth for bilayer profiles |
| phosphate_density_grid_nm | [-4.0, 4.0, 161] | membrane-relative grid, 0.05 nm spacing, odd count |
| leaflet_hysteresis_nm | 0.2 | leaflet-assignment hysteresis (avoid z-oscillation churn) |
| hydrophobic_core_half_thickness_nm | 1.4971 | from constructed geometry (POPC acyl-carbon z p025/p975 vs membrane centre) |
| protein_xy_exclusion_nm | 4.4164 | from constructed geometry (protein heavy-atom XY r99 around membrane axis) |
| water_cluster_cutoff_nm | 0.35 | standard 3D clustering cutoff (Smith et al. 2019) |
| maximum_cumulative_leaflet_flip_events | 0 | strict: no unexplained phospholipid leaflet identity change |
| water_defect_largest_cluster_threshold | 10 | prespecified positive integer (10 waters in one persistent membrane-core cluster) |
| water_defect_persistence_frames | 5 | 5 consecutive frames |

## 8. External membrane metrics

`protein_aware_area_per_lipid` and `popc_deuterium_order_parameters` remain
`status = not_available` (hard pre-production NO-GO) because no source-hashed
validated APL@Voro v3.3/FATSLiM or gorder route with exact CHARMM36 POPC
mappings exists in this project. This record does NOT fabricate that
validation. Consequence: the membrane COMPLETE status is
`blocked_external_membrane_metrics` and the primary validation claim gate is
`blocked_inconclusive` until a validated route is supplied (or the limitation
is reported explicitly, which is the intended scientific outcome here).

## 9. Inputs bound to this freeze

| Input | SHA-256 |
|---|---|
| step5_input.pdb (trajectory atom identity) | 2cc08d30b08bd5a724eaa9538b767c943ab5ebda47f0d63cc0237d69afee4dc5 |
| 8KCT_protonated.pdb (PLIP reference) | 9b7abae2b6649c79f9688b233def9dd15ff321656fa1b8fcf98b3c5bd8a6ed38 |
| O6U_CCD_CGENFF_ATOM_CORRESPONDENCE.tsv | 62b5a9500a0c5e0c2d85eb3fa51fb4e4cb82881dafd3d0b2b38f02231d6935f5 |
| PLIP XML (run1) | fd0bedd2a3ae53fc9ffbcc49ea5dbd554b09798fa96060d70b5ecb3ec972cd61 |
| PLIP normalized JSON | 69c57604a13cb454ee55bd98b4614eefaac0527fc0a2636351393aefb874969d |
| minimized.gro | 394b53a9754dd380925f281f8f0aa013b9ea38ca9261340c132bebdf8d836829 |
