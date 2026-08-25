# Automated native-structure figure protocol

## Scientific source

The production start is an experimental 2.60 Angstrom cryo-EM heavy-atom pose, not a docking prediction. The structural figure therefore uses 8KCT/O6U directly and does not contain a docking score, redocked pose, ligand-deleted control, or second inhibitor.

The figure legend must state:

> Native nirogacestat heavy-atom pose in the 2.60 Angstrom cryo-EM structure 8KCT. Dashed annotations denote rule-based putative geometric contacts identified on the prepared static structure; they are not measured forces, interaction energies, occupancies, or binding-affinity estimates.

## Panel design

### Panel A: complex and membrane orientation

- Show the full four-subunit biological assembly as a clean cartoon.
- Show the approved membrane boundaries as two subtle planes; bulk lipids are omitted for clarity and this omission is stated.
- Show native O6U as sticks/spheres and retain its location in the PSEN1 substrate-binding pocket.
- Display resolved glycans and structural CLR/PC1 only when they aid orientation; do not create an unreadable all-atom scene.

### Panel B: native binding pocket

- Align the camera to the frozen view matrix and show O6U plus pocket residues within the preregistered distance shell.
- Highlight the Guo et al. supplementary-table residue set as literature-reported structural contacts: Tyr77, Val261, Leu268, Leu271, Val272, Leu282, Ile287, Val379, Lys380, Leu381, Leu425, Ala431, and Leu432.
- Identify Lys380 and Leu432 as literature-reported hydrogen-bond contacts only when the exact wording matches the source. Do not invent atom endpoints that the paper did not report.
- Draw a dashed segment only when an endpoint-resolved PLIP record from the prepared static structure supports it. A literature residue label and a PLIP segment are separate provenance fields.
- For the pocket panel, show only the PSEN1 cartoon context, retain every approved PLIP segment, and label only the two source-reported hydrogen-bond residues (Lys380 and Leu432). This avoids label selection after viewing while keeping the chemically most interpretable annotations legible.

### Optional post-MD panel

After production, a separate data panel may show pocket-aligned O6U RMSD and frozen native-contact occupancy for all three realizations. It uses the complete unsmoothed traces and realization-level summaries. It is not merged into the static structural-contact graphic and does not show a single favorable snapshot as representative without quantitative support.

## Contact source and limitations

PLIP is the only automated endpoint-level contact detector in this package. It identifies putative hydrogen bonds, hydrophobic contacts, pi interactions, salt bridges, water bridges, halogen contacts, and metal coordination from chemical-group and geometric rules.

- PLIP output is not a force field and does not calculate van der Waals energy, electrostatic energy, interaction strength, stability, residence time, or affinity.
- The prepared structure, hydrogen-placement method, PLIP version, command, XML/JSON output, and all endpoint atom IDs are retained and hashed.
- 8KCT contains no resolved water. The native cryo-EM figure must not show a water bridge created from a modeled solvation water.
- All detected approved interaction classes are rendered by fixed rules. Manual addition, deletion for aesthetics, or changing cutoffs after viewing the figure is prohibited. Any exclusion must be machine-recorded with a prespecified reason.

## Frozen visual language

Use a colorblind-safe Okabe-Ito-derived palette consistently across every panel:

| Element | Hex | Display |
|---|---|---|
| PSEN1 | `#0072B2` | cartoon; pocket residues as sticks |
| Nicastrin | `#009E73` | cartoon |
| APH-1A | `#E69F00` | cartoon |
| PEN-2 | `#D55E00` | cartoon |
| O6U carbon scaffold | `#CC79A7` | sticks with conventional heteroatom colors |
| Background/non-pocket protein | `#B3B3B3` | light surface or cartoon, low visual weight |
| Hydrogen bond | `#0072B2` | dashed, fixed width |
| Hydrophobic contact | `#7A7A7A` | dashed, distinct spacing from hydrogen bonds |
| Halogen/other approved geometry | `#E69F00` | dotted/dashed geometry defined in `config/figure_style.json` |

Use conventional atom colors for oxygen red, nitrogen blue, fluorine pale green, sulfur yellow, and phosphorus orange. Do not use rainbow coloring, glossy gradients, arbitrary per-panel colors, or red/green as the sole distinction. Dash spacing and width must distinguish interaction classes even in grayscale.

Labels use a sans-serif font, residue format `PSEN1 Lys380`, and no more decimal precision than the coordinate model justifies. Distances, if printed, are labeled `geometric distance` and shown consistently in Angstrom.

## Deterministic PyMOL workflow

Local executable currently available:

`C:\Users\fsy\AppData\Local\Schrodinger\PyMOL2\Scripts\pymol.exe`

Observed local version: PyMOL 3.1.1. The real production record must query and retain the complete build/version tuple again at rendering time.

Use the console entry point with `-cq`. The Windows GUI launcher `PyMOLWin.exe` detaches before the ray-traced child process has completed and is therefore prohibited for automated output validation.

`scripts/render_plip_pymol.py` must receive:

- exact 8KCT prepared-coordinate path and SHA-256;
- native O6U identity and atom mapping;
- PLIP XML or strictly validated PLIP-derived JSON;
- literature-contact record with source/table citation;
- frozen camera/view matrix and `config/figure_style.json` hash;
- output directory that does not already contain an undeclared file.

The renderer writes the generated PML, PSE, ray-traced PNG, camera matrix, contact-to-segment map, stdout/stderr, PyMOL/PLIP versions, input/output hashes, and a machine-readable manifest. `scripts/validate_figure_outputs.py` reconstructs the PML from inputs and requires exact equality.

## Output specification

- Master raster: 16 x 12 cm, 600 dpi, non-interlaced RGB/RGBA PNG with embedded physical-resolution metadata.
- White or transparent uniform background as frozen in the style file.
- Ray tracing, antialiasing, ambient occlusion, fog, depth cue, shadows, stick radii, dash widths, and camera are fixed in code.
- Retain vector text/layout separately if multi-panel assembly is performed outside PyMOL; never resample the master below target resolution before final export.
- Validate pixel dimensions, DPI metadata, CRC, channel mode, legibility at final size, clipping, color consistency, and exact interaction count.

## Hard stops

Do not publish the panel if source identity, pose provenance, PLIP parsing, endpoint mapping, camera reconstruction, version capture, image dimensions, or input/output hashes fail. Do not substitute a manually beautified image when automation fails.
