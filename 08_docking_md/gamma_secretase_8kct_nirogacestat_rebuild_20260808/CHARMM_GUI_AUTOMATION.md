# CHARMM-GUI Quick Bilayer automation

## Scope

The official Quick Bilayer API automates one PPM orientation and membrane packing/solvation job after the curated structure has passed PDB Reader. It does not decide missing residues, topology segments, termini, disulfides, glycans, protonation, ligand parameters, or structural-lipid retention.

Endpoint: `POST https://www.charmm-gui.org/api/quick_bilayer`

Official documentation: <https://www.charmm-gui.org/?doc=api&module=quickb>

The `junepark6/cgui_api` repository is retained only as a public implementation example cited by the Quick Bilayer work. Its code is not vendored. `scripts/charmmgui_api.py` is an independent, fail-closed client for the documented endpoint.

Bearer-authenticated traffic is pinned to the exact origin/base `https://www.charmm-gui.org/api`; HTTPS lookalike/custom hosts and redirects are rejected. A custom base is allowed only as an explicit HTTP `localhost`/`127.0.0.1` test and the client then prohibits a real token file/environment token. Build acceptance parses the retained submit record and requires the exact official Quick Bilayer endpoint, frozen payload, returned job ID, and current client SHA-256; a generic hash-bound JSON record is insufficient.

## Two scientific sign-offs around one automated build

1. Upload the curated 8KCT biological assembly to PDB Reader.
2. Apply and verify the decisions in `STRUCTURE_CURATION_PROTOCOL.md` and `PROTONATION_MODEL.md`.
3. Attach the independently validated neutral O6U topology/parameters and verify their hashes in the generated topology.
4. Confirm that the PDB Reader coordinates have not already been transformed by PPM; the automated build must apply orientation exactly once.
5. Save the PDB Reader report, screenshots, curated coordinates, job ID, versions, and SHA-256 manifest.

The API job is not submitted until the structure and ligand-parameter records have two-person sign-off. After the job returns, inspect and sign the PPM transformation, membrane boundaries, retained O6U/CLR/PC1/glycan components, packing, and periodic-image separation before accepting `build01`.

## Frozen Quick Bilayer payload

- PDB Reader job: one approved job ID;
- upper leaflet: `POPC=1`;
- lower leaflet: `POPC=1`;
- XY margin: 20 Angstrom;
- water boundary: 22.5 Angstrom;
- salt: 0.15 M NaCl;
- bulk membrane: symmetric pure POPC;
- PPM/orientation: `ppm=true` applies one transformation during Quick Bilayer; no PPM-transformed coordinate is submitted as input;
- `heteroatoms=true` so O6U, 3 CLR, 2 PC1, and resolved glycans are retained;
- `topologyIn=true` unless the live API response documents a changed meaning;
- `clone_job=false`;
- exactly one Quick Bilayer construction.

The 20 Angstrom margin and 22.5 Angstrom water boundary are starting values supported by the official API defaults/workflow. The packed structure must still show at least 2-3 complete lipid shells between the protein and its periodic image. If it does not, the common build is rejected before any realization exists; the margin is revised once in a new manifest version and the build is repeated.

## Safe client workflow

Store the bearer token outside the repository:

```bash
export CHARMMGUI_TOKEN_FILE=/secure/path/session.token
```

Review a redacted dry run first:

```bash
python scripts/charmmgui_api.py submit \
  --pdb-reader-jobid PDB_READER_JOB_ID \
  --system-id 8kct_nirogacestat_native \
  --build-id build01 \
  --upper 'POPC=1' \
  --lower 'POPC=1' \
  --margin 20 \
  --wdist 22.5 \
  --ion-conc 0.15 \
  --ion-type NaCl \
  --ppm \
  --heteroatoms \
  --dry-run
```

Submit only after the redacted payload matches the frozen manifest, then poll and download atomically:

```bash
python scripts/charmmgui_api.py status --jobid QUICK_BILAYER_JOB_ID
python scripts/charmmgui_api.py download \
  --jobid QUICK_BILAYER_JOB_ID \
  --output builds/8kct_nirogacestat_native/build01/charmm-gui.tgz
```

The exact client flags shown by `--help` are authoritative; preflight rejects a command transcript that does not reproduce the signed payload.

## Archive validation

A successful HTTP response is not scientific approval. Before minimization, the validator must confirm:

- safe archive paths, expected size/type, and SHA-256;
- complete GROMACS topology/includes, coordinates, staged equilibration inputs, force-field release, and generator metadata;
- exact protein segment, disulfide, glycan, O6U, CLR, PC1, POPC, water, and ion counts;
- neutral O6U topology charge and nonzero partial charges;
- exact approved ligand parameter hashes and CHARMM/GROMACS single-point energy agreement;
- no duplicated native structural lipid, ring penetration, steric packing defect, vacuum gap, or insufficient periodic-image separation;
- valid box vectors, neutral total system charge after ions, and the requested NaCl concentration within integer-ion constraints;
- visual inspection of the pocket, full protein, glycans, both leaflets, and periodic images.

The accepted build record must retain the single PPM transformation matrix, pre/post-orientation coordinate hashes, membrane boundaries, screenshots, Quick Bilayer job ID, and the SHA-256 of the orientation record. A missing record or evidence of a second orientation rejects the build.

After safe extraction, run `scripts/hash_tree_manifest.py` over the entire GROMACS directory. The manifest must enumerate every file, byte count, and SHA-256. Then run `scripts/validate_charmmgui_archive.py --manifest config/study_manifest.json --archive ... --output ... --strict`; its schema-2 report binds the build contract, archive, PDB Reader/Quick Bilayer job IDs, orientation record, extracted-tree manifest, coordinate/topology/index/MDP hashes, exact six-stage equilibration hashes, and the no-HMR audit. A PASS report from another archive or job cannot satisfy preflight.

Three stochastic realizations are created only after this single archive passes. They do not trigger three additional Quick Bilayer jobs.

## Security and provenance

The token, password, cookie, email, SSH credential, and server key are never stored in the package or logs. The client must redact authorization headers, use atomic temporary downloads, preserve the API response and payload, and hash every accepted artifact. Historical plaintext server credentials are considered compromised and must not be reused.
