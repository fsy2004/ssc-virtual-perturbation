#!/usr/bin/env python
"""Build a skin-only, donor-resolved metadata table for the integrated atlas.

The panel atlas retains library/sample identifiers, but GSE195452 contains
multiple sorted libraries per patient and several non-skin tissues. GSE249279
also contains paired epidermal/dermal libraries and reprocessed GSE181957
samples. This script resolves those structures from deposited GEO SOFT
metadata before any donor-level inference is performed.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import re
import urllib.request
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
ATLAS = ROOT / "server_archive" / "fig_atlas" / "atlas_fig.h5ad"
SAMPLE_CONDITION = ROOT / "server_archive" / "tier1" / "powered" / "sample_condition.csv"
OUT = ROOT / "04_manuscript" / "revision_20260722" / "donor_metadata"
REFERENCES = OUT / "references"
OUTPUTS = OUT / "outputs"

SOFT_SOURCES = {
    "GSE138669": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE138nnn/GSE138669/soft/GSE138669_family.soft.gz",
        "3DAC75E62ED030DC33E2D9F352EE901AEEC7A6F37AF07C40695F2D02ECE74419",
    ),
    "GSE181957": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE181nnn/GSE181957/soft/GSE181957_family.soft.gz",
        "25A2F0AA8CF96E2C89B158E1C02365565968ABAFFE6A9D04EFDF02AB887E84B9",
    ),
    "GSE195452": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE195nnn/GSE195452/soft/GSE195452_family.soft.gz",
        "8DB12B4838D31A78AE2AE6516D167C69A43D931D371994C47641C71AD86FE577",
    ),
    "GSE210395": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE210nnn/GSE210395/soft/GSE210395_family.soft.gz",
        "2E03B00DCDB56783F6E5CE977783643045BB3B901E91067ED3AA31FEE41CEE9D",
    ),
    "GSE236111": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE236nnn/GSE236111/soft/GSE236111_family.soft.gz",
        "DC69412E0E29BD8E442F9F1DD4AFDCB5C3ED4291ABDB2B488B7FCA3D710F4E1A",
    ),
    "GSE249279": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE249nnn/GSE249279/soft/GSE249279_family.soft.gz",
        "25D24ADDDDB1E2DD2512A06EE58388173DE65B3B7D1506146C63B471F69D7CD1",
    ),
    "GSE292979": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE292nnn/GSE292979/soft/GSE292979_family.soft.gz",
        "F6C8C949B58E8DD67C11F04155F4BBD34FEA2FB52C3BF947F7C745308B65C11B",
    ),
}


def log(message: str) -> None:
    print(f"[donor-metadata] {message}", flush=True)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest().upper()


def fetch_soft(accession: str, url: str, expected_hash: str) -> Path:
    REFERENCES.mkdir(parents=True, exist_ok=True)
    destination = REFERENCES / f"{accession}_family.soft.gz"
    if not destination.exists() or sha256(destination) != expected_hash:
        log(f"downloading {accession} GEO SOFT metadata")
        urllib.request.urlretrieve(url, destination)
    observed_hash = sha256(destination)
    if observed_hash != expected_hash:
        raise ValueError(
            f"SHA-256 mismatch for {accession}: expected {expected_hash}, observed {observed_hash}"
        )
    return destination


def parse_soft(path: Path, accession: str) -> pd.DataFrame:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        text = handle.read()
    rows: list[dict] = []
    for part in text.split("^SAMPLE = ")[1:]:
        lines = part.splitlines()
        gsm = lines[0].strip()
        title = ""
        source = ""
        characteristics: list[str] = []
        for line in lines[1:]:
            if line.startswith("!Sample_title = "):
                title = line.split(" = ", 1)[1].strip()
            elif line.startswith("!Sample_source_name_ch1 = "):
                source = line.split(" = ", 1)[1].strip()
            elif line.startswith("!Sample_characteristics_ch1 = "):
                characteristics.append(line.split(" = ", 1)[1].strip())

        def characteristic(prefix: str) -> str:
            prefix_lower = prefix.lower() + ":"
            for value in characteristics:
                if value.lower().startswith(prefix_lower):
                    return value.split(":", 1)[1].strip()
            return ""

        rows.append(
            {
                "accession": accession,
                "gsm": gsm,
                "title": title,
                "source_name": source,
                "tissue_soft": characteristic("tissue"),
                "patient_id_soft": characteristic("patient id"),
                "condition_soft": characteristic("condition")
                or characteristic("disease")
                or characteristic("disease state")
                or characteristic("health_status"),
            }
        )
    return pd.DataFrame(rows)


def gse249279_donor(library: str) -> str:
    if not isinstance(library, str) or not library:
        return "unknown"
    if library.startswith("SCD-"):
        donor = library.removeprefix("SCD-")
        if re.fullmatch(r"\d+YJ[12]", donor):
            donor = re.sub(r"J[12]$", "", donor)
        return donor
    return library


def build_mapping() -> tuple[pd.DataFrame, dict]:
    if not ATLAS.exists() or not SAMPLE_CONDITION.exists():
        raise FileNotFoundError("atlas_fig.h5ad or sample_condition.csv is missing")

    soft_tables = []
    reference_manifest = {}
    for accession, (url, expected_hash) in SOFT_SOURCES.items():
        path = fetch_soft(accession, url, expected_hash)
        soft_tables.append(parse_soft(path, accession))
        reference_manifest[accession] = {
            "url": url,
            "path": str(path),
            "sha256": sha256(path),
        }
    soft = pd.concat(soft_tables, ignore_index=True)

    atlas = ad.read_h5ad(ATLAS, backed="r")
    atlas_samples = pd.DataFrame({"sample": atlas.obs["sample"].astype(str).unique()})
    cells_per_sample = atlas.obs["sample"].astype(str).value_counts().rename("n_cells_atlas")
    atlas.file.close()

    condition = pd.read_csv(SAMPLE_CONDITION, dtype=str)
    if condition["sample"].duplicated().any():
        raise ValueError("sample_condition.csv contains duplicate sample identifiers")
    mapping = atlas_samples.merge(condition, on="sample", how="left", validate="one_to_one")
    mapping["n_cells_atlas"] = mapping["sample"].map(cells_per_sample).astype(int)
    mapping["tissue"] = "skin"
    mapping["donor_local"] = mapping["sample"]
    mapping["donor_mapping_source"] = "one deposited library per donor"
    mapping["library_title"] = mapping["sample"]

    # GSE195452 contains repeated CD45+/CD90+ libraries and non-skin tissues.
    gur_soft = soft[soft["accession"].eq("GSE195452")].set_index("gsm")
    is_gur = mapping["cohort"].eq("Gur_GSE195452")
    mapping.loc[is_gur, "tissue"] = (
        mapping.loc[is_gur, "sample"].map(gur_soft["tissue_soft"]).fillna("unknown").str.lower()
    )
    mapping.loc[is_gur, "donor_local"] = (
        mapping.loc[is_gur, "sample"].map(gur_soft["patient_id_soft"]).fillna("unknown")
    )
    mapping.loc[is_gur, "library_title"] = (
        mapping.loc[is_gur, "sample"].map(gur_soft["title"]).fillna(mapping.loc[is_gur, "sample"])
    )
    mapping.loc[is_gur, "donor_mapping_source"] = "GSE195452 deposited patient id"

    # GSE210395 is a PBMC study and is not part of the dermal donor analysis.
    mapping.loc[mapping["cohort"].eq("GSE210395"), "tissue"] = "pbmc"

    # GSE249279 uses internal aggregate-library identifiers in atlas.obs.
    is_249 = mapping["cohort"].eq("GSE249279")
    libraries = mapping.loc[is_249, "cond_source"].str.extract(
        r"library '([^']+)'", expand=False
    )
    mapping.loc[is_249, "library_title"] = libraries.to_numpy()
    mapping.loc[is_249, "donor_local"] = libraries.map(gse249279_donor).to_numpy()
    mapping.loc[is_249, "donor_mapping_source"] = "GSE249279 deposited library title"

    mapping["donor_id"] = mapping["cohort"].fillna("unknown") + ":" + mapping["donor_local"]
    shared_sc = mapping["donor_local"].str.fullmatch(r"SC\d+", na=False)
    mapping.loc[shared_sc, "donor_id"] = "shared:" + mapping.loc[shared_sc, "donor_local"]

    gse249_donors = set(mapping.loc[is_249, "donor_id"])
    mapping["duplicate_reprocessing"] = (
        mapping["cohort"].eq("GSE181957") & mapping["donor_id"].isin(gse249_donors)
    )
    known_condition = mapping["condition"].isin(["SSc", "HC"])
    valid_donor = ~mapping["donor_local"].isin(["", "-", "unknown"])
    mapping["analysis_keep"] = (
        known_condition
        & mapping["tissue"].eq("skin")
        & valid_donor
        & ~mapping["duplicate_reprocessing"]
    )
    mapping["exclusion_reason"] = ""
    mapping.loc[~known_condition, "exclusion_reason"] = "condition not SSc or HC"
    mapping.loc[known_condition & ~mapping["tissue"].eq("skin"), "exclusion_reason"] = (
        "non-skin tissue"
    )
    mapping.loc[known_condition & mapping["tissue"].eq("skin") & ~valid_donor, "exclusion_reason"] = (
        "deposited donor identifier unavailable"
    )
    mapping.loc[mapping["duplicate_reprocessing"], "exclusion_reason"] = (
        "GSE181957 sample reprocessed within GSE249279"
    )

    kept = mapping[mapping["analysis_keep"]].copy()
    donor_condition_n = kept.groupby("donor_id")["condition"].nunique()
    if (donor_condition_n > 1).any():
        raise ValueError("At least one resolved donor has conflicting condition labels")
    donor_cohort_n = kept.groupby("donor_id")["cohort"].nunique()
    if (donor_cohort_n > 1).any():
        raise ValueError("At least one retained donor spans multiple cohorts after deduplication")

    mapping = mapping.sort_values(["analysis_keep", "cohort", "donor_id", "sample"], ascending=[False, True, True, True])
    manifest = {
        "analysis_unit": "deposited donor/patient identifier",
        "dermal_filter": "SSc/HC skin libraries only",
        "technical_library_handling": "sum cells within donor after tissue filtering",
        "duplicate_handling": "exclude 15 GSE181957 libraries reprocessed in GSE249279",
        "atlas_shape": [int(cells_per_sample.sum()), int(len(atlas_samples))],
        "retained_samples": int(mapping["analysis_keep"].sum()),
        "retained_donors": int(kept["donor_id"].nunique()),
        "retained_donors_by_condition": kept.drop_duplicates("donor_id")["condition"].value_counts().to_dict(),
        "retained_donors_by_cohort_condition": (
            kept.drop_duplicates("donor_id")
            .groupby(["cohort", "condition"], observed=True)
            .size()
            .rename("n_donors")
            .reset_index()
            .to_dict(orient="records")
        ),
        "references": reference_manifest,
    }
    return mapping, manifest


def write_outputs(mapping: pd.DataFrame, manifest: dict) -> None:
    OUTPUTS.mkdir(parents=True, exist_ok=True)
    mapping.to_csv(OUTPUTS / "sample_to_donor.csv", index=False)

    retained = mapping[mapping["analysis_keep"]]
    donor_summary = (
        retained.groupby(["donor_id", "cohort", "condition"], observed=True)
        .agg(
            n_libraries=("sample", "nunique"),
            n_cells_atlas=("n_cells_atlas", "sum"),
            sample_ids=("sample", lambda x: ";".join(sorted(x.astype(str)))),
            library_titles=("library_title", lambda x: ";".join(sorted(x.astype(str)))),
        )
        .reset_index()
    )
    donor_summary.to_csv(OUTPUTS / "donor_summary.csv", index=False)

    audit = (
        mapping.groupby(["cohort", "tissue", "condition", "analysis_keep", "exclusion_reason"], dropna=False)
        .agg(n_libraries=("sample", "nunique"), n_cells=("n_cells_atlas", "sum"))
        .reset_index()
    )
    audit.to_csv(OUTPUTS / "sample_exclusion_audit.csv", index=False)

    duplicate = mapping[mapping["duplicate_reprocessing"]][
        ["sample", "cohort", "condition", "donor_id", "n_cells_atlas", "exclusion_reason"]
    ].copy()
    duplicate.to_csv(OUTPUTS / "duplicate_reprocessing_audit.csv", index=False)
    (OUTPUTS / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8"
    )


def main() -> None:
    mapping, manifest = build_mapping()
    write_outputs(mapping, manifest)
    log(
        f"retained {manifest['retained_donors']} donors "
        f"({manifest['retained_donors_by_condition']})"
    )
    log(f"outputs: {OUTPUTS}")


if __name__ == "__main__":
    main()
