#!/usr/bin/env python3
"""Fail-closed PLIP-to-PyMOL renderer for audited protein-ligand figures.

PLIP XML or PLIP-derived JSON is the sole source of displayed putative
contacts. The endpoint map is a reviewed copy of those immutable records, not
an independent interaction list.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path
from typing import Any, Iterable, Sequence


ALLOWED_INTERACTION_TYPES = (
    "hydrogen_bond",
    "hydrophobic_interaction",
    "pi_stack",
    "pi_cation_interaction",
    "salt_bridge",
    "water_bridge",
    "halogen_bond",
    "metal_complex",
)

INTERACTION_TAG_ALIASES = {
    "hydrogen_bond": "hydrogen_bond",
    "hbond": "hydrogen_bond",
    "hydrophobic_interaction": "hydrophobic_interaction",
    "hydrophobic_contact": "hydrophobic_interaction",
    "hydrophobic": "hydrophobic_interaction",
    "pi_stack": "pi_stack",
    "pistack": "pi_stack",
    "pi_stacking": "pi_stack",
    "pi_cation_interaction": "pi_cation_interaction",
    "pi_cation": "pi_cation_interaction",
    "pication": "pi_cation_interaction",
    "salt_bridge": "salt_bridge",
    "saltbridge": "salt_bridge",
    "water_bridge": "water_bridge",
    "waterbridge": "water_bridge",
    "halogen_bond": "halogen_bond",
    "halogenbond": "halogen_bond",
    "metal_complex": "metal_complex",
    "metalcomplex": "metal_complex",
}

PROVENANCE_PREFIX = {
    "native": "Native experimental pose",
    "modelled": "Modelled pose",
    "transferred": "Transferred pose",
}

OKABE_ITO = {
    "#000000",
    "#E69F00",
    "#56B4E9",
    "#009E73",
    "#F0E442",
    "#0072B2",
    "#D55E00",
    "#CC79A7",
}

SIMPLE_IDENTIFIER = re.compile(r"^[A-Za-z0-9_.+\-]+$")


class ValidationError(RuntimeError):
    """Raised when an input fails a scientific or provenance gate."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def write_json(path: Path, value: Any) -> None:
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True, ensure_ascii=True) + "\n",
        encoding="utf-8",
    )


def load_json_object(path: Path, label: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValidationError(f"{label} is not readable JSON: {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise ValidationError(f"{label} must contain one JSON object: {path}")
    return value


def require_file(path: Path, label: str) -> Path:
    resolved = path.expanduser().resolve()
    if not resolved.is_file() or resolved.stat().st_size == 0:
        raise ValidationError(f"{label} is missing or empty: {resolved}")
    return resolved


def require_nonempty_string(value: Any, label: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValidationError(f"{label} must be a non-empty string")
    return value.strip()


def require_finite_number(value: Any, label: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValidationError(f"{label} must be numeric, got {value!r}") from exc
    if not math.isfinite(result):
        raise ValidationError(f"{label} must be finite")
    return result


def local_tag(element: ET.Element) -> str:
    return str(element.tag).rsplit("}", 1)[-1].strip().lower()


def first_descendant_text(node: ET.Element, name: str) -> str | None:
    wanted = name.lower()
    for element in node.iter():
        if local_tag(element) == wanted and element.text is not None:
            text = element.text.strip()
            if text:
                return text
    return None


def parse_coordinate_value(value: Any, label: str) -> list[float]:
    if isinstance(value, (list, tuple)) and len(value) == 3:
        result = [require_finite_number(item, label) for item in value]
        return result
    if isinstance(value, str):
        numbers = re.findall(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?", value)
        if len(numbers) == 3:
            return [float(item) for item in numbers]
    raise ValidationError(f"{label} must contain exactly three finite coordinates")


def coordinate_from_xml(node: ET.Element, candidates: Sequence[str], label: str) -> list[float]:
    wanted = {item.lower() for item in candidates}
    for element in node.iter():
        if local_tag(element) not in wanted:
            continue
        axes: dict[str, float] = {}
        for child in element.iter():
            axis = local_tag(child)
            if axis in {"x", "y", "z"} and child.text is not None:
                axes[axis] = require_finite_number(child.text.strip(), f"{label}.{axis}")
        if set(axes) == {"x", "y", "z"}:
            return [axes["x"], axes["y"], axes["z"]]
        text = " ".join(part.strip() for part in element.itertext() if part.strip())
        try:
            return parse_coordinate_value(text, label)
        except ValidationError:
            continue
    raise ValidationError(f"PLIP interaction lacks {label} coordinates")


def canonical_interaction_type(value: Any) -> str:
    text = require_nonempty_string(value, "interaction type").lower()
    text = text.replace("-", "_").replace(" ", "_")
    result = INTERACTION_TAG_ALIASES.get(text)
    if result is None:
        raise ValidationError(f"Unsupported PLIP interaction type: {value!r}")
    return result


def interaction_identity(interaction: dict[str, Any]) -> dict[str, Any]:
    return {
        "type": interaction["type"],
        "binding_site": interaction["binding_site"],
        "protein": interaction["protein"],
        "ligand": interaction["ligand"],
        "mediators": interaction.get("mediators", []),
        "geometry": interaction.get("geometry", {}),
    }


def interaction_display_segments(interaction: dict[str, Any]) -> list[dict[str, Any]]:
    """Return only PLIP-defined endpoint segments for one interaction."""

    if interaction["type"] != "water_bridge":
        return [
            {
                "segment_id": "direct",
                "source_role": "protein",
                "target_role": "ligand",
                "source_coordinate": interaction["protein"]["coordinate"],
                "target_coordinate": interaction["ligand"]["coordinate"],
            }
        ]
    mediators = interaction.get("mediators")
    if not isinstance(mediators, list) or len(mediators) != 1:
        raise ValidationError("A PLIP water bridge must contain exactly one water mediator")
    mediator = mediators[0]
    if not isinstance(mediator, dict) or mediator.get("role") != "water":
        raise ValidationError("A PLIP water-bridge mediator must have role=water")
    coordinate = parse_coordinate_value(mediator.get("coordinate"), "water mediator")
    return [
        {
            "segment_id": "protein_water",
            "source_role": "protein",
            "target_role": "water",
            "source_coordinate": interaction["protein"]["coordinate"],
            "target_coordinate": coordinate,
        },
        {
            "segment_id": "water_ligand",
            "source_role": "water",
            "target_role": "ligand",
            "source_coordinate": coordinate,
            "target_coordinate": interaction["ligand"]["coordinate"],
        },
    ]


def assign_interaction_ids(interactions: Iterable[dict[str, Any]]) -> list[dict[str, Any]]:
    ordered = sorted((dict(item) for item in interactions), key=lambda item: stable_json(interaction_identity(item)))
    base_ids: list[str] = []
    for item in ordered:
        digest = hashlib.sha256(stable_json(interaction_identity(item)).encode("utf-8")).hexdigest()[:12]
        base_ids.append(f"{item['type']}-{digest}")
    totals = Counter(base_ids)
    seen: Counter[str] = Counter()
    output: list[dict[str, Any]] = []
    for item, base_id in zip(ordered, base_ids):
        seen[base_id] += 1
        interaction_id = base_id
        if totals[base_id] > 1:
            interaction_id = f"{base_id}-{seen[base_id]:02d}"
        item["interaction_id"] = interaction_id
        p = item["protein"]["coordinate"]
        l = item["ligand"]["coordinate"]
        item["endpoint_distance_angstrom"] = math.sqrt(sum((a - b) ** 2 for a, b in zip(p, l)))
        item["display_segment_distances_angstrom"] = {
            segment["segment_id"]: math.sqrt(
                sum(
                    (a - b) ** 2
                    for a, b in zip(segment["source_coordinate"], segment["target_coordinate"])
                )
            )
            for segment in interaction_display_segments(item)
        }
        output.append(item)
    return output


def xml_interaction(node: ET.Element, site: dict[str, str]) -> dict[str, Any]:
    interaction_type = canonical_interaction_type(local_tag(node))
    protein_coordinate = coordinate_from_xml(node, ("protcoo", "targetcoo"), "protein endpoint")
    ligand_coordinate = coordinate_from_xml(node, ("ligcoo", "metalcoo"), "ligand endpoint")

    protein = {
        "chain": first_descendant_text(node, "reschain") or "",
        "residue_number": first_descendant_text(node, "resnr") or "",
        "residue_name": first_descendant_text(node, "restype") or "",
        "atom_ids": [],
        "coordinate": protein_coordinate,
    }
    ligand = {
        "chain": first_descendant_text(node, "reschain_lig") or site["chain"],
        "residue_number": first_descendant_text(node, "resnr_lig") or site["position"],
        "residue_name": first_descendant_text(node, "restype_lig") or site["hetid"],
        "atom_ids": [],
        "coordinate": ligand_coordinate,
    }

    atom_fields = (
        "protcarbonidx",
        "ligcarbonidx",
        "donoridx",
        "acceptoridx",
        "don_idx",
        "acc_idx",
        "metal_idx",
        "target_idx",
    )
    for field in atom_fields:
        value = first_descendant_text(node, field)
        if value:
            destination = ligand["atom_ids"] if field.startswith("lig") or field == "metal_idx" else protein["atom_ids"]
            destination.append(value)
    protein["atom_ids"] = sorted(set(protein["atom_ids"]))
    ligand["atom_ids"] = sorted(set(ligand["atom_ids"]))

    geometry: dict[str, Any] = {}
    for field in (
        "dist",
        "dist_h-a",
        "dist_d-a",
        "cent_dist",
        "angle",
        "offset",
        "don_angle",
        "acc_angle",
        "water_angle",
        "coordination",
        "geometry",
        "type",
    ):
        value = first_descendant_text(node, field)
        if value:
            geometry[field] = value

    mediators: list[dict[str, Any]] = []
    if interaction_type == "water_bridge":
        mediators.append(
            {
                "role": "water",
                "coordinate": coordinate_from_xml(node, ("watercoo",), "water mediator"),
                "source": "PLIP",
            }
        )

    return {
        "type": interaction_type,
        "binding_site": site["key"],
        "protein": protein,
        "ligand": ligand,
        "mediators": mediators,
        "geometry": geometry,
        "endpoint_source": "PLIP",
    }


def select_binding_site(sites: list[dict[str, Any]], requested: str) -> dict[str, Any]:
    matches: list[dict[str, Any]] = []
    for site in sites:
        aliases = {site["id"], site["key"]}
        if requested in aliases:
            matches.append(site)
    if not matches:
        hetid_matches = [site for site in sites if site["hetid"] == requested]
        if len(hetid_matches) == 1:
            matches = hetid_matches
    if len(matches) != 1:
        available = ", ".join(sorted({site["key"] for site in sites})) or "none"
        raise ValidationError(
            f"Binding site {requested!r} is absent or ambiguous; available binding sites: {available}"
        )
    return matches[0]


def parse_plip_xml(path: Path, binding_site: str) -> dict[str, Any]:
    try:
        root = ET.parse(path).getroot()
    except (OSError, ET.ParseError) as exc:
        raise ValidationError(f"PLIP XML is malformed: {path}: {exc}") from exc

    plip_version = first_descendant_text(root, "plipversion") or root.attrib.get("plipversion")
    if not plip_version:
        raise ValidationError("PLIP XML does not record a PLIP version")

    sites: list[dict[str, Any]] = []
    for element in root.iter():
        if local_tag(element) != "bindingsite":
            continue
        site_id = element.attrib.get("id", "").strip()
        hetid = first_descendant_text(element, "hetid") or ""
        chain = first_descendant_text(element, "chain") or ""
        position = first_descendant_text(element, "position") or ""
        if not hetid:
            raise ValidationError("A PLIP XML binding site lacks its ligand HET identifier")
        key = f"{hetid}:{chain}:{position}"
        interactions_parent = next(
            (child for child in element.iter() if local_tag(child) == "interactions"),
            None,
        )
        if interactions_parent is None:
            raise ValidationError(f"PLIP binding site {key} lacks an interactions element")
        sites.append(
            {
                "id": site_id or key,
                "key": key,
                "hetid": hetid,
                "chain": chain,
                "position": position,
                "xml_element": element,
                "interactions_element": interactions_parent,
            }
        )

    selected = select_binding_site(sites, binding_site)
    raw_interactions: list[dict[str, Any]] = []
    for element in selected["interactions_element"].iter():
        tag = local_tag(element)
        if tag in INTERACTION_TAG_ALIASES:
            raw_interactions.append(xml_interaction(element, selected))

    interactions = assign_interaction_ids(raw_interactions)
    site_record = {key: selected[key] for key in ("id", "key", "hetid", "chain", "position")}
    return {
        "schema_version": "1.0",
        "source": "PLIP",
        "source_format": "xml",
        "plip_version": str(plip_version),
        "report_sha256": sha256_file(path),
        "binding_site": site_record,
        "interactions": interactions,
    }


def json_interaction(value: Any, site: dict[str, str]) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ValidationError("Each PLIP JSON interaction must be an object")
    interaction_type = canonical_interaction_type(value.get("type"))
    protein_coordinate = parse_coordinate_value(
        value.get("protcoo", value.get("targetcoo")), "protein endpoint"
    )
    ligand_coordinate = parse_coordinate_value(
        value.get("ligcoo", value.get("metalcoo")), "ligand endpoint"
    )

    def text(name: str, fallback: str = "") -> str:
        raw = value.get(name, fallback)
        return "" if raw is None else str(raw).strip()

    protein_ids = [
        text(field)
        for field in ("protcarbonidx", "donoridx", "acceptoridx", "don_idx", "acc_idx", "target_idx")
        if text(field)
    ]
    ligand_ids = [text(field) for field in ("ligcarbonidx", "metal_idx") if text(field)]
    geometry = {
        key: value[key]
        for key in (
            "dist",
            "dist_h-a",
            "dist_d-a",
            "cent_dist",
            "angle",
            "offset",
            "don_angle",
            "acc_angle",
            "water_angle",
            "coordination",
            "geometry",
        )
        if key in value
    }
    mediators: list[dict[str, Any]] = []
    if interaction_type == "water_bridge":
        mediators.append(
            {
                "role": "water",
                "coordinate": parse_coordinate_value(value.get("watercoo"), "water mediator"),
                "source": "PLIP",
            }
        )
    return {
        "type": interaction_type,
        "binding_site": site["key"],
        "protein": {
            "chain": text("reschain"),
            "residue_number": text("resnr"),
            "residue_name": text("restype"),
            "atom_ids": sorted(set(protein_ids)),
            "coordinate": protein_coordinate,
        },
        "ligand": {
            "chain": text("reschain_lig", site["chain"]),
            "residue_number": text("resnr_lig", site["position"]),
            "residue_name": text("restype_lig", site["hetid"]),
            "atom_ids": sorted(set(ligand_ids)),
            "coordinate": ligand_coordinate,
        },
        "mediators": mediators,
        "geometry": geometry,
        "endpoint_source": "PLIP",
    }


def parse_plip_json(path: Path, binding_site: str) -> dict[str, Any]:
    data = load_json_object(path, "PLIP JSON")
    source = require_nonempty_string(data.get("source"), "PLIP JSON source")
    if not source.upper().startswith("PLIP"):
        raise ValidationError("JSON interaction report is not explicitly attributable to PLIP")
    plip_version = require_nonempty_string(data.get("plip_version"), "PLIP JSON version")
    values = data.get("bindingsites")
    if not isinstance(values, list):
        raise ValidationError("PLIP JSON bindingsites must be a list")

    sites: list[dict[str, Any]] = []
    for value in values:
        if not isinstance(value, dict):
            raise ValidationError("Each PLIP JSON binding site must be an object")
        ligand = value.get("ligand")
        if not isinstance(ligand, dict):
            raise ValidationError("Each PLIP JSON binding site must contain a ligand object")
        hetid = require_nonempty_string(ligand.get("hetid"), "binding-site ligand hetid")
        chain = str(ligand.get("chain", "")).strip()
        position = str(ligand.get("position", "")).strip()
        key = f"{hetid}:{chain}:{position}"
        interactions = value.get("interactions")
        if not isinstance(interactions, list):
            raise ValidationError(f"PLIP JSON binding site {key} interactions must be a list")
        sites.append(
            {
                "id": str(value.get("id", key)).strip() or key,
                "key": key,
                "hetid": hetid,
                "chain": chain,
                "position": position,
                "json_interactions": interactions,
            }
        )

    selected = select_binding_site(sites, binding_site)
    interactions = assign_interaction_ids(
        json_interaction(value, selected) for value in selected["json_interactions"]
    )
    site_record = {key: selected[key] for key in ("id", "key", "hetid", "chain", "position")}
    return {
        "schema_version": "1.0",
        "source": "PLIP",
        "source_format": "json",
        "plip_version": plip_version,
        "report_sha256": sha256_file(path),
        "binding_site": site_record,
        "interactions": interactions,
    }


def parse_plip_report(path: Path, binding_site: str) -> dict[str, Any]:
    suffix = path.suffix.lower()
    if suffix == ".xml":
        return parse_plip_xml(path, binding_site)
    if suffix == ".json":
        return parse_plip_json(path, binding_site)
    raise ValidationError("PLIP report must have an .xml or .json extension")


def endpoint_template(normalized: dict[str, Any], structure_hash: str) -> dict[str, Any]:
    records = []
    for interaction in normalized["interactions"]:
        records.append(
            {
                "interaction_id": interaction["interaction_id"],
                "type": interaction["type"],
                "endpoint_source": "PLIP",
                "protein": interaction["protein"],
                "ligand": interaction["ligand"],
                "mediators": interaction.get("mediators", []),
                "review_status": "pending",
                "review_note": "",
            }
        )
    return {
        "schema_version": "1.0",
        "review_status": "pending",
        "endpoint_source": "PLIP",
        "structure_sha256": structure_hash,
        "plip_report_sha256": normalized["report_sha256"],
        "binding_site": normalized["binding_site"],
        "interactions": records,
    }


def immutable_endpoint_record(interaction: dict[str, Any]) -> dict[str, Any]:
    return {
        "interaction_id": interaction["interaction_id"],
        "type": interaction["type"],
        "endpoint_source": "PLIP",
        "protein": interaction["protein"],
        "ligand": interaction["ligand"],
        "mediators": interaction.get("mediators", []),
    }


def validate_endpoint_map(
    path: Path,
    normalized: dict[str, Any],
    structure_hash: str,
) -> dict[str, Any]:
    data = load_json_object(path, "endpoint map")
    if data.get("schema_version") != "1.0":
        raise ValidationError("endpoint map schema_version must be 1.0")
    if data.get("review_status") != "approved":
        raise ValidationError("endpoint map review_status must be approved")
    if data.get("endpoint_source") != "PLIP":
        raise ValidationError("endpoint map endpoint_source must be PLIP")
    if str(data.get("structure_sha256", "")).lower() != structure_hash:
        raise ValidationError("endpoint map structure SHA-256 does not match the coordinate file")
    if str(data.get("plip_report_sha256", "")).lower() != normalized["report_sha256"]:
        raise ValidationError("endpoint map PLIP-report SHA-256 does not match")
    if stable_json(data.get("binding_site")) != stable_json(normalized["binding_site"]):
        raise ValidationError("endpoint map binding-site identity differs from PLIP")
    values = data.get("interactions")
    if not isinstance(values, list):
        raise ValidationError("endpoint map interactions must be a list")

    expected = {
        item["interaction_id"]: immutable_endpoint_record(item)
        for item in normalized["interactions"]
    }
    observed: dict[str, dict[str, Any]] = {}
    for value in values:
        if not isinstance(value, dict):
            raise ValidationError("Each endpoint-map interaction must be an object")
        interaction_id = require_nonempty_string(value.get("interaction_id"), "endpoint interaction_id")
        if interaction_id in observed:
            raise ValidationError(f"Duplicate endpoint-map interaction ID: {interaction_id}")
        if value.get("review_status") != "approved":
            raise ValidationError(f"Endpoint {interaction_id} is not approved")
        immutable = {
            key: value.get(key)
            for key in (
                "interaction_id",
                "type",
                "endpoint_source",
                "protein",
                "ligand",
                "mediators",
            )
        }
        observed[interaction_id] = immutable

    if set(observed) != set(expected):
        missing = sorted(set(expected) - set(observed))
        extras = sorted(set(observed) - set(expected))
        raise ValidationError(f"Endpoint-map interaction set differs from PLIP; missing={missing}, extras={extras}")
    for interaction_id in sorted(expected):
        if stable_json(observed[interaction_id]) != stable_json(expected[interaction_id]):
            raise ValidationError(f"Endpoint {interaction_id} differs from its immutable PLIP record")
    return data


def resolve_audit_record(provenance_path: Path, record_value: Any, label: str) -> Path:
    record = Path(require_nonempty_string(record_value, label))
    if not record.is_absolute():
        record = provenance_path.parent / record
    return require_file(record, label)


def validate_provenance(
    path: Path,
    structure_hash: str,
    normalized: dict[str, Any],
) -> dict[str, Any]:
    data = load_json_object(path, "pose provenance")
    if data.get("schema_version") != "1.0":
        raise ValidationError("pose provenance schema_version must be 1.0")
    if data.get("review_status") != "approved":
        raise ValidationError("pose provenance review_status must be approved")
    provenance = require_nonempty_string(data.get("pose_provenance"), "pose_provenance").lower()
    if provenance not in PROVENANCE_PREFIX:
        raise ValidationError("pose_provenance must be native, modelled, or transferred")
    display_label = require_nonempty_string(data.get("display_label"), "display_label")
    structure_id = require_nonempty_string(data.get("structure_id"), "structure_id")
    ligand_id = require_nonempty_string(data.get("ligand_id"), "ligand_id")
    ligand_formal_charge = require_finite_number(
        data.get("ligand_formal_charge"), "ligand_formal_charge"
    )
    if not ligand_formal_charge.is_integer():
        raise ValidationError("ligand_formal_charge must be an integer")
    ligand_formal_charge_int = int(ligand_formal_charge)
    if provenance != "native" or structure_id.upper() != "8KCT" or ligand_id.upper() != "O6U":
        raise ValidationError("this frozen renderer accepts only the native 8KCT-O6U structure")
    if structure_id.upper() == "8KCT" and ligand_id.upper() == "O6U" and ligand_formal_charge_int != 0:
        raise ValidationError("The frozen 8KCT-O6U workflow accepts only neutral O6U (formal charge 0)")
    ligand_chain = str(data.get("ligand_chain", "")).strip()
    ligand_position = str(data.get("ligand_residue_number", "")).strip()
    if ligand_id != normalized["binding_site"]["hetid"]:
        raise ValidationError("provenance ligand_id differs from the PLIP binding site")
    if ligand_chain != normalized["binding_site"]["chain"]:
        raise ValidationError("provenance ligand_chain differs from the PLIP binding site")
    if ligand_position != normalized["binding_site"]["position"]:
        raise ValidationError("provenance ligand_residue_number differs from the PLIP binding site")
    if str(data.get("structure_file_sha256", "")).lower() != structure_hash:
        raise ValidationError("provenance structure SHA-256 does not match the coordinate file")
    if str(data.get("plip_report_sha256", "")).lower() != normalized["report_sha256"]:
        raise ValidationError("provenance PLIP-report SHA-256 does not match")
    if provenance != "native" and re.search(
        r"\b(?:native|experimental|crystal|cryo[- ]?em)\b", display_label, re.IGNORECASE
    ):
        raise ValidationError("A modelled or transferred display_label contains an experimental-pose claim")

    if provenance == "native":
        source_structure_id = require_nonempty_string(
            data.get("source_structure_id"), "native source_structure_id"
        )
        if source_structure_id.upper() != structure_id.upper():
            raise ValidationError("native source_structure_id must equal structure_id")
        experimental_method = require_nonempty_string(
            data.get("experimental_method"), "native experimental_method"
        )
        if structure_id.upper() == "8KCT" and not re.search(
            r"(?:electron\s+microscopy|cryo[- ]?em)", experimental_method, re.IGNORECASE
        ):
            raise ValidationError("8KCT provenance must identify electron microscopy or cryo-EM")
        resolution = require_finite_number(data.get("resolution_angstrom"), "native resolution_angstrom")
        if abs(resolution - 2.60) > 0.01:
            raise ValidationError("8KCT provenance resolution must be 2.60 Angstrom")
        if not all(token in display_label.lower() for token in ("8kct", "native")) or not re.search(
            r"cryo[- ]?em", display_label, re.IGNORECASE
        ):
            raise ValidationError("8KCT display_label must visibly identify the native cryo-EM pose")
        require_nonempty_string(data.get("entry_version"), "native entry_version")
    elif provenance == "modelled":
        require_nonempty_string(data.get("model_method"), "model_method")
        require_nonempty_string(data.get("model_method_version"), "model_method_version")
        record = resolve_audit_record(path, data.get("model_record"), "model_record")
        if sha256_file(record) != str(data.get("model_record_sha256", "")).lower():
            raise ValidationError("model_record SHA-256 does not match provenance")
    else:
        source_structure_id = require_nonempty_string(
            data.get("source_structure_id"), "transferred source_structure_id"
        )
        target_structure_id = require_nonempty_string(
            data.get("target_structure_id"), "transferred target_structure_id"
        )
        if target_structure_id.upper() != structure_id.upper():
            raise ValidationError("transferred target_structure_id must equal structure_id")
        if source_structure_id.upper() == target_structure_id.upper():
            raise ValidationError("transferred source and target structure IDs must differ")
        require_nonempty_string(data.get("alignment_selection"), "alignment_selection")
        if require_finite_number(data.get("alignment_rmsd_angstrom"), "alignment_rmsd_angstrom") < 0:
            raise ValidationError("alignment_rmsd_angstrom cannot be negative")
        record = resolve_audit_record(path, data.get("transformation_record"), "transformation_record")
        if sha256_file(record) != str(data.get("transformation_record_sha256", "")).lower():
            raise ValidationError("transformation_record SHA-256 does not match provenance")

    data = dict(data)
    data["pose_provenance"] = provenance
    data["ligand_formal_charge"] = ligand_formal_charge_int
    data["visible_label"] = f"{PROVENANCE_PREFIX[provenance]}: {display_label}"
    return data


def validate_hex_color(value: Any, label: str) -> str:
    color = require_nonempty_string(value, label).upper()
    if not re.fullmatch(r"#[0-9A-F]{6}", color):
        raise ValidationError(f"{label} must be a #RRGGBB color")
    return color


def validate_style(path: Path) -> dict[str, Any]:
    data = load_json_object(path, "figure style")
    if data.get("schema_version") != "1.0":
        raise ValidationError("figure style schema_version must be 1.0")
    if data.get("palette_name") != "Okabe-Ito":
        raise ValidationError("figure style palette_name must be Okabe-Ito")
    if data.get("display_all_plip_contacts") is not True:
        raise ValidationError("figure style must display all PLIP contacts")
    palette = data.get("palette")
    if not isinstance(palette, dict):
        raise ValidationError("figure style palette must be an object")
    palette_values = {validate_hex_color(value, f"palette.{key}") for key, value in palette.items()}
    if palette_values != OKABE_ITO:
        raise ValidationError("figure style must retain the complete canonical Okabe-Ito palette")

    styles = data.get("interaction_styles")
    if not isinstance(styles, dict) or set(styles) != set(ALLOWED_INTERACTION_TYPES):
        raise ValidationError("figure style must define every supported PLIP interaction type exactly once")
    pattern_names: set[str] = set()
    pattern_geometries: set[tuple[float, float, float]] = set()
    for interaction_type in ALLOWED_INTERACTION_TYPES:
        item = styles[interaction_type]
        if not isinstance(item, dict):
            raise ValidationError(f"interaction style {interaction_type} must be an object")
        color = validate_hex_color(item.get("color"), f"{interaction_type}.color")
        if color not in palette_values:
            raise ValidationError(f"interaction style {interaction_type} is outside the Okabe-Ito palette")
        pattern_name = require_nonempty_string(item.get("pattern_name"), f"{interaction_type}.pattern_name")
        dash_gap = require_finite_number(item.get("dash_gap"), f"{interaction_type}.dash_gap")
        dash_length = require_finite_number(item.get("dash_length"), f"{interaction_type}.dash_length")
        dash_width = require_finite_number(item.get("dash_width"), f"{interaction_type}.dash_width")
        if min(dash_gap, dash_length, dash_width) <= 0:
            raise ValidationError(f"interaction style {interaction_type} dash settings must be positive")
        if pattern_name in pattern_names:
            raise ValidationError("interaction pattern_name values must be unique for dual encoding")
        geometry_signature = (dash_gap, dash_length, dash_width)
        if geometry_signature in pattern_geometries:
            raise ValidationError("rendered interaction dash geometries must be unique for dual encoding")
        pattern_names.add(pattern_name)
        pattern_geometries.add(geometry_signature)

    render = data.get("render")
    if not isinstance(render, dict):
        raise ValidationError("figure style render must be an object")
    if int(render.get("dpi", 0)) != 600:
        raise ValidationError("formal master render must be exactly 600 dpi")
    if render.get("ray_trace") is not True:
        raise ValidationError("formal master render must use ray tracing")
    for key in ("width_cm", "height_cm"):
        if require_finite_number(render.get(key), f"render.{key}") <= 0:
            raise ValidationError(f"render.{key} must be positive")

    structure = data.get("structure")
    labels = data.get("labels")
    camera = data.get("camera")
    if not isinstance(structure, dict) or not isinstance(labels, dict) or not isinstance(camera, dict):
        raise ValidationError("figure style structure, labels, and camera must be objects")
    subunit_colors = structure.get("subunit_colors")
    if not isinstance(subunit_colors, dict) or set(subunit_colors) != {"A", "B", "C", "D"}:
        raise ValidationError("structure.subunit_colors must define exact 8KCT author chains A-D")
    for chain_id, color in subunit_colors.items():
        if validate_hex_color(color, f"structure.subunit_colors.{chain_id}") not in palette_values:
            raise ValidationError(f"8KCT chain {chain_id} color is outside the Okabe-Ito palette")
    if int(labels.get("max_residue_labels", -1)) < 0:
        raise ValidationError("labels.max_residue_labels cannot be negative")
    if camera.get("mode") != "deterministic_orient_then_zoom":
        raise ValidationError("camera.mode must be deterministic_orient_then_zoom")
    validate_hex_color(data.get("background"), "background")
    return data


def validate_camera_input(path: Path, structure_hash: str, style_hash: str) -> list[float]:
    data = load_json_object(path, "camera record")
    values = data.get("view_matrix")
    if not isinstance(values, list) or len(values) != 18:
        raise ValidationError("camera record view_matrix must contain 18 values")
    matrix = [require_finite_number(value, "camera view value") for value in values]
    if str(data.get("structure_sha256", "")).lower() != structure_hash:
        raise ValidationError("camera record structure SHA-256 does not match")
    if str(data.get("style_sha256", "")).lower() != style_hash:
        raise ValidationError("camera record style SHA-256 does not match")
    return matrix


def sanitize_object_name(value: str) -> str:
    name = re.sub(r"[^A-Za-z0-9_]", "_", value)
    if not name or not name[0].isalpha():
        name = f"obj_{name}"
    return name[:120]


def hex_to_rgb(value: str) -> list[float]:
    color = value.lstrip("#")
    return [int(color[index : index + 2], 16) / 255.0 for index in (0, 2, 4)]


def safe_selection_token(value: str, label: str, allow_empty: bool = False) -> str:
    if allow_empty and value == "":
        return value
    if not SIMPLE_IDENTIFIER.fullmatch(value):
        raise ValidationError(f"{label} contains unsupported selection characters: {value!r}")
    return value


def render_object_name(interaction_id: str, segment_id: str = "direct") -> str:
    suffix = "" if segment_id == "direct" else f"_{segment_id}"
    return sanitize_object_name(f"contact_{interaction_id}{suffix}")


def build_pml(
    structure: Path,
    normalized: dict[str, Any],
    provenance: dict[str, Any],
    style: dict[str, Any],
    style_hash: str,
    structure_hash: str,
    pml_path: Path,
    pse_path: Path,
    png_path: Path,
    camera_output: Path,
    camera_matrix: list[float] | None,
) -> str:
    site = normalized["binding_site"]
    hetid = safe_selection_token(site["hetid"], "binding-site hetid")
    chain = safe_selection_token(site["chain"], "binding-site chain", allow_empty=True)
    position = safe_selection_token(site["position"], "binding-site position")
    ligand_selection = f"complex_input and resn {hetid} and resi {position}"
    if chain:
        ligand_selection += f" and chain {chain}"

    render = style["render"]
    width_px = int(round(float(render["width_cm"]) / 2.54 * int(render["dpi"])))
    height_px = int(round(float(render["height_cm"]) / 2.54 * int(render["dpi"])))
    colors = {
        "background": style["background"],
        "protein": style["structure"]["protein_cartoon_color"],
        "pocket": style["structure"]["pocket_carbon_color"],
        "ligand": style["structure"]["ligand_carbon_color"],
        "oxygen": style["structure"]["oxygen_color"],
        "nitrogen": style["structure"]["nitrogen_color"],
        "sulfur": style["structure"]["sulfur_color"],
        "halogen": style["structure"]["halogen_color"],
        "label": style["labels"]["color"],
        "label_outline": style["labels"]["outline_color"],
    }
    color_definitions = {name: hex_to_rgb(value) for name, value in colors.items()}
    for chain_id, value in style["structure"]["subunit_colors"].items():
        color_definitions[f"subunit_{chain_id}"] = hex_to_rgb(value)
    for interaction_type, settings in style["interaction_styles"].items():
        color_definitions[f"interaction_{interaction_type}"] = hex_to_rgb(settings["color"])

    lines = [
        "python",
        "from pymol import cmd",
        "from pathlib import Path",
        "import json",
        "cmd.reinitialize()",
    ]
    for name in sorted(color_definitions):
        lines.append(f"cmd.set_color({name!r}, {color_definitions[name]!r})")
    lines.extend(
        [
            f"cmd.load({structure.as_posix()!r}, 'complex_input')",
            "cmd.hide('everything', 'all')",
            "cmd.show('cartoon', 'complex_input and polymer.protein and chain B')",
            "cmd.color('subunit_B', 'complex_input and polymer.protein and chain B')",
            f"cmd.set('cartoon_transparency', {float(style['structure']['protein_cartoon_transparency'])!r}, 'complex_input and polymer.protein and chain B')",
            f"cmd.select('figure_ligand', {ligand_selection!r})",
            "if cmd.count_atoms('figure_ligand') == 0: raise RuntimeError('Ligand selection resolved to zero atoms')",
            f"cmd.select('figure_pocket', 'byres ((complex_input and polymer.protein) within {float(style['structure']['pocket_radius_angstrom']):.6f} of figure_ligand)')",
            "cmd.show('sticks', 'figure_pocket or figure_ligand')",
            f"cmd.set('stick_radius', {float(style['structure']['stick_radius'])!r})",
            "cmd.color('pocket', 'figure_pocket and elem C')",
            "cmd.color('ligand', 'figure_ligand and elem C')",
            "cmd.color('oxygen', '(figure_pocket or figure_ligand) and elem O')",
            "cmd.color('nitrogen', '(figure_pocket or figure_ligand) and elem N')",
            "cmd.color('sulfur', '(figure_pocket or figure_ligand) and elem S')",
            "cmd.color('halogen', '(figure_pocket or figure_ligand) and (elem F or elem Cl or elem Br or elem I)')",
            "cmd.set('dash_as_cylinders', 1)",
            "cmd.set('dash_round_ends', 1)",
        ]
    )

    for interaction in normalized["interactions"]:
        interaction_id = interaction["interaction_id"]
        endpoint_names = {
            "protein": sanitize_object_name(f"endpoint_protein_{interaction_id}"),
            "ligand": sanitize_object_name(f"endpoint_ligand_{interaction_id}"),
        }
        settings = style["interaction_styles"][interaction["type"]]
        color_name = f"interaction_{interaction['type']}"
        endpoint_coordinates = {
            "protein": interaction["protein"]["coordinate"],
            "ligand": interaction["ligand"]["coordinate"],
        }
        if interaction["type"] == "water_bridge":
            mediator = interaction["mediators"][0]
            endpoint_names["water"] = sanitize_object_name(f"endpoint_water_{interaction_id}")
            endpoint_coordinates["water"] = mediator["coordinate"]
        for role in sorted(endpoint_names):
            endpoint_name = endpoint_names[role]
            lines.extend(
                [
                    f"cmd.pseudoatom({endpoint_name!r}, pos={endpoint_coordinates[role]!r})",
                    f"cmd.show('spheres', {endpoint_name!r})",
                    f"cmd.set('sphere_scale', {float(style['structure']['endpoint_sphere_scale'])!r}, {endpoint_name!r})",
                    f"cmd.color({color_name!r}, {endpoint_name!r})",
                ]
            )
        for segment in interaction_display_segments(interaction):
            object_name = render_object_name(interaction_id, segment["segment_id"])
            source_name = endpoint_names[segment["source_role"]]
            target_name = endpoint_names[segment["target_role"]]
            lines.extend(
                [
                    f"cmd.distance({object_name!r}, {source_name!r}, {target_name!r})",
                    f"cmd.hide('labels', {object_name!r})",
                    f"cmd.set('dash_color', {color_name!r}, {object_name!r})",
                    f"cmd.set('dash_gap', {float(settings['dash_gap'])!r}, {object_name!r})",
                    f"cmd.set('dash_length', {float(settings['dash_length'])!r}, {object_name!r})",
                    f"cmd.set('dash_width', {float(settings['dash_width'])!r}, {object_name!r})",
                ]
            )

    unique_residues: dict[tuple[str, str, str], list[float]] = {}
    prioritized_interactions = sorted(
        normalized["interactions"],
        key=lambda item: (item["type"] != "hydrogen_bond", item["interaction_id"]),
    )
    for interaction in prioritized_interactions:
        protein = interaction["protein"]
        key = (protein["chain"], protein["residue_name"], protein["residue_number"])
        unique_residues.setdefault(key, protein["coordinate"])
    label_limit = int(style["labels"]["max_residue_labels"])
    ligand_contact_center = [
        sum(item["ligand"]["coordinate"][axis] for item in normalized["interactions"])
        / len(normalized["interactions"])
        for axis in range(3)
    ]
    for index, (key, coordinate) in enumerate(list(unique_residues.items())[:label_limit], start=1):
        chain_value, residue_name, residue_number = key
        label_text = f"{chain_value or '?'}:{residue_name}{residue_number}"
        label_name = f"residue_label_{index:02d}"
        direction = [coordinate[axis] - ligand_contact_center[axis] for axis in range(3)]
        norm = math.sqrt(sum(value * value for value in direction)) or 1.0
        label_coordinate = [
            coordinate[axis] + 2.4 * direction[axis] / norm for axis in range(3)
        ]
        lines.extend(
            [
                f"cmd.pseudoatom({label_name!r}, pos={label_coordinate!r}, label={label_text!r})",
                f"cmd.set('label_font_id', {int(style['labels']['font_id'])}, {label_name!r})",
                f"cmd.set('label_size', {float(style['labels']['point_size'])!r}, {label_name!r})",
                f"cmd.set('label_color', 'label', {label_name!r})",
                f"cmd.set('label_outline_color', 'label_outline', {label_name!r})",
            ]
        )

    lines.extend(
        [
            "banner_position = list(cmd.centerofmass('figure_ligand'))",
            f"cmd.pseudoatom('provenance_banner', pos=banner_position, label={provenance['visible_label']!r})",
            f"cmd.set('label_font_id', {int(style['labels']['font_id'])}, 'provenance_banner')",
            f"cmd.set('label_size', {float(style['labels']['provenance_point_size'])!r}, 'provenance_banner')",
            "cmd.set('label_color', 'label', 'provenance_banner')",
            "cmd.set('label_outline_color', 'label_outline', 'provenance_banner')",
            f"cmd.set('label_position', {list(style['labels']['provenance_offset_angstrom'])!r}, 'provenance_banner')",
            "cmd.bg_color('background')",
            f"cmd.set('orthoscopic', {1 if style['camera']['orthoscopic'] else 0})",
            f"cmd.set('antialias', {int(render['antialias'])})",
            f"cmd.set('ray_shadows', {1 if render['ray_shadows'] else 0})",
            f"cmd.set('ray_opaque_background', {1 if render['opaque_background'] else 0})",
            "cmd.set('spec_reflect', 0.0)",
        ]
    )
    if camera_matrix is None:
        lines.extend(
            [
                "cmd.orient('figure_pocket or figure_ligand')",
                f"cmd.zoom('figure_pocket or figure_ligand', buffer={float(style['camera']['zoom_buffer_angstrom'])!r})",
            ]
        )
    else:
        lines.append(f"cmd.set_view({tuple(camera_matrix)!r})")

    camera_payload_literal = {
        "schema_version": "1.0",
        "camera_source": "provided" if camera_matrix is not None else "deterministic_orient_then_zoom",
        "structure_sha256": structure_hash,
        "style_sha256": style_hash,
        "width_px": width_px,
        "height_px": height_px,
        "dpi": int(render["dpi"]),
    }
    lines.extend(
        [
            f"camera_payload = {camera_payload_literal!r}",
            "camera_payload['view_matrix'] = [float(value) for value in cmd.get_view()]",
            "camera_payload['pymol_version'] = list(cmd.get_version())",
            f"Path({camera_output.as_posix()!r}).write_text(json.dumps(camera_payload, indent=2, sort_keys=True) + '\\n', encoding='utf-8')",
            f"cmd.viewport({width_px}, {height_px})",
            f"cmd.save({pse_path.as_posix()!r})",
            f"cmd.png({png_path.as_posix()!r}, width={width_px}, height={height_px}, dpi={int(render['dpi'])}, ray=1, quiet=0)",
            "cmd.quit()",
            "python end",
            "",
        ]
    )
    return "\n".join(lines)


def resolve_executable(value: str) -> Path:
    candidate = Path(value).expanduser()
    if candidate.is_file():
        return candidate.resolve()
    located = shutil.which(value)
    if located:
        return Path(located).resolve()
    raise ValidationError(f"PyMOL executable was not found: {value}")


def output_paths(output_dir: Path, stem: str) -> dict[str, Path]:
    return {
        "normalized_interactions": output_dir / f"{stem}.interactions.normalized.json",
        "pml": output_dir / f"{stem}.pml",
        "pse": output_dir / f"{stem}.pse",
        "png": output_dir / f"{stem}.600dpi.png",
        "camera": output_dir / f"{stem}.camera.json",
        "stdout_log": output_dir / f"{stem}.pymol.stdout.log",
        "stderr_log": output_dir / f"{stem}.pymol.stderr.log",
        "manifest": output_dir / f"{stem}.manifest.json",
    }


def ensure_outputs_available(paths: dict[str, Path], overwrite: bool) -> None:
    existing = [path for path in paths.values() if path.exists()]
    if existing and not overwrite:
        joined = ", ".join(str(path) for path in existing)
        raise ValidationError(f"Output files already exist; use --overwrite to replace this exact set: {joined}")
    if overwrite:
        for path in existing:
            if not path.is_file():
                raise ValidationError(f"Declared output path exists but is not a file: {path}")
            path.unlink()


def build_manifest(
    structure: Path,
    report: Path,
    provenance_path: Path,
    endpoint_path: Path,
    style_path: Path,
    normalized: dict[str, Any],
    provenance: dict[str, Any],
    binding_site: str,
    outputs: dict[str, Path],
    pymol_executable: Path,
    command: list[str],
) -> dict[str, Any]:
    output_records = {
        key: {
            "path": str(path.resolve()),
            "sha256": sha256_file(path),
            "bytes": path.stat().st_size,
        }
        for key, path in sorted(outputs.items())
        if key != "manifest"
    }
    style = load_json_object(style_path, "figure style")
    render = style["render"]
    width_px = int(round(float(render["width_cm"]) / 2.54 * int(render["dpi"])))
    height_px = int(round(float(render["height_cm"]) / 2.54 * int(render["dpi"])))
    return {
        "schema_version": "figure-render-manifest-1.0",
        "binding_site_request": binding_site,
        "binding_site": normalized["binding_site"],
        "plip_version": normalized["plip_version"],
        "pose_provenance": provenance["pose_provenance"],
        "structure_id": provenance["structure_id"],
        "ligand_id": provenance["ligand_id"],
        "ligand_formal_charge": provenance["ligand_formal_charge"],
        "visible_provenance_label": provenance["visible_label"],
        "interaction_count": len(normalized["interactions"]),
        "interaction_ids": [item["interaction_id"] for item in normalized["interactions"]],
        "render": {"dpi": int(render["dpi"]), "width_px": width_px, "height_px": height_px},
        "inputs": {
            "structure": {"path": str(structure), "sha256": sha256_file(structure)},
            "plip_report": {"path": str(report), "sha256": sha256_file(report)},
            "provenance": {"path": str(provenance_path), "sha256": sha256_file(provenance_path)},
            "endpoint_map": {"path": str(endpoint_path), "sha256": sha256_file(endpoint_path)},
            "style": {"path": str(style_path), "sha256": sha256_file(style_path)},
        },
        "pymol": {"executable": str(pymol_executable), "command": command, "returncode": 0},
        "outputs": output_records,
    }


def parser() -> argparse.ArgumentParser:
    package_root = Path(__file__).resolve().parent.parent
    result = argparse.ArgumentParser(
        description=(
            "Render an audited PLIP contact panel with deterministic PyMOL settings. "
            "PLIP XML/JSON is the only contact source."
        )
    )
    result.add_argument("--structure", required=True, type=Path, help="Local PDB or mmCIF coordinate file")
    result.add_argument("--plip-report", required=True, type=Path, help="Archived PLIP XML or PLIP-derived JSON")
    result.add_argument("--binding-site", required=True, help="Exact PLIP binding-site ID or HET:CHAIN:POSITION key")
    result.add_argument("--provenance", type=Path, help="Approved pose-provenance JSON")
    result.add_argument("--endpoint-map", type=Path, help="Approved endpoint-map JSON")
    result.add_argument(
        "--style",
        type=Path,
        default=package_root / "config" / "figure_style.json",
        help="Versioned figure-style JSON",
    )
    result.add_argument("--output-dir", type=Path, help="Directory for PML, PSE, PNG, camera, logs, and manifest")
    result.add_argument("--output-stem", default="plip_contacts", help="Deterministic output filename stem")
    result.add_argument("--camera", type=Path, help="Previously saved compatible camera JSON")
    result.add_argument("--pymol-executable", default="pymol", help="PyMOL executable name or local path")
    result.add_argument(
        "--emit-endpoint-template",
        type=Path,
        help="Write a pending endpoint-review template from PLIP and exit",
    )
    result.add_argument("--dry-run", action="store_true", help="Validate and print the render plan without writing or invoking PyMOL")
    result.add_argument("--overwrite", action="store_true", help="Replace only the exact declared output files")
    return result


def main(argv: Sequence[str] | None = None) -> int:
    args = parser().parse_args(argv)
    try:
        structure = require_file(args.structure, "coordinate file")
        report = require_file(args.plip_report, "PLIP report")
        normalized = parse_plip_report(report, args.binding_site)
        structure_hash = sha256_file(structure)

        if args.emit_endpoint_template is not None:
            if args.dry_run:
                raise ValidationError("--emit-endpoint-template and --dry-run cannot be combined")
            target = args.emit_endpoint_template.expanduser().resolve()
            if target.exists() and not args.overwrite:
                raise ValidationError(f"Endpoint template already exists; use --overwrite: {target}")
            target.parent.mkdir(parents=True, exist_ok=True)
            write_json(target, endpoint_template(normalized, structure_hash))
            print(f"Wrote pending endpoint-review template: {target}")
            print(f"PLIP interactions requiring review: {len(normalized['interactions'])}")
            return 0

        for name in ("provenance", "endpoint_map", "output_dir"):
            if getattr(args, name) is None:
                raise ValidationError(f"--{name.replace('_', '-')} is required for render preflight")
        if not SIMPLE_IDENTIFIER.fullmatch(args.output_stem):
            raise ValidationError("--output-stem may contain only letters, digits, period, underscore, plus, or hyphen")

        provenance_path = require_file(args.provenance, "pose provenance")
        endpoint_path = require_file(args.endpoint_map, "endpoint map")
        style_path = require_file(args.style, "figure style")
        style = validate_style(style_path)
        style_hash = sha256_file(style_path)
        provenance = validate_provenance(provenance_path, structure_hash, normalized)
        validate_endpoint_map(endpoint_path, normalized, structure_hash)
        camera_matrix = None
        camera_path = None
        if args.camera is not None:
            camera_path = require_file(args.camera, "camera record")
            camera_matrix = validate_camera_input(camera_path, structure_hash, style_hash)

        output_dir = args.output_dir.expanduser().resolve()
        outputs = output_paths(output_dir, args.output_stem)
        plan = {
            "status": "preflight_passed",
            "dry_run": bool(args.dry_run),
            "structure": str(structure),
            "plip_report": str(report),
            "plip_version": normalized["plip_version"],
            "binding_site": normalized["binding_site"],
            "pose_provenance": provenance["pose_provenance"],
            "structure_id": provenance["structure_id"],
            "ligand_id": provenance["ligand_id"],
            "ligand_formal_charge": provenance["ligand_formal_charge"],
            "visible_provenance_label": provenance["visible_label"],
            "interaction_count": len(normalized["interactions"]),
            "interaction_ids": [item["interaction_id"] for item in normalized["interactions"]],
            "camera_source": str(camera_path) if camera_path else "deterministic_orient_then_zoom",
            "output_dir": str(output_dir),
            "dpi": int(style["render"]["dpi"]),
        }
        if args.dry_run:
            print(json.dumps(plan, indent=2, sort_keys=True))
            return 0

        ensure_outputs_available(outputs, args.overwrite)
        output_dir.mkdir(parents=True, exist_ok=True)
        write_json(outputs["normalized_interactions"], normalized)
        pml_text = build_pml(
            structure=structure,
            normalized=normalized,
            provenance=provenance,
            style=style,
            style_hash=style_hash,
            structure_hash=structure_hash,
            pml_path=outputs["pml"],
            pse_path=outputs["pse"],
            png_path=outputs["png"],
            camera_output=outputs["camera"],
            camera_matrix=camera_matrix,
        )
        outputs["pml"].write_text(pml_text, encoding="utf-8", newline="\n")

        pymol_executable = resolve_executable(args.pymol_executable)
        command = [str(pymol_executable), "-cq", str(outputs["pml"])]
        completed = subprocess.run(
            command,
            cwd=output_dir,
            text=True,
            encoding="utf-8",
            errors="replace",
            capture_output=True,
            check=False,
        )
        outputs["stdout_log"].write_text(completed.stdout, encoding="utf-8")
        outputs["stderr_log"].write_text(completed.stderr, encoding="utf-8")
        if completed.returncode != 0:
            raise ValidationError(
                f"PyMOL exited with code {completed.returncode}; inspect {outputs['stderr_log']}"
            )
        for key in ("pse", "png", "camera"):
            if not outputs[key].is_file() or outputs[key].stat().st_size == 0:
                raise ValidationError(f"PyMOL completed without a non-empty {key} output: {outputs[key]}")

        manifest = build_manifest(
            structure=structure,
            report=report,
            provenance_path=provenance_path,
            endpoint_path=endpoint_path,
            style_path=style_path,
            normalized=normalized,
            provenance=provenance,
            binding_site=args.binding_site,
            outputs=outputs,
            pymol_executable=pymol_executable,
            command=command,
        )
        write_json(outputs["manifest"], manifest)
        print(f"Render completed: {outputs['png']}")
        print(f"Manifest: {outputs['manifest']}")
        return 0
    except ValidationError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    except OSError as exc:
        print(f"ERROR: operating-system failure: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
