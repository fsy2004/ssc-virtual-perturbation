#!/usr/bin/env python3
"""Validate an audited PLIP-to-PyMOL figure bundle without modifying it.

The validator independently reparses the archived PLIP report, rechecks the
approved endpoint and pose-provenance records, verifies every declared file
hash, and confirms that the PML displays exactly the PLIP interaction set.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import struct
import sys
import zlib
from collections import Counter
from pathlib import Path
from typing import Any, Sequence

from render_plip_pymol import (
    ValidationError,
    build_pml,
    interaction_display_segments,
    load_json_object,
    parse_plip_report,
    render_object_name,
    require_file,
    require_finite_number,
    require_nonempty_string,
    sanitize_object_name,
    sha256_file,
    stable_json,
    validate_endpoint_map,
    validate_provenance,
    validate_style,
)


EXPECTED_INPUT_KEYS = {
    "structure",
    "plip_report",
    "provenance",
    "endpoint_map",
    "style",
}

EXPECTED_OUTPUT_KEYS = {
    "normalized_interactions",
    "pml",
    "pse",
    "png",
    "camera",
    "stdout_log",
    "stderr_log",
}

SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")
DISTANCE_CALL_PATTERN = re.compile(r"cmd\.distance\(\s*(['\"])([^'\"]+)\1")
FATAL_LOG_PATTERNS = (
    re.compile(r"Traceback \(most recent call last\)", re.IGNORECASE),
    re.compile(r"Segmentation fault", re.IGNORECASE),
    re.compile(r"Fatal Error", re.IGNORECASE),
    re.compile(r"Unhandled exception", re.IGNORECASE),
)


def require_mapping(value: Any, label: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ValidationError(f"{label} must be a JSON object")
    return value


def manifest_file_record(
    container: dict[str, Any],
    key: str,
    label: str,
    require_size: bool,
    allow_empty: bool = False,
) -> tuple[Path, dict[str, Any]]:
    record = require_mapping(container.get(key), f"{label}.{key}")
    raw_path = require_nonempty_string(record.get("path"), f"{label}.{key}.path")
    path = Path(raw_path).expanduser()
    if not path.is_absolute():
        raise ValidationError(f"{label}.{key}.path must be absolute")
    path = path.resolve()
    if not path.is_file():
        raise ValidationError(f"{label}.{key} is missing: {path}")
    if not allow_empty and path.stat().st_size == 0:
        raise ValidationError(f"{label}.{key} is empty: {path}")
    recorded_hash = str(record.get("sha256", "")).lower()
    if not SHA256_PATTERN.fullmatch(recorded_hash):
        raise ValidationError(f"{label}.{key}.sha256 must be a lowercase SHA-256 digest")
    if require_size:
        try:
            recorded_bytes = int(record.get("bytes"))
        except (TypeError, ValueError) as exc:
            raise ValidationError(f"{label}.{key}.bytes must be an integer") from exc
        minimum = 0 if allow_empty else 1
        if recorded_bytes < minimum:
            qualifier = "non-negative" if allow_empty else "positive"
            raise ValidationError(f"{label}.{key}.bytes must be {qualifier}")
    return path, record


def validate_declared_hash(path: Path, record: dict[str, Any], label: str) -> None:
    observed_hash = sha256_file(path)
    expected_hash = str(record["sha256"]).lower()
    if observed_hash != expected_hash:
        raise ValidationError(f"{label} SHA-256 differs from the render manifest")
    if "bytes" in record and path.stat().st_size != int(record["bytes"]):
        raise ValidationError(f"{label} byte count differs from the render manifest")


def validate_camera_record(
    path: Path,
    manifest_render: dict[str, Any],
    structure_hash: str,
    style_hash: str,
) -> dict[str, Any]:
    data = load_json_object(path, "camera output")
    if data.get("schema_version") != "1.0":
        raise ValidationError("camera output schema_version must be 1.0")
    if data.get("camera_source") not in {"provided", "deterministic_orient_then_zoom"}:
        raise ValidationError("camera output camera_source is unsupported")
    values = data.get("view_matrix")
    if not isinstance(values, list) or len(values) != 18:
        raise ValidationError("camera output view_matrix must contain exactly 18 values")
    matrix = [require_finite_number(value, "camera output view_matrix value") for value in values]
    rotation = [matrix[0:3], matrix[3:6], matrix[6:9]]
    for index, row in enumerate(rotation, start=1):
        norm = math.sqrt(sum(value * value for value in row))
        if not math.isclose(norm, 1.0, rel_tol=0.0, abs_tol=0.02):
            raise ValidationError(f"camera rotation row {index} is not unit length")
    for left, right in ((0, 1), (0, 2), (1, 2)):
        dot_product = sum(a * b for a, b in zip(rotation[left], rotation[right]))
        if not math.isclose(dot_product, 0.0, rel_tol=0.0, abs_tol=0.02):
            raise ValidationError("camera rotation rows are not mutually orthogonal")
    determinant = (
        rotation[0][0] * (rotation[1][1] * rotation[2][2] - rotation[1][2] * rotation[2][1])
        - rotation[0][1] * (rotation[1][0] * rotation[2][2] - rotation[1][2] * rotation[2][0])
        + rotation[0][2] * (rotation[1][0] * rotation[2][1] - rotation[1][1] * rotation[2][0])
    )
    if not math.isclose(abs(determinant), 1.0, rel_tol=0.0, abs_tol=0.03):
        raise ValidationError("camera rotation matrix determinant is not approximately +/-1")
    pymol_version = data.get("pymol_version")
    if not isinstance(pymol_version, list) or not pymol_version:
        raise ValidationError("camera output must record a non-empty PyMOL version sequence")
    if str(data.get("structure_sha256", "")).lower() != structure_hash:
        raise ValidationError("camera output structure SHA-256 differs from the manifest input")
    if str(data.get("style_sha256", "")).lower() != style_hash:
        raise ValidationError("camera output style SHA-256 differs from the manifest input")
    for key in ("width_px", "height_px", "dpi"):
        try:
            observed = int(data.get(key))
            expected = int(manifest_render.get(key))
        except (TypeError, ValueError) as exc:
            raise ValidationError(f"camera/manifest {key} must be an integer") from exc
        if observed != expected:
            raise ValidationError(f"camera output {key} differs from the render manifest")
    return data


def read_png_metadata(path: Path) -> dict[str, Any]:
    """Decode the PNG stream and verify dimensions, resolution, and chunk CRCs."""

    signature = b"\x89PNG\r\n\x1a\n"
    width: int | None = None
    height: int | None = None
    bit_depth: int | None = None
    color_type: int | None = None
    channels: int | None = None
    scanline_bytes: int | None = None
    physical: tuple[int, int, int] | None = None
    seen_ihdr = 0
    seen_phys = 0
    seen_idat = 0
    seen_iend = 0
    idat_ended = False
    decompressor: zlib.Decompress | None = None
    decoded_bytes = 0

    def consume_decoded(payload: bytes) -> None:
        nonlocal decoded_bytes
        if scanline_bytes is None or height is None:
            raise ValidationError(f"PNG IDAT appeared before a valid IHDR: {path}")
        expected_bytes = scanline_bytes * height
        if decoded_bytes + len(payload) > expected_bytes:
            raise ValidationError(f"PNG decompressed data exceed its declared dimensions: {path}")
        start = 0
        while start < len(payload):
            row_offset = decoded_bytes % scanline_bytes
            take = min(len(payload) - start, scanline_bytes - row_offset)
            if row_offset == 0 and payload[start] > 4:
                raise ValidationError(f"PNG contains an invalid scanline filter byte: {path}")
            decoded_bytes += take
            start += take

    with path.open("rb") as handle:
        if handle.read(8) != signature:
            raise ValidationError(f"PNG signature is invalid: {path}")
        chunk_index = 0
        while True:
            header = handle.read(8)
            if len(header) != 8:
                raise ValidationError(f"PNG ended before a valid IEND chunk: {path}")
            length, chunk_type = struct.unpack(">I4s", header)
            chunk_index += 1
            if chunk_index == 1 and chunk_type != b"IHDR":
                raise ValidationError(f"PNG first chunk must be IHDR: {path}")
            if length > 2**31:
                raise ValidationError(f"PNG contains an implausibly large chunk: {path}")
            payload = handle.read(length)
            checksum_raw = handle.read(4)
            if len(payload) != length or len(checksum_raw) != 4:
                raise ValidationError(f"PNG chunk is truncated: {path}")
            expected_crc = struct.unpack(">I", checksum_raw)[0]
            observed_crc = zlib.crc32(chunk_type)
            observed_crc = zlib.crc32(payload, observed_crc) & 0xFFFFFFFF
            if observed_crc != expected_crc:
                label = chunk_type.decode("ascii", errors="replace")
                raise ValidationError(f"PNG chunk CRC failed for {label}: {path}")

            if chunk_type == b"IHDR":
                seen_ihdr += 1
                if seen_ihdr != 1 or length != 13:
                    raise ValidationError(f"PNG must contain one valid IHDR chunk: {path}")
                (
                    width,
                    height,
                    bit_depth,
                    color_type,
                    compression_method,
                    filter_method,
                    interlace_method,
                ) = struct.unpack(">IIBBBBB", payload)
                if width <= 0 or height <= 0:
                    raise ValidationError(f"PNG dimensions must be positive: {path}")
                if bit_depth != 8 or color_type not in {2, 6}:
                    raise ValidationError(
                        f"PNG master must be 8-bit RGB or RGBA, got bit_depth={bit_depth}, "
                        f"color_type={color_type}: {path}"
                    )
                if compression_method != 0 or filter_method != 0 or interlace_method != 0:
                    raise ValidationError(
                        f"PNG must use standard compression/filtering and be non-interlaced: {path}"
                    )
                channels = 3 if color_type == 2 else 4
                scanline_bytes = 1 + width * channels
            elif chunk_type == b"pHYs":
                seen_phys += 1
                if seen_phys != 1 or length != 9:
                    raise ValidationError(f"PNG must contain one valid pHYs chunk: {path}")
                if seen_idat:
                    raise ValidationError(f"PNG pHYs metadata must precede image data: {path}")
                physical = struct.unpack(">IIB", payload)
            elif chunk_type == b"IDAT":
                if seen_ihdr != 1 or idat_ended:
                    raise ValidationError(f"PNG IDAT chunks are out of order: {path}")
                seen_idat += 1
                if decompressor is None:
                    decompressor = zlib.decompressobj()
                compressed = payload
                try:
                    while compressed:
                        decoded = decompressor.decompress(compressed, 1024 * 1024)
                        consume_decoded(decoded)
                        compressed = decompressor.unconsumed_tail
                        if not decoded and not compressed:
                            break
                except zlib.error as exc:
                    raise ValidationError(f"PNG IDAT zlib stream is invalid: {path}: {exc}") from exc
            elif chunk_type == b"IEND":
                seen_iend += 1
                if length != 0 or seen_iend != 1:
                    raise ValidationError(f"PNG must contain one valid IEND chunk: {path}")
                if decompressor is None:
                    raise ValidationError(f"PNG lacks a decodable IDAT stream: {path}")
                try:
                    consume_decoded(decompressor.flush())
                except zlib.error as exc:
                    raise ValidationError(f"PNG IDAT zlib stream cannot be finalized: {path}: {exc}") from exc
                if not decompressor.eof or decompressor.unused_data:
                    raise ValidationError(f"PNG IDAT zlib stream is incomplete or contains extras: {path}")
                if scanline_bytes is None or height is None or decoded_bytes != scanline_bytes * height:
                    raise ValidationError(f"PNG decoded byte count differs from its declared dimensions: {path}")
                if handle.read(1):
                    raise ValidationError(f"PNG contains trailing bytes after IEND: {path}")
                break
            elif seen_idat:
                idat_ended = True

    if seen_ihdr != 1 or seen_idat == 0 or seen_iend != 1:
        raise ValidationError(f"PNG lacks a required IHDR, IDAT, or IEND chunk: {path}")
    if physical is None:
        raise ValidationError(f"PNG lacks embedded pHYs resolution metadata: {path}")
    x_pixels_per_metre, y_pixels_per_metre, unit = physical
    if unit != 1:
        raise ValidationError(f"PNG pHYs resolution unit is not metres: {path}")
    return {
        "width_px": width,
        "height_px": height,
        "bit_depth": bit_depth,
        "color_type": color_type,
        "channels": channels,
        "dpi_x": x_pixels_per_metre * 0.0254,
        "dpi_y": y_pixels_per_metre * 0.0254,
        "x_pixels_per_metre": x_pixels_per_metre,
        "y_pixels_per_metre": y_pixels_per_metre,
    }


def validate_png(
    path: Path,
    manifest_render: dict[str, Any],
) -> dict[str, Any]:
    metadata = read_png_metadata(path)
    for key in ("width_px", "height_px"):
        if int(metadata[key]) != int(manifest_render[key]):
            raise ValidationError(f"PNG {key} differs from the render manifest")
    expected_dpi = int(manifest_render["dpi"])
    for axis in ("dpi_x", "dpi_y"):
        observed = float(metadata[axis])
        if not math.isclose(observed, expected_dpi, rel_tol=0.0, abs_tol=0.5):
            raise ValidationError(
                f"PNG embedded {axis} is {observed:.3f}, expected {expected_dpi} +/- 0.5"
            )
    return metadata


def require_pml_fragment(pml_text: str, fragment: str, label: str) -> None:
    if fragment not in pml_text:
        raise ValidationError(f"PML is missing the required {label}")


def validate_pml(
    path: Path,
    expected_pml: str,
    normalized: dict[str, Any],
    provenance: dict[str, Any],
    style: dict[str, Any],
    manifest_render: dict[str, Any],
    output_paths: dict[str, Path],
) -> dict[str, Any]:
    try:
        pml_text = path.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError) as exc:
        raise ValidationError(f"PML is not readable UTF-8: {path}: {exc}") from exc
    if pml_text != expected_pml:
        raise ValidationError(
            "PML differs from the deterministic script reconstructed from the approved inputs"
        )

    observed_names = [match.group(2) for match in DISTANCE_CALL_PATTERN.finditer(pml_text)]
    expected_names = [
        render_object_name(interaction["interaction_id"], segment["segment_id"])
        for interaction in normalized["interactions"]
        for segment in interaction_display_segments(interaction)
    ]
    if Counter(observed_names) != Counter(expected_names) or len(observed_names) != len(expected_names):
        missing = sorted((Counter(expected_names) - Counter(observed_names)).elements())
        extras = sorted((Counter(observed_names) - Counter(expected_names)).elements())
        raise ValidationError(
            f"PML contact objects differ from PLIP; missing={missing}, extras={extras}"
        )

    for interaction in normalized["interactions"]:
        interaction_id = interaction["interaction_id"]
        endpoint_names = {
            "protein": sanitize_object_name(f"endpoint_protein_{interaction_id}"),
            "ligand": sanitize_object_name(f"endpoint_ligand_{interaction_id}"),
        }
        endpoint_coordinates = {
            "protein": interaction["protein"]["coordinate"],
            "ligand": interaction["ligand"]["coordinate"],
        }
        if interaction["type"] == "water_bridge":
            endpoint_names["water"] = sanitize_object_name(f"endpoint_water_{interaction_id}")
            endpoint_coordinates["water"] = interaction["mediators"][0]["coordinate"]
        settings = style["interaction_styles"][interaction["type"]]
        color_name = f"interaction_{interaction['type']}"
        for role in sorted(endpoint_names):
            require_pml_fragment(
                pml_text,
                f"cmd.pseudoatom({endpoint_names[role]!r}, pos={endpoint_coordinates[role]!r})",
                f"{role} endpoint for {interaction_id}",
            )
        for segment in interaction_display_segments(interaction):
            object_name = render_object_name(interaction_id, segment["segment_id"])
            source_name = endpoint_names[segment["source_role"]]
            target_name = endpoint_names[segment["target_role"]]
            fragments = (
                (
                    f"cmd.distance({object_name!r}, {source_name!r}, {target_name!r})",
                    f"distance mapping for {interaction_id}:{segment['segment_id']}",
                ),
                (
                    f"cmd.set('dash_color', {color_name!r}, {object_name!r})",
                    f"color encoding for {interaction_id}:{segment['segment_id']}",
                ),
                (
                    f"cmd.set('dash_gap', {float(settings['dash_gap'])!r}, {object_name!r})",
                    f"dash-gap encoding for {interaction_id}:{segment['segment_id']}",
                ),
                (
                    f"cmd.set('dash_length', {float(settings['dash_length'])!r}, {object_name!r})",
                    f"dash-length encoding for {interaction_id}:{segment['segment_id']}",
                ),
                (
                    f"cmd.set('dash_width', {float(settings['dash_width'])!r}, {object_name!r})",
                    f"dash-width encoding for {interaction_id}:{segment['segment_id']}",
                ),
            )
            for fragment, label in fragments:
                require_pml_fragment(pml_text, fragment, label)

    width_px = int(manifest_render["width_px"])
    height_px = int(manifest_render["height_px"])
    dpi = int(manifest_render["dpi"])
    require_pml_fragment(pml_text, "cmd.reinitialize()", "deterministic reset")
    require_pml_fragment(
        pml_text,
        f"cmd.viewport({width_px}, {height_px})",
        "declared viewport",
    )
    require_pml_fragment(
        pml_text,
        f"cmd.save({output_paths['pse'].as_posix()!r})",
        "declared PSE destination",
    )
    require_pml_fragment(
        pml_text,
        (
            f"cmd.png({output_paths['png'].as_posix()!r}, width={width_px}, "
            f"height={height_px}, dpi={dpi}, ray=1, quiet=0)"
        ),
        "600-dpi ray-traced PNG command",
    )
    require_pml_fragment(
        pml_text,
        f"Path({output_paths['camera'].as_posix()!r}).write_text",
        "declared camera destination",
    )
    require_pml_fragment(
        pml_text,
        f"label={provenance['visible_label']!r}",
        "visible pose-provenance label",
    )
    return {
        "contact_object_count": len(observed_names),
        "contact_object_names": observed_names,
    }


def validate_logs(stdout_path: Path, stderr_path: Path) -> dict[str, Any]:
    logs: dict[str, str] = {}
    for key, path in (("stdout", stdout_path), ("stderr", stderr_path)):
        try:
            logs[key] = path.read_text(encoding="utf-8", errors="strict")
        except (OSError, UnicodeDecodeError) as exc:
            raise ValidationError(f"PyMOL {key} log is not readable UTF-8: {path}: {exc}") from exc
    combined = "\n".join(logs.values())
    for pattern in FATAL_LOG_PATTERNS:
        if pattern.search(combined):
            raise ValidationError(f"PyMOL logs contain a fatal marker: {pattern.pattern}")
    return {"stdout_bytes": len(logs["stdout"].encode("utf-8")), "stderr_bytes": len(logs["stderr"].encode("utf-8"))}


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        description=(
            "Validate a PLIP/PyMOL render bundle against its manifest, archived "
            "PLIP report, approved endpoints, provenance, and 600-dpi PNG metadata."
        )
    )
    result.add_argument("--manifest", required=True, type=Path, help="Render manifest produced by render_plip_pymol.py")
    result.add_argument(
        "--dry-run",
        action="store_true",
        help="Validate inputs and list output checks without hashing or decoding render outputs",
    )
    result.add_argument("--json", action="store_true", help="Print the result as JSON")
    return result


def execute_validation(manifest_path: Path, dry_run: bool) -> dict[str, Any]:
    manifest_path = require_file(manifest_path, "render manifest")
    manifest = load_json_object(manifest_path, "render manifest")
    if manifest.get("schema_version") != "figure-render-manifest-1.0":
        raise ValidationError("render manifest schema_version must be figure-render-manifest-1.0")

    inputs = require_mapping(manifest.get("inputs"), "manifest.inputs")
    outputs = require_mapping(manifest.get("outputs"), "manifest.outputs")
    if set(inputs) != EXPECTED_INPUT_KEYS:
        raise ValidationError(
            f"manifest input set must be exactly {sorted(EXPECTED_INPUT_KEYS)}, got {sorted(inputs)}"
        )
    if set(outputs) != EXPECTED_OUTPUT_KEYS:
        raise ValidationError(
            f"manifest output set must be exactly {sorted(EXPECTED_OUTPUT_KEYS)}, got {sorted(outputs)}"
        )

    input_paths: dict[str, Path] = {}
    input_records: dict[str, dict[str, Any]] = {}
    for key in sorted(EXPECTED_INPUT_KEYS):
        input_paths[key], input_records[key] = manifest_file_record(inputs, key, "inputs", False)
        validate_declared_hash(input_paths[key], input_records[key], f"input {key}")

    output_paths: dict[str, Path] = {}
    output_records: dict[str, dict[str, Any]] = {}
    for key in sorted(EXPECTED_OUTPUT_KEYS):
        output_paths[key], output_records[key] = manifest_file_record(
            outputs,
            key,
            "outputs",
            True,
            allow_empty=key in {"stdout_log", "stderr_log"},
        )

    binding_site_request = require_nonempty_string(
        manifest.get("binding_site_request"), "manifest binding_site_request"
    )
    normalized = parse_plip_report(input_paths["plip_report"], binding_site_request)
    structure_hash = sha256_file(input_paths["structure"])
    style = validate_style(input_paths["style"])
    provenance = validate_provenance(input_paths["provenance"], structure_hash, normalized)
    validate_endpoint_map(input_paths["endpoint_map"], normalized, structure_hash)

    if stable_json(manifest.get("binding_site")) != stable_json(normalized["binding_site"]):
        raise ValidationError("manifest binding-site record differs from the PLIP report")
    if str(manifest.get("plip_version", "")) != normalized["plip_version"]:
        raise ValidationError("manifest PLIP version differs from the PLIP report")
    if manifest.get("pose_provenance") != provenance["pose_provenance"]:
        raise ValidationError("manifest pose provenance differs from the approved provenance record")
    if manifest.get("structure_id") != provenance["structure_id"]:
        raise ValidationError("manifest structure_id differs from the approved provenance record")
    if manifest.get("ligand_id") != provenance["ligand_id"]:
        raise ValidationError("manifest ligand_id differs from the approved provenance record")
    if manifest.get("ligand_formal_charge") != provenance["ligand_formal_charge"]:
        raise ValidationError("manifest ligand formal charge differs from the approved provenance record")
    if manifest.get("visible_provenance_label") != provenance["visible_label"]:
        raise ValidationError("manifest visible provenance label differs from the approved record")

    expected_ids = [item["interaction_id"] for item in normalized["interactions"]]
    if manifest.get("interaction_ids") != expected_ids:
        raise ValidationError("manifest interaction IDs differ from the reparsed PLIP report")
    try:
        interaction_count = int(manifest.get("interaction_count"))
    except (TypeError, ValueError) as exc:
        raise ValidationError("manifest interaction_count must be an integer") from exc
    if interaction_count != len(expected_ids):
        raise ValidationError("manifest interaction_count differs from the reparsed PLIP report")

    manifest_render = require_mapping(manifest.get("render"), "manifest.render")
    for key in ("dpi", "width_px", "height_px"):
        try:
            if int(manifest_render.get(key)) <= 0:
                raise ValueError
        except (TypeError, ValueError) as exc:
            raise ValidationError(f"manifest.render.{key} must be a positive integer") from exc
    expected_width = int(round(float(style["render"]["width_cm"]) / 2.54 * int(style["render"]["dpi"])))
    expected_height = int(round(float(style["render"]["height_cm"]) / 2.54 * int(style["render"]["dpi"])))
    if int(manifest_render["dpi"]) != 600:
        raise ValidationError("manifest render must be exactly 600 dpi")
    if int(manifest_render["width_px"]) != expected_width or int(manifest_render["height_px"]) != expected_height:
        raise ValidationError("manifest pixel dimensions differ from the approved style")

    base_result: dict[str, Any] = {
        "manifest": str(manifest_path),
        "binding_site": normalized["binding_site"],
        "plip_version": normalized["plip_version"],
        "pose_provenance": provenance["pose_provenance"],
        "structure_id": provenance["structure_id"],
        "ligand_id": provenance["ligand_id"],
        "ligand_formal_charge": provenance["ligand_formal_charge"],
        "visible_provenance_label": provenance["visible_label"],
        "interaction_count": len(expected_ids),
        "interaction_ids": expected_ids,
        "render": {
            "dpi": int(manifest_render["dpi"]),
            "width_px": int(manifest_render["width_px"]),
            "height_px": int(manifest_render["height_px"]),
        },
    }
    if dry_run:
        base_result.update(
            {
                "status": "validation_plan_ready",
                "dry_run": True,
                "verified_now": [
                    "manifest_schema",
                    "input_hashes",
                    "plip_reparse",
                    "approved_endpoint_map",
                    "approved_pose_provenance",
                    "okabe_ito_dual_encoding_style",
                    "output_presence",
                ],
                "deferred_until_full_validation": [
                    "output_hashes_and_sizes",
                    "normalized_interaction_equality",
                    "pml_contact_and_provenance_scan",
                    "camera_record",
                    "png_crc_dimensions_and_embedded_dpi",
                    "pymol_fatal_log_scan",
                ],
                "outputs": {key: str(output_paths[key]) for key in sorted(output_paths)},
            }
        )
        return base_result

    for key in sorted(EXPECTED_OUTPUT_KEYS):
        validate_declared_hash(output_paths[key], output_records[key], f"output {key}")

    saved_normalized = load_json_object(
        output_paths["normalized_interactions"], "normalized interactions output"
    )
    if stable_json(saved_normalized) != stable_json(normalized):
        raise ValidationError("normalized interactions output differs from the reparsed PLIP report")

    style_hash = sha256_file(input_paths["style"])
    camera_result = validate_camera_record(
        output_paths["camera"],
        manifest_render,
        structure_hash,
        style_hash,
    )
    camera_matrix = (
        [float(value) for value in camera_result["view_matrix"]]
        if camera_result["camera_source"] == "provided"
        else None
    )
    expected_pml = build_pml(
        structure=input_paths["structure"],
        normalized=normalized,
        provenance=provenance,
        style=style,
        style_hash=style_hash,
        structure_hash=structure_hash,
        pml_path=output_paths["pml"],
        pse_path=output_paths["pse"],
        png_path=output_paths["png"],
        camera_output=output_paths["camera"],
        camera_matrix=camera_matrix,
    )
    pml_result = validate_pml(
        output_paths["pml"],
        expected_pml,
        normalized,
        provenance,
        style,
        manifest_render,
        output_paths,
    )
    png_result = validate_png(output_paths["png"], manifest_render)
    log_result = validate_logs(output_paths["stdout_log"], output_paths["stderr_log"])

    pymol = require_mapping(manifest.get("pymol"), "manifest.pymol")
    require_nonempty_string(pymol.get("executable"), "manifest.pymol.executable")
    if pymol.get("returncode") != 0:
        raise ValidationError("manifest PyMOL returncode must be exactly zero")
    command = pymol.get("command")
    if not isinstance(command, list) or not command or not all(isinstance(item, str) and item for item in command):
        raise ValidationError("manifest.pymol.command must be a non-empty string list")
    if str(output_paths["pml"]) not in command:
        raise ValidationError("manifest PyMOL command does not reference the declared PML output")
    if "-cq" not in command:
        raise ValidationError("manifest PyMOL command must record quiet headless execution with -cq")

    base_result.update(
        {
            "status": "validation_passed",
            "dry_run": False,
            "pml": pml_result,
            "camera_source": camera_result["camera_source"],
            "png": png_result,
            "logs": log_result,
            "verified_output_hashes": sorted(EXPECTED_OUTPUT_KEYS),
        }
    )
    return base_result


def main(argv: Sequence[str] | None = None) -> int:
    args = parser().parse_args(argv)
    try:
        result = execute_validation(args.manifest, args.dry_run)
        if args.json:
            print(json.dumps(result, indent=2, sort_keys=True, ensure_ascii=True))
        elif args.dry_run:
            print("Figure validation preflight passed.")
            print(f"Manifest: {result['manifest']}")
            print(f"PLIP interactions: {result['interaction_count']}")
            print("Render-output hashing and decoding are deferred until full validation.")
        else:
            print("Figure output validation passed.")
            print(f"Manifest: {result['manifest']}")
            print(f"PLIP interactions rendered: {result['interaction_count']}")
            print(
                "PNG: "
                f"{result['png']['width_px']} x {result['png']['height_px']} px at "
                f"{result['png']['dpi_x']:.3f} x {result['png']['dpi_y']:.3f} dpi"
            )
        return 0
    except ValidationError as exc:
        if args.json:
            print(json.dumps({"status": "validation_failed", "error": str(exc)}, ensure_ascii=True), file=sys.stderr)
        else:
            print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    except OSError as exc:
        message = f"operating-system failure: {exc}"
        if args.json:
            print(json.dumps({"status": "validation_failed", "error": message}, ensure_ascii=True), file=sys.stderr)
        else:
            print(f"ERROR: {message}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
