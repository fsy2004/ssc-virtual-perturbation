#!/usr/bin/env python3
"""Minimal, redacted CHARMM-GUI Quick Bilayer API client.

PDB Reader remains a manual scientific checkpoint. This client submits only an
already-approved PDB Reader job ID. Account login uses a hidden password prompt,
saves the bearer token outside the package, and never writes either secret to a
project record.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import ssl
import sys
import tarfile
import tempfile
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from getpass import getpass
from pathlib import Path, PurePosixPath
from typing import Any


DEFAULT_BASE = "https://www.charmm-gui.org/api"
LOCAL_TEST_TOKEN = "LOCALHOST_TEST_ONLY_NONSECRET"
JOB_RE = re.compile(r"^[A-Za-z0-9._-]+$")
ALLOWED_API_PATHS = {"/api/quick_bilayer", "/api/check_status", "/api/download"}
LOGIN_PATH = "/api/login"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def package_root() -> Path:
    return Path(__file__).resolve().parent.parent


def validate_base_url(value: str, allow_http_local: bool) -> str:
    normalized = value.rstrip("/")
    parsed = urllib.parse.urlparse(normalized)
    if normalized == DEFAULT_BASE and parsed.username is None and parsed.password is None:
        return normalized
    if (
        allow_http_local and parsed.scheme == "http" and parsed.hostname in {"127.0.0.1", "localhost"}
        and parsed.username is None and parsed.password is None and not parsed.query and not parsed.fragment
        and parsed.path.rstrip("/") == "/api"
    ):
        return normalized
    raise ValueError(
        f"bearer-authenticated requests are pinned to {DEFAULT_BASE}; custom bases are allowed only for explicit HTTP localhost tests"
    )


def token_for_base(args: argparse.Namespace, base: str) -> str:
    if base == DEFAULT_BASE:
        return read_token(args.token_file)
    if args.token_file or os.environ.get("CHARMMGUI_TOKEN_FILE"):
        raise ValueError("a real token file/environment token is prohibited with a localhost test endpoint")
    return LOCAL_TEST_TOKEN


def validate_request_destination(url: str, token: str) -> None:
    """Pin every bearer-bearing request, including calls made outside ``main``."""
    parsed = urllib.parse.urlparse(url)
    common_ok = (
        parsed.username is None and parsed.password is None and not parsed.fragment
        and parsed.path in ALLOWED_API_PATHS
    )
    official = (
        common_ok and parsed.scheme == "https" and parsed.netloc == "www.charmm-gui.org"
        and token != LOCAL_TEST_TOKEN
    )
    local_test = (
        common_ok and parsed.scheme == "http" and parsed.hostname in {"127.0.0.1", "localhost"}
        and token == LOCAL_TEST_TOKEN
    )
    if not (official or local_test):
        raise ValueError("refusing to send an Authorization header outside the pinned CHARMM-GUI API destination")


def read_token(path_arg: str | None) -> str:
    candidate = path_arg or os.environ.get("CHARMMGUI_TOKEN_FILE")
    if not candidate:
        raise ValueError("provide --token-file or set CHARMMGUI_TOKEN_FILE")
    path = Path(candidate).expanduser().resolve()
    root = package_root()
    if path == root or root in path.parents:
        raise ValueError("token file must be outside the project package")
    if not path.is_file():
        raise ValueError(f"token file does not exist: {path}")
    token = path.read_text(encoding="utf-8").strip()
    if not token or any(char.isspace() for char in token):
        raise ValueError("token is empty or contains whitespace")
    return token


def request(
    url: str,
    token: str,
    data: bytes | None = None,
    timeout: float = 60.0,
) -> urllib.response.addinfourl:
    validate_request_destination(url, token)
    headers = {
        "Authorization": f"Bearer {token}",
        "User-Agent": "ssc-gamma-secretase-rebuild/2.0",
        "Accept": "application/json, application/gzip, application/octet-stream",
    }
    if data is not None:
        headers["Content-Type"] = "application/x-www-form-urlencoded"
    req = urllib.request.Request(url, data=data, headers=headers, method="POST" if data is not None else "GET")

    class NoRedirectWithAuthorization(urllib.request.HTTPRedirectHandler):
        def redirect_request(self, req, fp, code, msg, headers, newurl):  # type: ignore[no-untyped-def]
            return None

    opener = urllib.request.build_opener(
        urllib.request.HTTPSHandler(context=ssl.create_default_context()),
        urllib.request.HTTPHandler(),
        NoRedirectWithAuthorization(),
    )
    return opener.open(req, timeout=timeout)


def login_request(base: str, email: str, password: str, timeout: float) -> dict[str, Any]:
    """Authenticate only to the pinned official login endpoint without logging secrets."""
    url = f"{base}/login"
    parsed = urllib.parse.urlparse(url)
    if not (
        parsed.scheme == "https"
        and parsed.netloc == "www.charmm-gui.org"
        and parsed.path == LOGIN_PATH
        and parsed.username is None
        and parsed.password is None
        and not parsed.query
        and not parsed.fragment
    ):
        raise ValueError("login is pinned to https://www.charmm-gui.org/api/login")
    payload = json.dumps({"email": email, "password": password}).encode("utf-8")
    req = urllib.request.Request(
        url,
        data=payload,
        headers={
            "Content-Type": "application/json",
            "Accept": "application/json",
            "User-Agent": "ssc-gamma-secretase-rebuild/2.0",
        },
        method="POST",
    )

    class NoRedirect(urllib.request.HTTPRedirectHandler):
        def redirect_request(self, req, fp, code, msg, headers, newurl):  # type: ignore[no-untyped-def]
            return None

    opener = urllib.request.build_opener(
        urllib.request.HTTPSHandler(context=ssl.create_default_context()), NoRedirect()
    )
    with opener.open(req, timeout=timeout) as response:
        return parse_json_response(response)


def save_token_outside_package(token: str, destination: Path) -> Path:
    root = package_root()
    destination = destination.expanduser().resolve()
    if destination == root or root in destination.parents:
        raise ValueError("token output must be outside the project package")
    if not token or any(char.isspace() for char in token):
        raise ValueError("login response token is empty or malformed")
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", dir=destination.parent, delete=False, prefix=".charmmgui-token."
    ) as handle:
        handle.write(token)
        temporary = Path(handle.name)
    os.chmod(temporary, 0o600)
    temporary.replace(destination)
    return destination


def record_json(directory: Path, stem: str, payload: dict[str, Any]) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    target = directory / f"{stem}.json"
    text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=directory, delete=False, prefix=f".{stem}.") as handle:
        handle.write(text)
        temporary = Path(handle.name)
    temporary.replace(target)
    return target


def parse_json_response(response: urllib.response.addinfourl) -> dict[str, Any]:
    raw = response.read()
    try:
        value = json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise RuntimeError("server response was not valid JSON") from exc
    if not isinstance(value, dict):
        raise RuntimeError("server response JSON was not an object")
    return value


def safe_tar_members(path: Path) -> tuple[int, list[str]]:
    names: list[str] = []
    with tarfile.open(path, mode="r:*") as archive:
        members = archive.getmembers()
        if not members:
            raise ValueError("downloaded archive is empty")
        for member in members:
            posix = PurePosixPath(member.name)
            if posix.is_absolute() or ".." in posix.parts:
                raise ValueError(f"unsafe archive path: {member.name}")
            if member.issym() or member.islnk():
                raise ValueError(f"archive links are not accepted: {member.name}")
            names.append(member.name)
    return len(names), names


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-url", default=DEFAULT_BASE)
    parser.add_argument("--allow-http-local-test", action="store_true")
    parser.add_argument("--token-file", help="path outside the project; alternatively set CHARMMGUI_TOKEN_FILE")
    parser.add_argument("--timeout", type=float, default=60.0)
    parser.add_argument("--record-dir", type=Path, default=package_root() / "api_records")
    sub = parser.add_subparsers(dest="command", required=True)

    login = sub.add_parser("login", help="obtain a bearer token using a hidden password prompt")
    login.add_argument("--email", required=True)
    login.add_argument("--output-token", required=True, type=Path)

    submit = sub.add_parser("submit", help="submit one Quick Bilayer job")
    submit.add_argument("--pdb-reader-jobid", required=True)
    submit.add_argument("--system-id", default="8kct_nirogacestat_native")
    submit.add_argument("--build-id", default="build01")
    submit.add_argument("--upper", default="POPC=1")
    submit.add_argument("--lower", default="POPC=1")
    submit.add_argument("--margin", type=float, required=True)
    submit.add_argument("--wdist", type=float, default=22.5)
    submit.add_argument("--ion-conc", type=float, default=0.15)
    submit.add_argument("--ion-type", default="NaCl")
    submit.add_argument("--ppm", action="store_true")
    submit.add_argument("--heteroatoms", action="store_true")
    submit.add_argument("--dry-run", action="store_true")

    status = sub.add_parser("status", help="check one job")
    status.add_argument("--jobid", required=True)

    download = sub.add_parser("download", help="download and safely validate a finished archive")
    download.add_argument("--jobid", required=True)
    download.add_argument("--output", required=True, type=Path)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    try:
        base = validate_base_url(args.base_url, args.allow_http_local_test)
        if args.command == "login":
            if base != DEFAULT_BASE:
                raise ValueError("real account login is allowed only at the pinned official API base")
            email = args.email.strip()
            if not re.fullmatch(r"[^@\s]+@[^@\s]+\.[^@\s]+", email):
                raise ValueError("invalid email address")
            password = getpass("CHARMM-GUI password (hidden): ")
            if not password:
                raise ValueError("password is empty")
            result = login_request(base, email, password, args.timeout)
            password = ""
            token = str(result.get("token", ""))
            destination = save_token_outside_package(token, args.output_token)
            print(
                json.dumps(
                    {
                        "authenticated": True,
                        "endpoint": f"{base}/login",
                        "token_file": str(destination),
                        "token_recorded_in_project": False,
                    },
                    indent=2,
                )
            )
            return 0
        if hasattr(args, "jobid") and not JOB_RE.fullmatch(args.jobid):
            raise ValueError("job ID contains unexpected characters")

        if args.command == "submit":
            if not JOB_RE.fullmatch(args.pdb_reader_jobid):
                raise ValueError("PDB Reader job ID contains unexpected characters")
            if args.system_id != "8kct_nirogacestat_native" or args.build_id != "build01":
                raise ValueError("this frozen package permits only system 8kct_nirogacestat_native and build01")
            if args.upper != "POPC=1" or args.lower != "POPC=1":
                raise ValueError("this frozen package permits only symmetric pure POPC leaflets (POPC=1)")
            if args.ion_type != "NaCl" or abs(args.ion_conc - 0.15) > 1e-12:
                raise ValueError("this frozen package permits only 0.15 M NaCl")
            if abs(args.margin - 20.0) > 1e-12 or abs(args.wdist - 22.5) > 1e-12:
                raise ValueError("this frozen package permits only margin=20 and wdist=22.5 Angstrom")
            if not args.ppm:
                raise ValueError("--ppm is mandatory; Quick Bilayer must apply one PPM orientation for post-build review")
            if not args.heteroatoms:
                raise ValueError("--heteroatoms is mandatory to retain native O6U, CLR, PC1, and glycans")
            if args.margin <= 0 or args.wdist <= 0 or args.ion_conc < 0:
                raise ValueError("margin/wdist must be positive and ion concentration non-negative")
            payload = {
                "jobid": args.pdb_reader_jobid,
                "upper": args.upper,
                "lower": args.lower,
                "margin": str(args.margin),
                "wdist": str(args.wdist),
                "Ion_conc": str(args.ion_conc),
                "Ion_type": args.ion_type,
                "clone_job": "false",
                "ppm": "true",
                "topologyIn": "true",
                "heteroatoms": "true",
            }
            audit_record: dict[str, Any] = {
                "recorded_at_utc": datetime.now(timezone.utc).isoformat(),
                "endpoint": f"{base}/quick_bilayer",
                "system_id": args.system_id,
                "build_id": args.build_id,
                "request": payload,
                "authentication": "bearer token supplied externally; not recorded",
                "client_sha256": sha256(Path(__file__).resolve()),
                "dry_run": args.dry_run,
            }
            if args.dry_run:
                print(json.dumps(audit_record, indent=2, sort_keys=True))
                return 0
            token = token_for_base(args, base)
            encoded = urllib.parse.urlencode(payload).encode("ascii")
            with request(f"{base}/quick_bilayer", token, encoded, args.timeout) as response:
                result = parse_json_response(response)
            jobid = str(result.get("jobid", ""))
            if not jobid or not JOB_RE.fullmatch(jobid):
                raise RuntimeError("submission response did not contain a valid job ID")
            if str(result.get("submitted", "")).lower() not in {"true", "1"}:
                raise RuntimeError(f"server did not confirm submission: {result}")
            audit_record["response"] = result
            audit_record["quick_bilayer_jobid"] = jobid
            path = record_json(args.record_dir.resolve(), f"{args.system_id}_{args.build_id}_{jobid}_submit", audit_record)
            print(json.dumps({"jobid": jobid, "record": str(path)}, indent=2))
            return 0

        token = token_for_base(args, base)
        if args.command == "status":
            query = urllib.parse.urlencode({"jobid": args.jobid})
            with request(f"{base}/check_status?{query}", token, timeout=args.timeout) as response:
                result = parse_json_response(response)
            path = record_json(
                args.record_dir.resolve(),
                f"{args.jobid}_status",
                {"recorded_at_utc": datetime.now(timezone.utc).isoformat(), "jobid": args.jobid, "response": result},
            )
            print(json.dumps({"status": result, "record": str(path)}, indent=2, sort_keys=True))
            return 0

        output = args.output.resolve()
        root = package_root()
        if output == root or root not in output.parents:
            raise ValueError("download output must be inside this package")
        if output.exists():
            raise ValueError(f"refusing to overwrite existing output: {output}")
        output.parent.mkdir(parents=True, exist_ok=True)
        query = urllib.parse.urlencode({"jobid": args.jobid})
        partial = output.with_suffix(output.suffix + ".part")
        if partial.exists():
            partial.unlink()
        try:
            with request(f"{base}/download?{query}", token, timeout=args.timeout) as response, partial.open("xb") as handle:
                content_type = response.headers.get("Content-Type", "").lower()
                if "json" in content_type or "text/html" in content_type or "text/plain" in content_type:
                    body = response.read(8192).decode("utf-8", errors="replace")
                    raise RuntimeError(f"server did not return an archive: {body[:500]}")
                while True:
                    block = response.read(1024 * 1024)
                    if not block:
                        break
                    handle.write(block)
            if partial.stat().st_size == 0:
                raise RuntimeError("downloaded archive is empty")
            member_count, names = safe_tar_members(partial)
            partial.replace(output)
        finally:
            if partial.exists():
                partial.unlink()
        digest = sha256(output)
        record = {
            "recorded_at_utc": datetime.now(timezone.utc).isoformat(),
            "jobid": args.jobid,
            "output": str(output.relative_to(root)),
            "bytes": output.stat().st_size,
            "sha256": digest,
            "safe_archive_members": member_count,
            "sample_members": names[:20],
        }
        path = record_json(args.record_dir.resolve(), f"{args.jobid}_download", record)
        print(json.dumps({"archive": str(output), "sha256": digest, "record": str(path)}, indent=2))
        return 0
    except (ValueError, RuntimeError, OSError, tarfile.TarError, urllib.error.URLError, urllib.error.HTTPError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
