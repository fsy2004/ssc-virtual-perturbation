#!/usr/bin/env python3
"""Credential-safe connector for the current user-designated AutoDL MD node."""

from __future__ import annotations

import re
from pathlib import Path

import paramiko


CREDENTIAL_FILE = Path(r"C:\Users\fsy\Desktop\密码.txt")
SSH_RE = re.compile(r"ssh\s+-p\s+(\d+)\s+([^@\s]+)@([^\s]+)", re.IGNORECASE)


def connection_material() -> tuple[str, int, str, str]:
    lines = CREDENTIAL_FILE.read_text(encoding="utf-8-sig").splitlines()
    for index, raw in enumerate(lines):
        match = SSH_RE.search(raw)
        if not match:
            continue
        port, username, host = int(match.group(1)), match.group(2), match.group(3)
        for candidate in lines[index + 1 :]:
            value = candidate.strip()
            if not value:
                continue
            value = re.sub(r"^(?:密码|password|passwd|pwd)\s*[:：=]\s*", "", value, flags=re.I)
            if not value or SSH_RE.search(value):
                break
            return host, port, username, value
        raise RuntimeError("No password line found after the current SSH endpoint")
    raise RuntimeError("No SSH endpoint found in the designated credential file")


def connect() -> paramiko.SSHClient:
    host, port, username, password = connection_material()
    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect(host, port=port, username=username, password=password, timeout=30)
    return client


def run(client: paramiko.SSHClient, command: str, timeout: int = 300) -> tuple[int, str, str]:
    _, stdout, stderr = client.exec_command(command, timeout=timeout)
    out = stdout.read().decode(errors="replace")
    err = stderr.read().decode(errors="replace")
    return stdout.channel.recv_exit_status(), out, err


def endpoint_label() -> str:
    host, port, _, _ = connection_material()
    return f"{host}:{port}"

