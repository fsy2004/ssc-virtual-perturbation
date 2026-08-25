from __future__ import annotations

import hashlib
import importlib.util
import io
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "_sync_analysis_to_server.py"


def load_module():
    spec = importlib.util.spec_from_file_location("sync_analysis_under_test", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


class FakeSftp:
    def __init__(self, payload: bytes):
        self.payload = payload

    def open(self, path: str, mode: str):
        self.last_open = (path, mode)
        return io.BytesIO(self.payload)


class SyncAnalysisTests(unittest.TestCase):
    def test_remote_sha256_reads_uploaded_bytes(self):
        module = load_module()
        payload = b"hash-bound-control-plane\n"
        sftp = FakeSftp(payload)
        observed = module.sha256_remote_file(sftp, "/release/scripts/example.py")
        self.assertEqual(observed, hashlib.sha256(payload).hexdigest())
        self.assertEqual(sftp.last_open, ("/release/scripts/example.py", "rb"))


if __name__ == "__main__":
    unittest.main()
