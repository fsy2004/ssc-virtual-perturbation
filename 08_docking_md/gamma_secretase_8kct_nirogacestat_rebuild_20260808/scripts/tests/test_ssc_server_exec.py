from __future__ import annotations

import contextlib
import importlib.util
import io
import tempfile
import unittest
from pathlib import Path
from unittest import mock


SCRIPT = Path(__file__).resolve().parents[1] / "_ssc_server_exec.py"


def load_module():
    spec = importlib.util.spec_from_file_location("ssc_server_exec_under_test", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


class SscServerExecTests(unittest.TestCase):
    def test_returns_remote_exit_code_and_closes_connection(self):
        module = load_module()
        connection = mock.Mock()
        with mock.patch.object(module, "connect", return_value=connection), mock.patch.object(
            module, "run", return_value=(17, "remote-out\n", "remote-err\n")
        ):
            stdout = io.StringIO()
            with contextlib.redirect_stdout(stdout):
                code = module.main(["printf ok"])
        self.assertEqual(code, 17)
        self.assertIn("remote-out", stdout.getvalue())
        self.assertIn("remote-err", stdout.getvalue())
        connection.close.assert_called_once_with()

    def test_at_file_reads_command_verbatim(self):
        module = load_module()
        connection = mock.Mock()
        with tempfile.TemporaryDirectory() as tmpdir:
            command_file = Path(tmpdir) / "command.sh"
            command_file.write_text("set -euo pipefail\nprintf 'safe'\n", encoding="utf-8")
            with mock.patch.object(module, "connect", return_value=connection), mock.patch.object(
                module, "run", return_value=(0, "", "")
            ) as remote_run:
                code = module.main([f"@{command_file}"])
        self.assertEqual(code, 0)
        self.assertEqual(remote_run.call_args.args[1], "set -euo pipefail\nprintf 'safe'\n")


if __name__ == "__main__":
    unittest.main()
