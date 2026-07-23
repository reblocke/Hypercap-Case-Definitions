from __future__ import annotations

import json
import stat
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import stata_run  # noqa: E402


def write_fake_stata(
    path: Path,
    *,
    write_status: bool = True,
    write_guarded_controls: bool = True,
) -> None:
    xlsx = repr(list(stata_run.EXPECTED_XLSX))
    png = repr(list(stata_run.EXPECTED_PNG))
    gph = repr(list(stata_run.EXPECTED_GPH))
    source = f"""#!/usr/bin/env python3
import pathlib
import sys

def clean(value):
    return value.strip('"')

output_root = pathlib.Path(clean(sys.argv[7]))
run_id = clean(sys.argv[8])
status_path = pathlib.Path(clean(sys.argv[9]))
dependency_path = pathlib.Path(clean(sys.argv[10]))
mode = clean(sys.argv[11])
artifact_root = output_root / run_id
artifact_root.mkdir(parents=True, exist_ok=True)
(artifact_root / "graph-temp").mkdir(exist_ok=True)
(artifact_root / "Logs").mkdir(exist_ok=True)
for name in {xlsx}:
    (artifact_root / name).write_bytes(b"x")
for name in {png}:
    (artifact_root / name).write_bytes(b"x")
for name in {gph}:
    (artifact_root / "graph-temp" / name).write_bytes(b"x")
(artifact_root / "Logs" / "analysis.log").write_text("log\\n")
(artifact_root / "Logs" / "analysis.do").write_bytes(b"x")
(artifact_root / "ANALYSIS_COMPLETE").write_text("analysis_complete=true\\n")
dependency_path.write_text(
    "name\\tkind\\trequired\\tfound\\tresolved_path\\n"
    "test\\tcommand\\ttrue\\ttrue\\t/opt/stata/test.ado\\n"
)
"""
    if write_guarded_controls:
        source += """(artifact_root / "input_validation.tsv").write_text(
    "check\\tstatus\\n"
)
(artifact_root / "run_metrics.tsv").write_text("metric\\tvalue\\n")
"""
    if write_status:
        source += """status_path.write_text(
    "status\\tsuccess\\n"
    "dependency_rc\\t0\\n"
    "analysis_rc\\t0\\n"
    "completion_marker_found\\t1\\n"
    "stata_version\\t18\\n"
    "stata_flavor\\tIC\\n"
    "operating_system\\tMacOSX\\n"
    "machine_type\\tTest\\n"
)
"""
    path.write_text(source, encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IXUSR)


class StataRunUnitTests(unittest.TestCase):
    def test_run_id_rejects_path_characters(self) -> None:
        self.assertIsNone(stata_run.RUN_ID_PATTERN.fullmatch("../unsafe"))

    def test_inventory_reports_missing_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            inventory = stata_run.validate_artifact_inventory(Path(tmp))
        self.assertEqual(40, inventory["expected_legacy_count"])
        self.assertTrue(inventory["missing"])

    def test_macos_app_uses_exit_batch_mode(self) -> None:
        path = Path("/Applications/Stata/StataBE.app/Contents/MacOS/StataBE")
        self.assertEqual(("-e", "macos-e"), stata_run.stata_mode(path, "auto"))

    def test_zero_process_exit_without_status_fails(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            input_root = base / "input"
            output_root = base / "output"
            input_root.mkdir()
            (input_root / "full_db.dta").write_bytes(b"not-patient-data")
            fake = base / "fake-stata"
            write_fake_stata(fake, write_status=False)
            args = stata_run.parser().parse_args(
                [
                    "--input-root",
                    str(input_root),
                    "--output-root",
                    str(output_root),
                    "--stata-bin",
                    str(fake),
                    "--run-id",
                    "missing-status",
                ]
            )
            code, run_dir = stata_run.run(args)
            manifest = json.loads((run_dir / "run_manifest.json").read_text())
            success_exists = (run_dir / "SUCCESS").exists()
        self.assertEqual(1, code)
        self.assertEqual("missing_status", manifest["failure_category"])
        self.assertFalse(success_exists)

    def test_success_manifest_contains_no_absolute_paths(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            input_root = base / "input"
            output_root = base / "output"
            input_root.mkdir()
            (input_root / "full_db.dta").write_bytes(b"not-patient-data")
            fake = base / "fake-stata"
            write_fake_stata(fake)
            args = stata_run.parser().parse_args(
                [
                    "--input-root",
                    str(input_root),
                    "--output-root",
                    str(output_root),
                    "--stata-bin",
                    str(fake),
                    "--run-id",
                    "successful-run",
                ]
            )
            code, run_dir = stata_run.run(args)
            manifest_text = (run_dir / "run_manifest.json").read_text()
            manifest = json.loads(manifest_text)
            success_exists = (run_dir / "SUCCESS").is_file()
        self.assertEqual(0, code)
        self.assertEqual("success", manifest["status"])
        self.assertNotIn(str(base), manifest_text)
        self.assertNotIn("/opt/stata", manifest_text)
        self.assertTrue(success_exists)

    def test_success_is_removed_when_control_artifacts_are_incomplete(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            input_root = base / "input"
            output_root = base / "output"
            input_root.mkdir()
            (input_root / "full_db.dta").write_bytes(b"not-patient-data")
            fake = base / "fake-stata"
            write_fake_stata(fake, write_guarded_controls=False)
            args = stata_run.parser().parse_args(
                [
                    "--input-root",
                    str(input_root),
                    "--output-root",
                    str(output_root),
                    "--stata-bin",
                    str(fake),
                    "--run-id",
                    "incomplete-controls",
                ]
            )
            code, run_dir = stata_run.run(args)
            manifest = json.loads((run_dir / "run_manifest.json").read_text())
            success_exists = (run_dir / "SUCCESS").exists()
        self.assertEqual(1, code)
        self.assertEqual("incomplete_artifacts", manifest["failure_category"])
        self.assertFalse(success_exists)

    def test_existing_run_directory_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            input_root = base / "input"
            output_root = base / "output"
            input_root.mkdir()
            output_root.mkdir()
            (input_root / "full_db.dta").write_bytes(b"not-patient-data")
            (output_root / "collision").mkdir()
            fake = base / "fake-stata"
            write_fake_stata(fake)
            args = stata_run.parser().parse_args(
                [
                    "--input-root",
                    str(input_root),
                    "--output-root",
                    str(output_root),
                    "--stata-bin",
                    str(fake),
                    "--run-id",
                    "collision",
                ]
            )
            with self.assertRaises(stata_run.RunFailure) as context:
                stata_run.run(args)
        self.assertEqual("run_collision", context.exception.category)


class AnalysisBodyLockTests(unittest.TestCase):
    def test_scientific_body_matches_hcd_000a_except_row_listing(self) -> None:
        import hashlib

        text = (ROOT / stata_run.MAIN_DO).read_text(encoding="utf-8")
        start_marker = "/* ------------------\n   Pre-processing"
        end_marker = (
            "graph export \"`outdir'/Location - Figure S3 Prob Hypercap ICD.png\", "
            'name("Graph") width(3200) replace'
        )
        start = text.index(start_marker)
        end = text.index(end_marker, start) + len(end_marker)
        body = text[start:end] + "\n"
        observed = hashlib.sha256(body.encode("utf-8")).hexdigest()
        expected = "9df22eae3edbb3bb369218550452d3023b86bec96143da5260906031eaa493b6"
        self.assertEqual(expected, observed)
        self.assertNotIn("in 1/200", body)


if __name__ == "__main__":
    unittest.main()
