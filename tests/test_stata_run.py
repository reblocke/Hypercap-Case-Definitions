from __future__ import annotations

import json
import stat
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import stata_run  # noqa: E402
import input_manifest  # noqa: E402


def write_approved_input(input_root: Path) -> dict[str, object]:
    (input_root / "full_db.dta").write_bytes(b"not-patient-data")
    manifest_path = input_manifest.approve_input(
        input_root,
        ROOT,
        approval="YES",
    )
    return json.loads(manifest_path.read_text(encoding="utf-8"))


def write_fake_stata(
    path: Path,
    *,
    write_status: bool = True,
    write_guarded_controls: bool = True,
    mutate_input: bool = False,
    exit_code: int = 0,
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
input_root = pathlib.Path(clean(sys.argv[6]))
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
    if mutate_input:
        source += """(input_root / "full_db.dta").write_bytes(
    b"changed-during-stata-run"
)
"""
    if exit_code:
        source += f"raise SystemExit({exit_code})\n"
    path.write_text(source, encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IXUSR)


class StataRunUnitTests(unittest.TestCase):
    def test_run_id_rejects_path_characters(self) -> None:
        self.assertIsNone(stata_run.RUN_ID_PATTERN.fullmatch("../unsafe"))

    def test_analysis_root_must_match_runner_checkout_before_run_creation(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            analysis_root = base / "other-checkout"
            output_root = base / "output"
            analysis_root.mkdir()
            args = stata_run.parser().parse_args(
                [
                    "--analysis-root",
                    str(analysis_root),
                    "--output-root",
                    str(output_root),
                ]
            )
            with self.assertRaises(stata_run.RunFailure) as raised:
                stata_run.run(args)
        self.assertEqual("analysis_root_mismatch", raised.exception.category)
        self.assertFalse(output_root.exists())

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
            write_approved_input(input_root)
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

    def test_nonzero_process_exit_rejects_successful_status_and_artifacts(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            input_root = base / "input"
            output_root = base / "output"
            input_root.mkdir()
            write_approved_input(input_root)
            fake = base / "fake-stata"
            write_fake_stata(fake, exit_code=7)
            args = stata_run.parser().parse_args(
                [
                    "--input-root",
                    str(input_root),
                    "--output-root",
                    str(output_root),
                    "--stata-bin",
                    str(fake),
                    "--run-id",
                    "nonzero-exit",
                ]
            )
            code, run_dir = stata_run.run(args)
            manifest = json.loads((run_dir / "run_manifest.json").read_text())
            success_exists = (run_dir / "SUCCESS").exists()
        self.assertEqual(1, code)
        self.assertEqual("stata_process_exit", manifest["failure_category"])
        self.assertEqual(7, manifest["stata"]["process_return_code"])
        self.assertFalse(success_exists)

    def test_negative_process_return_rejects_successful_status_and_artifacts(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            input_root = base / "input"
            output_root = base / "output"
            input_root.mkdir()
            write_approved_input(input_root)
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
                    "negative-exit",
                ]
            )
            real_subprocess_run = stata_run.subprocess.run

            def force_negative_return(command, **kwargs):
                result = real_subprocess_run(command, **kwargs)
                if len(command) > 2 and command[2] == "do":
                    result.returncode = -9
                return result

            with mock.patch.object(
                stata_run.subprocess,
                "run",
                side_effect=force_negative_return,
            ):
                code, run_dir = stata_run.run(args)
            manifest = json.loads((run_dir / "run_manifest.json").read_text())
            success_exists = (run_dir / "SUCCESS").exists()
        self.assertEqual(1, code)
        self.assertEqual("stata_process_exit", manifest["failure_category"])
        self.assertEqual(-9, manifest["stata"]["process_return_code"])
        self.assertFalse(success_exists)

    def test_success_manifest_contains_no_absolute_paths(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            input_root = base / "input"
            output_root = base / "output"
            input_root.mkdir()
            approval = write_approved_input(input_root)
            approval_manifest = input_root / "full_db.manifest.json"
            approval_reference = {
                "manifest_sha256": stata_run.sha256_file(approval_manifest),
                "upstream_repository": approval["upstream_repository"],
                "producer_commit": approval["producer_commit"],
                "input_schema_version": approval["input_schema_version"],
                "data_dictionary_sha256": approval["data_dictionary_sha256"],
            }
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
        self.assertEqual(2, manifest["schema_version"])
        self.assertEqual("success", manifest["status"])
        self.assertEqual(approval_reference, manifest["input"]["approval"])
        self.assertNotIn(str(base), manifest_text)
        self.assertNotIn("/opt/stata", manifest_text)
        self.assertTrue(success_exists)

    def test_input_mutation_during_stata_prevents_success(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            input_root = base / "input"
            output_root = base / "output"
            input_root.mkdir()
            write_approved_input(input_root)
            fake = base / "fake-stata"
            write_fake_stata(fake, mutate_input=True)
            args = stata_run.parser().parse_args(
                [
                    "--input-root",
                    str(input_root),
                    "--output-root",
                    str(output_root),
                    "--stata-bin",
                    str(fake),
                    "--run-id",
                    "mutated-input",
                ]
            )
            code, run_dir = stata_run.run(args)
            manifest = json.loads((run_dir / "run_manifest.json").read_text())
            success_exists = (run_dir / "SUCCESS").exists()
        self.assertEqual(1, code)
        self.assertEqual(
            "input_manifest_input_mismatch",
            manifest["failure_category"],
        )
        self.assertFalse(success_exists)

    def test_reapproved_input_during_stata_is_detected_as_run_drift(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            input_root = base / "input"
            output_root = base / "output"
            input_root.mkdir()
            write_approved_input(input_root)
            initial_approval = input_manifest.validate_approved_input(
                input_root,
                ROOT,
            )
            replacement_approval = {
                "input_file": initial_approval["input_file"],
                "input": {
                    **initial_approval["input"],
                    "sha256": "f" * 64,
                },
                "approval": {
                    **initial_approval["approval"],
                    "manifest_sha256": "e" * 64,
                },
            }
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
                    "reapproved-input",
                ]
            )
            with mock.patch.object(
                stata_run.input_manifest,
                "validate_approved_input",
                side_effect=[initial_approval, replacement_approval],
            ):
                code, run_dir = stata_run.run(args)
            manifest = json.loads((run_dir / "run_manifest.json").read_text())
            success_exists = (run_dir / "SUCCESS").exists()
        self.assertEqual(1, code)
        self.assertEqual("input_manifest_run_drift", manifest["failure_category"])
        self.assertFalse(success_exists)

    def test_success_is_removed_when_control_artifacts_are_incomplete(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            input_root = base / "input"
            output_root = base / "output"
            input_root.mkdir()
            write_approved_input(input_root)
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

    def test_missing_approval_fails_before_stata_resolution_in_all_modes(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            input_root = base / "input"
            output_root = base / "output"
            input_root.mkdir()
            (input_root / "full_db.dta").write_bytes(b"not-patient-data")
            for run_id, extra_args in (
                ("missing-approval", []),
                ("missing-approval-legacy", ["--legacy-two-arg"]),
            ):
                with self.subTest(run_id=run_id):
                    args = stata_run.parser().parse_args(
                        [
                            "--input-root",
                            str(input_root),
                            "--output-root",
                            str(output_root),
                            "--run-id",
                            run_id,
                            *extra_args,
                        ]
                    )
                    with mock.patch.object(
                        stata_run,
                        "resolve_stata",
                    ) as resolve_stata:
                        code, run_dir = stata_run.run(args)
                    manifest_text = (run_dir / "run_manifest.json").read_text()
                    manifest = json.loads(manifest_text)
                    self.assertEqual(1, code)
                    self.assertEqual(
                        "input_manifest_missing",
                        manifest["failure_category"],
                    )
                    self.assertNotIn(str(base), manifest_text)
                    resolve_stata.assert_not_called()

    def test_malformed_approval_fails_before_stata_resolution(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            input_root = base / "input"
            output_root = base / "output"
            input_root.mkdir()
            (input_root / "full_db.dta").write_bytes(b"not-patient-data")
            (input_root / "full_db.manifest.json").write_text("{not-json\n")
            args = stata_run.parser().parse_args(
                [
                    "--input-root",
                    str(input_root),
                    "--output-root",
                    str(output_root),
                    "--run-id",
                    "malformed-approval",
                ]
            )
            with mock.patch.object(stata_run, "resolve_stata") as resolve_stata:
                code, run_dir = stata_run.run(args)
            manifest_text = (run_dir / "run_manifest.json").read_text()
            manifest = json.loads(manifest_text)
        self.assertEqual(1, code)
        self.assertEqual("input_manifest_malformed", manifest["failure_category"])
        self.assertNotIn(str(base), manifest_text)
        resolve_stata.assert_not_called()

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


class InputValidationSourceTests(unittest.TestCase):
    def test_input_contract_import_error_is_saved_before_cleanup(self) -> None:
        text = (ROOT / "stata" / "validate_input.do").read_text(encoding="utf-8")
        import_marker = "capture frame hcd_contract: import delimited"
        save_marker = "local import_rc = _rc"
        failure_marker = "if `import_rc' {"
        cleanup_marker = "capture frame drop hcd_contract"
        exit_marker = "exit `import_rc'"

        import_at = text.index(import_marker)
        save_at = text.index(save_marker, import_at)
        failure_at = text.index(failure_marker, save_at)
        cleanup_at = text.index(cleanup_marker, failure_at)
        exit_at = text.index(exit_marker, cleanup_at)

        self.assertLess(import_at, save_at)
        self.assertIn(
            "varnames(1) stringcols(_all) clear\nlocal import_rc = _rc",
            text[import_at : save_at + len(save_marker)],
        )
        self.assertLess(save_at, failure_at)
        self.assertLess(failure_at, cleanup_at)
        self.assertLess(cleanup_at, exit_at)
        self.assertNotIn("exit _rc", text[failure_at:exit_at])

    def test_verified_upstream_aggregates_are_contract_checked(self) -> None:
        text = (ROOT / "stata" / "validate_input.do").read_text(encoding="utf-8")
        self.assertIn(
            "cond(missing(hypercap_on_abg), 0, hypercap_on_abg) != ///",
            text,
        )
        self.assertIn("(paco2 >= 45 & !missing(paco2))", text)
        self.assertIn(
            "cond(missing(hypercap_resp_failure), 0, ///",
            text,
        )
        self.assertIn(
            "(ohs_code == 1 | has_j9602 == 1 | has_j9612 == 1 | ///",
            text,
        )
        aggregate_block = text[
            text.index("quietly count if cond(missing(hypercap_on_abg)") :
            text.index("quietly count if has_abg == 0")
        ]
        self.assertEqual(2, aggregate_block.count('"aggregate_definition"'))
        self.assertEqual(
            2,
            aggregate_block.count("if `violations' > 0 local ++hard_failures"),
        )


class AnalysisBodyLockTests(unittest.TestCase):
    def test_diagnostic_program_is_locked(self) -> None:
        import hashlib

        text = (ROOT / stata_run.MAIN_DO).read_text(encoding="utf-8")
        start_marker = "capture program drop test_char_from_icd_list"
        definition_marker = "program define test_char_from_icd_list"
        start = text.index(start_marker)
        definition = text.index(definition_marker, start)
        end = text.index("\nend", definition) + len("\nend")
        program = text[start:end] + "\n"
        observed = hashlib.sha256(program.encode("utf-8")).hexdigest()
        expected = "9635569dae8be5253639a7abeb69e41de62b486a0fe3b7bfc067c58b92f19225"

        self.assertEqual(expected, observed)
        self.assertIn("quietly diagt `ref_std' `icdcode'", program)
        self.assertIn("quietly roctab `ref_std' `icdcode', nograph", program)
        mutated = program.replace("local auc = string(r(area)", "local auc = string(r(lb)", 1)
        self.assertNotEqual(
            expected,
            hashlib.sha256(mutated.encode("utf-8")).hexdigest(),
        )

    def test_scientific_body_matches_approved_hcd_001_correction(self) -> None:
        import hashlib

        text = (ROOT / stata_run.MAIN_DO).read_text(encoding="utf-8")
        start_marker = "/* ------------------\n   Pre-processing"
        end_marker = (
            "graph export \"`outdir'/Location - e-Figure 5 Prob Hypercap ICD.png\", "
            'name("Graph") width(3200) replace'
        )
        start = text.index(start_marker)
        end = text.index(end_marker, start) + len(end_marker)
        body = text[start:end] + "\n"
        observed = hashlib.sha256(body.encode("utf-8")).hexdigest()
        expected = "7eb9c366e27c1e146ca0fa09807fbf1fe9bf95ee8eed6b598b68c0e4a17484b2"
        self.assertEqual(expected, observed)
        self.assertNotIn("in 1/200", body)
        self.assertIn(
            "gen def4 = (paco2 > 45 & abg_ph < 7.35 & niv_proc == 1)",
            body,
        )
        self.assertIn("gen def5 = (paco2 >= 45) if !missing(paco2)", body)
        self.assertIn("gen def8 = (paco2 >= 50", body)
        self.assertIn("vbg_ph > 7.35", body)
        self.assertIn("abg_ph <= 7.45", body)
        self.assertIn("keep if abg_vbg_confusion_matrix == 1", body)
        self.assertIn("keep if abg_vbg_confusion_matrix == 2", body)
        self.assertIn("keep if abg_vbg_confusion_matrix == 3", body)
        self.assertEqual(
            3,
            body.count(
                "Kappa is undefined when either definition is constant "
                "in the subgroup."
            ),
        )
        self.assertEqual(3, body.count("local def_i_varies = r(N) > 0"))
        self.assertEqual(3, body.count("local def_j_varies = r(N) > 0"))
        self.assertEqual(2, body.count("assert r(N) == 4"))
        self.assertIn("Figure 3 Prob Hypercap ICD.png", body)
        self.assertIn("Location - e-Figure 5 Prob Hypercap ICD.png", body)
        self.assertIn("version 17.0", text)


if __name__ == "__main__":
    unittest.main()
