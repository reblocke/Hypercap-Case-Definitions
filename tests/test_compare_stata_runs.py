from __future__ import annotations

import contextlib
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from openpyxl import Workbook
from PIL import Image
from PIL.PngImagePlugin import PngInfo

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import compare_stata_runs as compare  # noqa: E402

APPROVAL_REFERENCE = {
    "manifest_sha256": "1" * 64,
    "upstream_repository": "reblocke/trinetx-hypercapnia-code",
    "producer_commit": "2" * 40,
    "input_schema_version": "hypercapnia-full-db-v1",
    "data_dictionary_sha256": "3" * 64,
}


def save_workbook(path: Path, value: object) -> None:
    workbook = Workbook()
    sheet = workbook.active
    sheet.title = "Results"
    sheet["A1"] = "value"
    sheet["A2"] = value
    workbook.save(path)
    workbook.close()


def write_run_manifest(
    path: Path,
    *,
    legacy: bool,
    commit: str,
    input_sha256: str,
    run_id: str,
    approval: dict[str, str] | None = None,
    include_approval: bool = True,
    artifact_relative_path: str = ".",
) -> None:
    path.mkdir()
    input_record: dict[str, object] = {"sha256": input_sha256}
    if include_approval:
        input_record["approval"] = (
            approval if approval is not None else APPROVAL_REFERENCE
        )
    (path / "run_manifest.json").write_text(
        json.dumps(
            {
                "run_id": run_id,
                "status": "success",
                "legacy_two_argument_mode": legacy,
                "analysis_commit": commit,
                "analysis_worktree_dirty": False,
                "input": input_record,
                "artifact_relative_path": artifact_relative_path,
                "artifacts": {
                    "missing": [],
                    "control_missing": [],
                    "analysis_log": "Logs/analysis.log",
                    "copied_do": "Logs/analysis.do",
                },
            }
        ),
        encoding="utf-8",
    )


def write_complete_run_evidence(path: Path, *, legacy: bool) -> None:
    (path / "Logs").mkdir()
    (path / "graph-temp").mkdir()
    for name in compare.stata_run.EXPECTED_XLSX:
        (path / name).write_bytes(b"x")
    for name in compare.stata_run.EXPECTED_PNG:
        observed_name = (
            compare.stata_run.LEGACY_PNG_ALIASES.get(name, name)
            if legacy
            else name
        )
        (path / observed_name).write_bytes(b"x")
    for name in compare.stata_run.EXPECTED_GPH:
        (path / "graph-temp" / name).write_bytes(b"x")
    (path / "Logs" / "analysis.log").write_text("log\n", encoding="utf-8")
    (path / "Logs" / "analysis.do").write_text("do\n", encoding="utf-8")
    (path / "dependency_report.tsv").write_text("dependency\n", encoding="utf-8")
    (path / "run_status.tsv").write_text("status\tsuccess\n", encoding="utf-8")
    (path / "SUCCESS").write_text("status=success\n", encoding="utf-8")
    if not legacy:
        (path / "input_validation.tsv").write_text("check\tstatus\n", encoding="utf-8")
        (path / "run_metrics.tsv").write_text("metric\tvalue\n", encoding="utf-8")
        (path / "ANALYSIS_COMPLETE").write_text(
            "analysis_complete=true\n",
            encoding="utf-8",
        )


def passing_comparison(
    _left: Path,
    _right: Path,
    label: str,
) -> dict[str, object]:
    return {"comparison": label, "status": "pass", "failures": []}


def write_approved_input_fixture(root: Path) -> tuple[Path, Path, dict[str, object]]:
    analysis_root = root / "analysis"
    metadata_root = analysis_root / "metadata"
    input_root = root / "input"
    metadata_root.mkdir(parents=True)
    input_root.mkdir()
    (metadata_root / "upstream_dependency.yml").write_text(
        (ROOT / "metadata/upstream_dependency.yml").read_text(encoding="utf-8"),
        encoding="utf-8",
    )
    (analysis_root / "data_dictionary.csv").write_text(
        (ROOT / "data_dictionary.csv").read_text(encoding="utf-8"),
        encoding="utf-8",
    )
    input_file = input_root / "full_db.dta"
    input_file.write_bytes(b"current input")
    compare.input_manifest.approve_input(
        input_root,
        analysis_root,
        approval="YES",
    )
    approved = compare.input_manifest.validate_approved_input(
        input_root,
        analysis_root,
    )
    return input_file, analysis_root, approved


def run_stubbed_comparator(
    root: Path,
    *,
    legacy_overrides: dict[str, bool] | None = None,
    commit_overrides: dict[str, str] | None = None,
    input_sha256: str | None = None,
    input_hash_overrides: dict[str, str] | None = None,
    input_file: Path | None = None,
    current_approval: dict[str, str] | None = None,
    approval_overrides: dict[str, dict[str, str]] | None = None,
    omit_approval_roles: set[str] | None = None,
    run_id_overrides: dict[str, str] | None = None,
    artifact_relative_overrides: dict[str, str] | None = None,
    missing_control: tuple[str, str] | None = None,
    comparison_mode: str = "equivalence",
    comparison_side_effect: object = passing_comparison,
) -> tuple[int, dict[str, object], int]:
    roles = ("baseline", "candidate_1", "candidate_2")
    if input_file is None:
        input_file = root / "full_db.dta"
        input_file.write_bytes(b"current input")
    if input_sha256 is None:
        input_sha256 = compare.stata_run.sha256_file(input_file)

    legacy_modes = {
        "baseline": True,
        "candidate_1": False,
        "candidate_2": False,
    }
    legacy_modes.update(legacy_overrides or {})
    commits = {
        "baseline": "baseline-commit",
        "candidate_1": "candidate-commit",
        "candidate_2": "candidate-commit",
    }
    commits.update(commit_overrides or {})
    input_hashes = {role: input_sha256 for role in roles}
    input_hashes.update(input_hash_overrides or {})
    current_approval = current_approval or APPROVAL_REFERENCE
    approvals = {role: current_approval for role in roles}
    approvals.update(approval_overrides or {})
    omit_approval_roles = omit_approval_roles or set()
    run_ids = {role: f"{role}-run" for role in roles}
    run_ids.update(run_id_overrides or {})
    artifact_relatives = {role: "." for role in roles}
    artifact_relatives.update(artifact_relative_overrides or {})

    runs = {role: root / role for role in roles}
    for role, path in runs.items():
        write_run_manifest(
            path,
            legacy=legacy_modes[role],
            commit=commits[role],
            input_sha256=input_hashes[role],
            run_id=run_ids[role],
            approval=approvals[role],
            include_approval=role not in omit_approval_roles,
            artifact_relative_path=artifact_relatives[role],
        )
        write_complete_run_evidence(path, legacy=legacy_modes[role])
    if missing_control is not None:
        role, relative_path = missing_control
        (runs[role] / relative_path).unlink()

    report_path = root / "comparison.json"
    arguments = [
        "--baseline-run",
        str(runs["baseline"]),
        "--candidate-run-1",
        str(runs["candidate_1"]),
        "--candidate-run-2",
        str(runs["candidate_2"]),
        "--report",
        str(report_path),
        "--input-file",
        str(input_file),
        "--comparison-mode",
        comparison_mode,
    ]
    with (
        patch.object(
            compare.input_manifest,
            "validate_approved_input",
            return_value={
                "input_file": input_file,
                "input": {
                    "logical_name": "full_db.dta",
                    "size_bytes": input_file.stat().st_size,
                    "sha256": input_sha256,
                },
                "approval": current_approval,
            },
        ),
        patch.object(
            compare,
            "compare_pair",
            side_effect=comparison_side_effect,
        ) as compare_pair,
    ):
        status = compare.main(arguments)
    report = json.loads(report_path.read_text(encoding="utf-8"))
    return status, report, compare_pair.call_count


class ComparatorUnitTests(unittest.TestCase):
    def test_identical_workbooks_pass(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            left = root / "left.xlsx"
            right = root / "right.xlsx"
            save_workbook(left, 1.25)
            save_workbook(right, 1.25)
            failures: list[dict[str, object]] = []
            compare.compare_workbook(left, right, failures)
        self.assertEqual([], failures)

    def test_numeric_workbook_difference_is_detected_without_values(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            left = root / "left.xlsx"
            right = root / "right.xlsx"
            save_workbook(left, 1.25)
            save_workbook(right, 1.5)
            failures: list[dict[str, object]] = []
            compare.compare_workbook(left, right, failures)
        self.assertTrue(failures)
        self.assertNotIn("1.25", repr(failures))
        self.assertNotIn("1.5", repr(failures))

    def test_png_metadata_is_ignored_but_pixels_are_compared(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            left = root / "left.png"
            right = root / "right.png"
            changed = root / "changed.png"
            image = Image.new("RGB", (3, 3), "white")
            image.save(left)
            metadata = PngInfo()
            metadata.add_text("created", "different")
            image.save(right, pnginfo=metadata)
            Image.new("RGB", (3, 3), "black").save(changed)
            equal_failures: list[dict[str, object]] = []
            changed_failures: list[dict[str, object]] = []
            compare.compare_png(left, right, equal_failures)
            compare.compare_png(left, changed, changed_failures)
        self.assertEqual([], equal_failures)
        self.assertTrue(changed_failures)

    def test_legacy_inventory_uses_explicit_png_aliases(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            write_complete_run_evidence(root, legacy=True)
            legacy_inventory = compare.stata_run.validate_artifact_inventory(
                root,
                legacy=True,
            )
            current_inventory = compare.stata_run.validate_artifact_inventory(
                root,
                legacy=False,
            )
        self.assertEqual([], legacy_inventory["missing"])
        for current_name, legacy_name in compare.stata_run.LEGACY_PNG_ALIASES.items():
            self.assertIn(current_name, current_inventory["missing"])
            self.assertNotEqual(current_name, legacy_name)

    def test_log_normalization_removes_restricted_row_listing(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            baseline = root / "baseline.log"
            candidate = root / "candidate.log"
            prefix = "/* ------------------\n   Pre-processing\n--------------------*/\n"
            suffix = (
                "/* Ones I can't do */\n"
                'graph export "`outdir\'/Location - Figure S3 Prob Hypercap ICD.png", '
                'name("Graph") width(3200) replace\n'
            )
            baseline.write_text(
                prefix
                + "list hypercap_resp_failure paco2 vbg_co2 in 1/200\n"
                + "restricted row output\n"
                + suffix,
                encoding="utf-8",
            )
            candidate.write_text(prefix + suffix, encoding="utf-8")
            left = compare.normalize_log(baseline, root, root)
            right = compare.normalize_log(candidate, root, root)
        self.assertEqual(left, right)
        self.assertNotIn("restricted row output", left)

    def test_log_normalization_excludes_driver_tail(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            baseline = root / "baseline.log"
            candidate = root / "candidate.log"
            analysis = (
                "/* ------------------\n"
                "   Pre-processing\n"
                "--------------------*/\n"
                ". display \"analysis\"\n"
            )
            baseline.write_text(
                analysis + "\nend of do-file\n. display \"driver\"\n",
                encoding="utf-8",
            )
            candidate.write_text(analysis, encoding="utf-8")
            left = compare.normalize_log(baseline, root, root)
            right = compare.normalize_log(candidate, root, root)
        self.assertEqual(left, right)
        self.assertNotIn("driver", left)

    def test_log_normalization_is_safe_after_run_relocation(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            relocated_left = root / "archive" / "candidate-one"
            relocated_right = root / "archive" / "candidate-two"
            relocated_left.mkdir(parents=True)
            relocated_right.mkdir(parents=True)
            left_log = relocated_left / "analysis.log"
            right_log = relocated_right / "analysis.log"
            prefix = "/* ------------------\n   Pre-processing\n--------------------*/\n"
            left_log.write_text(
                prefix
                + "file /original/run-one/Overall Cohort chars.xlsx saved\n"
                + "(file /original/run-one/graph-temp/"
                + "All_Encounters_Prob_Dx_spline.gph not found)\n",
                encoding="utf-8",
            )
            right_log.write_text(
                prefix
                + "file /original/run-two/Overall Cohort chars.xlsx saved\n"
                + "(file /original/run-two/graph-temp/"
                + "All_Encounters_Prob_Dx_spline.gph not found)\n",
                encoding="utf-8",
            )
            left = compare.normalize_log(
                left_log,
                relocated_left,
                relocated_left,
            )
            right = compare.normalize_log(
                right_log,
                relocated_right,
                relocated_right,
            )
        self.assertEqual(left, right)
        self.assertIn("file <OUTPUT>/Overall Cohort chars.xlsx saved", left)
        self.assertIn(
            "(file <OUTPUT>/graph-temp/All_Encounters_Prob_Dx_spline.gph "
            "not found)",
            left,
        )
        self.assertNotIn("/original/run-one", left)
        self.assertNotIn("/original/run-two", right)

    def test_log_normalization_does_not_hide_unknown_path_differences(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            left_log = root / "left.log"
            right_log = root / "right.log"
            prefix = "/* ------------------\n   Pre-processing\n--------------------*/\n"
            left_log.write_text(
                prefix + "file /original/run-one/unexpected.txt saved\n",
                encoding="utf-8",
            )
            right_log.write_text(
                prefix + "file /original/run-two/unexpected.txt saved\n",
                encoding="utf-8",
            )
            left = compare.normalize_log(left_log, root, root)
            right = compare.normalize_log(right_log, root, root)
        self.assertNotEqual(left, right)

    def test_log_normalization_rejects_mixed_expected_artifact_roots(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            left_log = root / "left.log"
            right_log = root / "right.log"
            prefix = "/* ------------------\n   Pre-processing\n--------------------*/\n"
            left_log.write_text(
                prefix
                + "file /original/run-one/Overall Cohort chars.xlsx saved\n"
                + "file /original/run-one/Definition Overlap HeatPlot.png "
                + "saved as PNG format\n",
                encoding="utf-8",
            )
            right_log.write_text(
                prefix
                + "file /original/run-two/Overall Cohort chars.xlsx saved\n"
                + "file /unexpected/run-three/Definition Overlap HeatPlot.png "
                + "saved as PNG format\n",
                encoding="utf-8",
            )
            left = compare.normalize_log(left_log, root, root)
            with self.assertRaises(compare.LogOutputRootFailure):
                compare.normalize_log(right_log, root, root)
        self.assertIn("<OUTPUT>/Overall Cohort chars.xlsx", left)

    def test_log_comparison_fails_for_identical_mixed_output_roots(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            left = root / "left"
            right = root / "right"
            for run in (left, right):
                (run / "Logs").mkdir(parents=True)
                (run / "Logs" / "analysis.log").write_text(
                    "/* ------------------\n"
                    "   Pre-processing\n"
                    "--------------------*/\n"
                    "file /shared/run-one/Overall Cohort chars.xlsx saved\n"
                    "file /shared/run-two/Definition Overlap HeatPlot.png "
                    "saved as PNG format\n",
                    encoding="utf-8",
                )
            manifest = {"artifacts": {"analysis_log": "Logs/analysis.log"}}
            failures: list[dict[str, object]] = []
            compare.compare_logs(
                manifest,
                left,
                left,
                manifest,
                right,
                right,
                failures,
            )
        self.assertEqual(
            [
                {
                    "category": "analysis_log_output_roots",
                    "artifact": "analysis_log",
                    "role": "left",
                },
                {
                    "category": "analysis_log_output_roots",
                    "artifact": "analysis_log",
                    "role": "right",
                },
            ],
            failures,
        )

    def test_known_artifact_paths_remain_sensitive_outside_file_messages(self) -> None:
        left = (
            'display "/original/run-one/Overall Cohort chars.xlsx"\n'
            "file /original/run-one/Overall Cohort chars.xlsx saved"
        )
        right = (
            'display "/original/run-two/Overall Cohort chars.xlsx"\n'
            "file /original/run-two/Overall Cohort chars.xlsx saved"
        )
        self.assertNotEqual(
            compare.redact_known_output_paths(left),
            compare.redact_known_output_paths(right),
        )
        self.assertIn(
            'display "/original/run-one/Overall Cohort chars.xlsx"',
            compare.redact_known_output_paths(left),
        )
        self.assertIn(
            "file <OUTPUT>/Overall Cohort chars.xlsx saved",
            compare.redact_known_output_paths(left),
        )

    def test_known_output_path_redaction_supports_windows_paths(self) -> None:
        text = (
            r"file C:\original\run-one\graph-temp"
            r"\All_Encounters_Prob_Dx_spline.gph saved"
        )
        self.assertEqual(
            "file <OUTPUT>/graph-temp/All_Encounters_Prob_Dx_spline.gph saved",
            compare.redact_known_output_paths(text),
        )

    def test_stata_wrapped_path_is_reassembled_before_redaction(self) -> None:
        lines = [
            "file /approved/output/hcd0",
            "> 00b/Overall Cohort chars.xlsx saved",
        ]
        self.assertEqual(
            ["file /approved/output/hcd000b/Overall Cohort chars.xlsx saved"],
            compare.unwrap_stata_continuations(lines),
        )

    def test_indented_stata_continuation_is_reassembled(self) -> None:
        lines = ["(file /approved/All_Prob", "    > _Dx.gph not found)"]
        self.assertEqual(
            ["(file /approved/All_Prob_Dx.gph not found)"],
            compare.unwrap_stata_continuations(lines),
        )

    def test_parenthesized_file_message_is_reassembled(self) -> None:
        lines = [
            "(file /approved/Research",
            "    Projects/output.gph not",
            "    found)",
            "",
        ]
        self.assertEqual(
            ["(file /approved/Research Projects/output.gph not found)", ""],
            compare.unwrap_stata_parenthesized_file_messages(lines),
        )

    def test_path_redaction_handles_unmarked_stata_line_wrap(self) -> None:
        path = Path("/approved/Research Projects/output/run-one")
        text = "file /approved/Research\nProjects/output/run-one/figure.png saved\n"
        self.assertEqual(
            "file <OUTPUT>/figure.png saved\n",
            compare.redact_wrapped_path(text, path, "<OUTPUT>"),
        )

    def test_wrapped_file_saved_message_is_reassembled(self) -> None:
        lines = [
            "file <OUTPUT>/Definition Overlap",
            "HeatPlot.png saved as PNG format",
            "",
        ]
        self.assertEqual(
            [
                "file <OUTPUT>/Definition Overlap HeatPlot.png saved as PNG format",
                "",
            ],
            compare.unwrap_stata_file_messages(lines),
        )

    def test_wrapped_png_format_suffix_is_reassembled(self) -> None:
        lines = [
            "file <OUTPUT>/figure.png saved as PNG",
            "format",
            "",
        ]
        self.assertEqual(
            ["file <OUTPUT>/figure.png saved as PNG format", ""],
            compare.unwrap_stata_file_messages(lines),
        )

    def test_wrapped_saved_as_suffix_is_reassembled(self) -> None:
        lines = [
            "file <OUTPUT>/figure.png saved as",
            "PNG format",
            "",
        ]
        self.assertEqual(
            ["file <OUTPUT>/figure.png saved as PNG format", ""],
            compare.unwrap_stata_file_messages(lines),
        )

    def test_wrapped_saved_then_as_suffix_is_reassembled(self) -> None:
        lines = [
            "file <OUTPUT>/figure.png saved",
            "as PNG format",
            "",
        ]
        self.assertEqual(
            ["file <OUTPUT>/figure.png saved as PNG format", ""],
            compare.unwrap_stata_file_messages(lines),
        )

    def test_post_run_input_hash_is_checked_against_every_manifest(self) -> None:
        roles = ("baseline", "candidate_1", "candidate_2")
        for mismatched_role in roles:
            with (
                self.subTest(role=mismatched_role),
                tempfile.TemporaryDirectory() as tmp,
            ):
                root = Path(tmp)
                input_file = root / "full_db.dta"
                input_file.write_bytes(b"current input")
                current_hash = compare.stata_run.sha256_file(input_file)
                status, report, comparison_count = run_stubbed_comparator(
                    root,
                    input_sha256=current_hash,
                    input_hash_overrides={mismatched_role: "0" * 64},
                    input_file=input_file,
                )
                self.assertEqual(1, status)
                self.assertEqual(
                    [
                        {
                            "category": "post_run_input_hash",
                            "role": mismatched_role,
                        }
                    ],
                    report["failures"],
                )
                self.assertEqual([], report["comparisons"])
                self.assertEqual(0, comparison_count)

    def test_expected_input_hash_remains_an_additional_pin(self) -> None:
        observed_hash = "a" * 64
        failures = compare.validate_post_run_input(
            observed_hash,
            {
                role: {"input": {"sha256": observed_hash}}
                for role in ("baseline", "candidate_1", "candidate_2")
            },
            "b" * 64,
        )
        self.assertEqual(
            [{"category": "post_run_input_hash", "role": "expected"}],
            failures,
        )

    def test_input_file_is_required_even_with_expected_hash(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            report_path = root / "comparison.json"
            with contextlib.redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit) as raised:
                    compare.main(
                        [
                            "--baseline-run",
                            str(root / "baseline"),
                            "--candidate-run-1",
                            str(root / "candidate-1"),
                            "--candidate-run-2",
                            str(root / "candidate-2"),
                            "--expected-input-sha256",
                            "0" * 64,
                            "--report",
                            str(report_path),
                        ]
                    )
            self.assertEqual(2, raised.exception.code)
            self.assertFalse(report_path.exists())

    def test_wrong_input_filename_writes_safe_failure_report(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            report_path = root / "comparison.json"
            restricted_name = root / "private-patient-export.dta"
            with patch.object(
                compare.input_manifest,
                "validate_approved_input",
            ) as validate:
                status = compare.main(
                    [
                        "--baseline-run",
                        str(root / "baseline"),
                        "--candidate-run-1",
                        str(root / "candidate-1"),
                        "--candidate-run-2",
                        str(root / "candidate-2"),
                        "--input-file",
                        str(restricted_name),
                        "--report",
                        str(report_path),
                    ]
                )
            report_text = report_path.read_text(encoding="utf-8")
            report = json.loads(report_text)
            self.assertEqual(1, status)
            self.assertEqual(2, report["schema_version"])
            self.assertEqual("fail", report["status"])
            self.assertEqual([], report["comparisons"])
            self.assertEqual(
                [
                    {
                        "category": "input_manifest_malformed",
                        "role": "current_input",
                    }
                ],
                report["failures"],
            )
            self.assertNotIn(str(restricted_name), report_text)
            self.assertNotIn("private-patient-export", report_text)
            validate.assert_not_called()

    def test_current_input_manifest_failures_write_safe_reports(self) -> None:
        categories = (
            "missing_input",
            "input_manifest_missing",
            "input_manifest_malformed",
            "input_manifest_unapproved",
            "input_manifest_input_mismatch",
            "input_manifest_provenance_mismatch",
            "input_manifest_contract_mismatch",
        )
        for category in categories:
            with self.subTest(category=category), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                input_file = root / "full_db.dta"
                report_path = root / "comparison.json"
                analysis_root = root / "repository-authority"
                secret = f"restricted path {root}/patient-value-123"
                failure = compare.input_manifest.InputManifestFailure(
                    category,
                    secret,
                )
                with patch.object(
                    compare.input_manifest,
                    "validate_approved_input",
                    side_effect=failure,
                ) as validate:
                    status = compare.main(
                        [
                            "--baseline-run",
                            str(root / "baseline"),
                            "--candidate-run-1",
                            str(root / "candidate-1"),
                            "--candidate-run-2",
                            str(root / "candidate-2"),
                            "--input-file",
                            str(input_file),
                            "--analysis-root",
                            str(analysis_root),
                            "--report",
                            str(report_path),
                        ]
                    )
                report_text = report_path.read_text(encoding="utf-8")
                report = json.loads(report_text)
                self.assertEqual(1, status)
                self.assertEqual(2, report["schema_version"])
                self.assertEqual("fail", report["status"])
                self.assertIsNone(report["input_sha256_after_runs"])
                self.assertEqual([], report["comparisons"])
                self.assertEqual(
                    [{"category": category, "role": "current_input"}],
                    report["failures"],
                )
                self.assertNotIn(secret, report_text)
                self.assertNotIn(str(root), report_text)
                self.assertNotIn("patient-value-123", report_text)
                validate.assert_called_once_with(
                    root,
                    analysis_root.resolve(),
                )

    def test_stale_pass_report_is_invalidated_before_run_loading(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            report_path = root / "comparison.json"
            report_path.write_text(
                '{"status":"pass","stale_secret":"must-disappear"}\n',
                encoding="utf-8",
            )
            input_file = root / "full_db.dta"
            current_input = {
                "input_file": input_file,
                "input": {
                    "logical_name": "full_db.dta",
                    "size_bytes": 0,
                    "sha256": "a" * 64,
                },
                "approval": APPROVAL_REFERENCE,
            }
            with patch.object(
                compare.input_manifest,
                "validate_approved_input",
                return_value=current_input,
            ):
                status = compare.main(
                    [
                        "--baseline-run",
                        str(root / "missing-baseline"),
                        "--candidate-run-1",
                        str(root / "missing-candidate-1"),
                        "--candidate-run-2",
                        str(root / "missing-candidate-2"),
                        "--input-file",
                        str(input_file),
                        "--report",
                        str(report_path),
                    ]
                )
            report_text = report_path.read_text(encoding="utf-8")
            report = json.loads(report_text)
            temporary_reports = list(root.glob(".comparison.json.*.tmp"))
        self.assertEqual(1, status)
        self.assertEqual("fail", report["status"])
        self.assertEqual(
            [{"category": "run_manifest_unavailable", "role": "run_set"}],
            report["failures"],
        )
        self.assertNotIn("stale_secret", report_text)
        self.assertEqual([], temporary_reports)

    def test_input_mutation_during_comparison_prevents_pass(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            input_file, analysis_root, approved = write_approved_input_fixture(root)
            input_hash = approved["input"]["sha256"]
            runs = {
                "baseline": root / "baseline",
                "candidate_1": root / "candidate-1",
                "candidate_2": root / "candidate-2",
            }
            write_run_manifest(
                runs["baseline"],
                legacy=True,
                commit="baseline",
                input_sha256=input_hash,
                run_id="baseline-run",
                include_approval=False,
            )
            write_complete_run_evidence(runs["baseline"], legacy=True)
            for role in ("candidate_1", "candidate_2"):
                write_run_manifest(
                    runs[role],
                    legacy=False,
                    commit="candidate",
                    input_sha256=input_hash,
                    run_id=f"{role}-run",
                    approval=approved["approval"],
                )
                write_complete_run_evidence(runs[role], legacy=False)
            report_path = root / "comparison.json"

            def mutate_then_pass(
                _left: Path,
                _right: Path,
                label: str,
            ) -> dict[str, object]:
                input_file.write_bytes(b"mutated during comparison")
                return passing_comparison(_left, _right, label)

            with patch.object(
                compare,
                "compare_pair",
                side_effect=mutate_then_pass,
            ):
                status = compare.main(
                    [
                        "--baseline-run",
                        str(runs["baseline"]),
                        "--candidate-run-1",
                        str(runs["candidate_1"]),
                        "--candidate-run-2",
                        str(runs["candidate_2"]),
                        "--input-file",
                        str(input_file),
                        "--analysis-root",
                        str(analysis_root),
                        "--report",
                        str(report_path),
                    ]
                )
            report = json.loads(report_path.read_text(encoding="utf-8"))
        self.assertEqual(1, status)
        self.assertEqual("fail", report["status"])
        self.assertIn(
            {
                "category": "input_manifest_input_mismatch",
                "role": "current_input",
            },
            report["failures"],
        )

    def test_reapproved_input_during_comparison_is_detected_as_drift(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            input_file, analysis_root, approved = write_approved_input_fixture(root)
            input_hash = approved["input"]["sha256"]
            runs = {
                "baseline": root / "baseline",
                "candidate_1": root / "candidate-1",
                "candidate_2": root / "candidate-2",
            }
            write_run_manifest(
                runs["baseline"],
                legacy=True,
                commit="baseline",
                input_sha256=input_hash,
                run_id="baseline-run",
                include_approval=False,
            )
            write_complete_run_evidence(runs["baseline"], legacy=True)
            for role in ("candidate_1", "candidate_2"):
                write_run_manifest(
                    runs[role],
                    legacy=False,
                    commit="candidate",
                    input_sha256=input_hash,
                    run_id=f"{role}-run",
                    approval=approved["approval"],
                )
                write_complete_run_evidence(runs[role], legacy=False)
            report_path = root / "comparison.json"
            replacement_done = False

            def reapprove_then_pass(
                _left: Path,
                _right: Path,
                label: str,
            ) -> dict[str, object]:
                nonlocal replacement_done
                if not replacement_done:
                    input_file.write_bytes(b"reapproved during comparison")
                    compare.input_manifest.approve_input(
                        input_file.parent,
                        analysis_root,
                        approval="YES",
                        replace="YES",
                    )
                    replacement_done = True
                return passing_comparison(_left, _right, label)

            with patch.object(
                compare,
                "compare_pair",
                side_effect=reapprove_then_pass,
            ):
                status = compare.main(
                    [
                        "--baseline-run",
                        str(runs["baseline"]),
                        "--candidate-run-1",
                        str(runs["candidate_1"]),
                        "--candidate-run-2",
                        str(runs["candidate_2"]),
                        "--input-file",
                        str(input_file),
                        "--analysis-root",
                        str(analysis_root),
                        "--report",
                        str(report_path),
                    ]
                )
            report = json.loads(report_path.read_text(encoding="utf-8"))
        self.assertEqual(1, status)
        self.assertEqual("fail", report["status"])
        self.assertIn(
            {
                "category": "input_manifest_comparison_drift",
                "role": "current_input",
            },
            report["failures"],
        )

    def test_wrong_run_roles_prevent_pair_comparisons(self) -> None:
        cases = (
            ("baseline", False),
            ("candidate_1", True),
            ("candidate_2", True),
        )
        for wrong_role, wrong_mode in cases:
            with self.subTest(role=wrong_role), tempfile.TemporaryDirectory() as tmp:
                status, report, comparison_count = run_stubbed_comparator(
                    Path(tmp),
                    legacy_overrides={wrong_role: wrong_mode},
                )
                self.assertEqual(1, status)
                self.assertIn(
                    {"category": "run_role", "role": wrong_role},
                    report["failures"],
                )
                self.assertEqual([], report["comparisons"])
                self.assertEqual(0, comparison_count)

    def test_filesystem_aliases_fail_run_isolation(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            baseline = root / "baseline"
            candidate = root / "candidate"
            alias = root / "candidate-alias"
            baseline.mkdir()
            candidate.mkdir()
            alias.symlink_to(candidate, target_is_directory=True)
            runs = {
                "baseline": baseline,
                "candidate_1": candidate,
                "candidate_2": alias,
            }
            manifests = {
                "baseline": {
                    "legacy_two_argument_mode": True,
                    "analysis_commit": "baseline",
                    "run_id": "baseline",
                },
                "candidate_1": {
                    "legacy_two_argument_mode": False,
                    "analysis_commit": "candidate",
                    "run_id": "candidate-one",
                },
                "candidate_2": {
                    "legacy_two_argument_mode": False,
                    "analysis_commit": "candidate",
                    "run_id": "candidate-two",
                },
            }
            failures = compare.validate_run_set(runs, manifests, runs)
        self.assertIn(
            {"category": "run_isolation", "role": "run_set"},
            failures,
        )

    def test_artifact_roots_cannot_escape_run_directories(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            run = root / "run"
            shared = root / "shared"
            run.mkdir()
            shared.mkdir()
            (run / "linked").symlink_to(shared, target_is_directory=True)
            for relative in (str(shared), "../shared", "linked"):
                with self.subTest(relative=relative):
                    with self.assertRaises(compare.ArtifactRootFailure):
                        compare.resolve_artifact_root(
                            run,
                            {"artifact_relative_path": relative},
                        )

    def test_artifact_roots_must_be_pairwise_distinct(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            runs = {
                role: root / role
                for role in ("baseline", "candidate_1", "candidate_2")
            }
            for path in runs.values():
                path.mkdir()
            manifests = {
                "baseline": {
                    "legacy_two_argument_mode": True,
                    "analysis_commit": "baseline",
                    "run_id": "baseline",
                },
                "candidate_1": {
                    "legacy_two_argument_mode": False,
                    "analysis_commit": "candidate",
                    "run_id": "candidate-one",
                },
                "candidate_2": {
                    "legacy_two_argument_mode": False,
                    "analysis_commit": "candidate",
                    "run_id": "candidate-two",
                },
            }
            artifact_roots = {
                "baseline": runs["baseline"],
                "candidate_1": runs["candidate_1"],
                "candidate_2": runs["candidate_1"],
            }
            failures = compare.validate_run_set(
                runs,
                manifests,
                artifact_roots,
            )
        self.assertIn(
            {"category": "artifact_isolation", "role": "run_set"},
            failures,
        )

    def test_candidate_run_ids_must_be_nonempty_and_distinct(self) -> None:
        cases = (
            (
                {"candidate_1": "", "candidate_2": "candidate-two"},
                {"category": "run_id", "role": "candidate_1"},
            ),
            (
                {"candidate_1": "candidate-one", "candidate_2": ""},
                {"category": "run_id", "role": "candidate_2"},
            ),
            (
                {"candidate_1": "same-run", "candidate_2": "same-run"},
                {
                    "category": "candidate_run_id_match",
                    "role": "candidate_set",
                },
            ),
        )
        for overrides, expected_failure in cases:
            with (
                self.subTest(expected_failure=expected_failure),
                tempfile.TemporaryDirectory() as tmp,
            ):
                status, report, comparison_count = run_stubbed_comparator(
                    Path(tmp),
                    run_id_overrides=overrides,
                )
                self.assertEqual(1, status)
                self.assertIn(expected_failure, report["failures"])
                self.assertEqual([], report["comparisons"])
                self.assertEqual(0, comparison_count)

    def test_live_guarded_controls_are_required_before_comparison(self) -> None:
        controls = (
            "SUCCESS",
            "run_status.tsv",
            "input_validation.tsv",
            "run_metrics.tsv",
            "ANALYSIS_COMPLETE",
        )
        for control in controls:
            with self.subTest(control=control), tempfile.TemporaryDirectory() as tmp:
                status, report, comparison_count = run_stubbed_comparator(
                    Path(tmp),
                    missing_control=("candidate_1", control),
                )
                self.assertEqual(1, status)
                self.assertIn(
                    {"category": "control_inventory", "role": "candidate_1"},
                    report["failures"],
                )
                self.assertEqual([], report["comparisons"])
                self.assertEqual(0, comparison_count)

    def test_unsafe_artifact_root_writes_path_safe_failure(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            secret = "../restricted-patient-artifacts"
            status, report, comparison_count = run_stubbed_comparator(
                root,
                artifact_relative_overrides={"candidate_1": secret},
            )
            report_text = json.dumps(report)
        self.assertEqual(1, status)
        self.assertEqual(
            [{"category": "artifact_root", "role": "run_set"}],
            report["failures"],
        )
        self.assertEqual([], report["comparisons"])
        self.assertEqual(0, comparison_count)
        self.assertNotIn(secret, report_text)

    def test_candidate_commits_must_be_nonempty_and_match(self) -> None:
        cases = (
            (
                "",
                "candidate-two",
                {"category": "analysis_commit", "role": "candidate_1"},
            ),
            (
                "candidate-one",
                "",
                {"category": "analysis_commit", "role": "candidate_2"},
            ),
            (
                "candidate-one",
                "candidate-two",
                {"category": "candidate_commit_match", "role": "candidate_set"},
            ),
        )
        for first_commit, second_commit, expected_failure in cases:
            with (
                self.subTest(expected_failure=expected_failure),
                tempfile.TemporaryDirectory() as tmp,
            ):
                status, report, comparison_count = run_stubbed_comparator(
                    Path(tmp),
                    commit_overrides={
                        "candidate_1": first_commit,
                        "candidate_2": second_commit,
                    },
                )
                self.assertEqual(1, status)
                self.assertIn(expected_failure, report["failures"])
                self.assertEqual([], report["comparisons"])
                self.assertEqual(0, comparison_count)

    def test_candidates_require_current_input_approval_reference(self) -> None:
        mismatched_approval = {
            **APPROVAL_REFERENCE,
            "manifest_sha256": "9" * 64,
        }
        cases = (
            (
                "candidate_1_missing",
                {"omit_approval_roles": {"candidate_1"}},
                "candidate_1",
            ),
            (
                "candidate_2_missing",
                {"omit_approval_roles": {"candidate_2"}},
                "candidate_2",
            ),
            (
                "candidate_1_mismatch",
                {"approval_overrides": {"candidate_1": mismatched_approval}},
                "candidate_1",
            ),
            (
                "candidate_2_mismatch",
                {"approval_overrides": {"candidate_2": mismatched_approval}},
                "candidate_2",
            ),
        )
        for label, overrides, failed_role in cases:
            with self.subTest(case=label), tempfile.TemporaryDirectory() as tmp:
                status, report, comparison_count = run_stubbed_comparator(
                    Path(tmp),
                    **overrides,
                )
                self.assertEqual(1, status)
                self.assertEqual(
                    [
                        {
                            "category": "run_input_manifest_reference",
                            "role": failed_role,
                        }
                    ],
                    report["failures"],
                )
                self.assertEqual([], report["comparisons"])
                self.assertEqual(0, comparison_count)

    def test_legacy_baseline_may_omit_input_approval_reference(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            status, report, comparison_count = run_stubbed_comparator(
                Path(tmp),
                omit_approval_roles={"baseline"},
            )
            self.assertEqual(0, status)
            self.assertEqual([], report["failures"])
            self.assertEqual(2, comparison_count)

    def test_present_baseline_input_approval_must_match_current(self) -> None:
        mismatched_approval = {
            **APPROVAL_REFERENCE,
            "producer_commit": "8" * 40,
        }
        with tempfile.TemporaryDirectory() as tmp:
            status, report, comparison_count = run_stubbed_comparator(
                Path(tmp),
                approval_overrides={"baseline": mismatched_approval},
            )
            self.assertEqual(1, status)
            self.assertEqual(
                [
                    {
                        "category": "run_input_manifest_reference",
                        "role": "baseline",
                    }
                ],
                report["failures"],
            )
            self.assertEqual([], report["comparisons"])
            self.assertEqual(0, comparison_count)

    def test_input_approval_failure_does_not_expose_reference_values(self) -> None:
        secret = "restricted-producer-secret"
        mismatched_approval = {
            **APPROVAL_REFERENCE,
            "producer_commit": secret,
        }
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            status, report, comparison_count = run_stubbed_comparator(
                root,
                approval_overrides={"candidate_1": mismatched_approval},
            )
            report_text = (root / "comparison.json").read_text(encoding="utf-8")
            self.assertEqual(1, status)
            self.assertNotIn(secret, report_text)
            self.assertEqual(
                {"category", "role"},
                set(report["failures"][0]),
            )
            self.assertEqual(0, comparison_count)

    def test_run_set_failure_details_do_not_expose_commit_values(self) -> None:
        runs = {
            "baseline": Path("/baseline"),
            "candidate_1": Path("/candidate-one"),
            "candidate_2": Path("/candidate-two"),
        }
        failures = compare.validate_run_set(
            runs,
            {
                "baseline": {
                    "legacy_two_argument_mode": False,
                    "analysis_commit": "baseline-secret",
                    "run_id": "baseline-run",
                },
                "candidate_1": {
                    "legacy_two_argument_mode": False,
                    "analysis_commit": "candidate-secret-one",
                    "run_id": "candidate-one-run",
                },
                "candidate_2": {
                    "legacy_two_argument_mode": False,
                    "analysis_commit": "candidate-secret-two",
                    "run_id": "candidate-two-run",
                },
            },
            runs,
        )
        self.assertTrue(failures)
        self.assertTrue(
            all(set(failure) == {"category", "role"} for failure in failures)
        )
        self.assertNotIn("secret", repr(failures))

    def test_valid_run_set_reaches_labeled_pair_comparisons(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            input_file = root / "full_db.dta"
            input_file.write_bytes(b"current input")
            current_hash = compare.stata_run.sha256_file(input_file)
            status, report, comparison_count = run_stubbed_comparator(
                root,
                input_sha256=current_hash,
                input_file=input_file,
            )
            self.assertEqual(0, status)
            self.assertEqual(2, report["schema_version"])
            self.assertEqual(
                [
                    "baseline_vs_candidate_1",
                    "candidate_repeatability",
                ],
                [item["comparison"] for item in report["comparisons"]],
            )
            self.assertEqual(2, comparison_count)

    def test_correction_mode_passes_with_historical_impact_and_repeatability(
        self,
    ) -> None:
        def correction_comparison(
            _left: Path,
            _right: Path,
            label: str,
        ) -> dict[str, object]:
            if label == "historical_impact":
                return {
                    "comparison": label,
                    "status": "fail",
                    "failures": [
                        {
                            "category": "workbook_value",
                            "artifact": "Def8-Summary.xlsx",
                            "location": "Results!A2",
                        }
                    ],
                }
            return passing_comparison(_left, _right, label)

        with tempfile.TemporaryDirectory() as tmp:
            status, report, comparison_count = run_stubbed_comparator(
                Path(tmp),
                comparison_mode="correction",
                comparison_side_effect=correction_comparison,
            )
        self.assertEqual(0, status)
        self.assertEqual("correction", report["comparison_mode"])
        self.assertEqual("pass", report["status"])
        self.assertEqual("changed", report["comparisons"][0]["status"])
        self.assertEqual("pass", report["comparisons"][1]["status"])
        self.assertEqual(2, comparison_count)
        self.assertNotIn("1.25", repr(report))

    def test_correction_mode_still_requires_candidate_repeatability(self) -> None:
        def nonrepeatable_comparison(
            _left: Path,
            _right: Path,
            label: str,
        ) -> dict[str, object]:
            if label == "candidate_repeatability":
                return {
                    "comparison": label,
                    "status": "fail",
                    "failures": [
                        {
                            "category": "png_pixels",
                            "artifact": "Definition Overlap HeatPlot.png",
                        }
                    ],
                }
            return passing_comparison(_left, _right, label)

        with tempfile.TemporaryDirectory() as tmp:
            status, report, comparison_count = run_stubbed_comparator(
                Path(tmp),
                comparison_mode="correction",
                comparison_side_effect=nonrepeatable_comparison,
            )
        self.assertEqual(1, status)
        self.assertEqual("fail", report["status"])
        self.assertEqual(2, comparison_count)

    def test_correction_report_cannot_pass_without_historical_impact(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            report = compare.write_report(
                Path(tmp) / "comparison.json",
                observed_input_hash="a" * 64,
                comparisons=[
                    {
                        "comparison": "candidate_repeatability",
                        "status": "pass",
                        "failures": [],
                    }
                ],
                failures=[],
                comparison_mode="correction",
                announce=False,
            )
        self.assertEqual("fail", report["status"])

    def test_correction_mode_preserves_evidence_integrity_failures(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            status, report, comparison_count = run_stubbed_comparator(
                Path(tmp),
                comparison_mode="correction",
                legacy_overrides={"candidate_1": True},
            )
        self.assertEqual(1, status)
        self.assertIn(
            {"category": "run_role", "role": "candidate_1"},
            report["failures"],
        )
        self.assertEqual([], report["comparisons"])
        self.assertEqual(0, comparison_count)

    def test_correction_mode_rejects_baseline_log_root_integrity_failure(
        self,
    ) -> None:
        def baseline_root_failure(
            _left: Path,
            _right: Path,
            label: str,
        ) -> dict[str, object]:
            if label == "historical_impact":
                return {
                    "comparison": label,
                    "status": "fail",
                    "failures": [
                        {
                            "category": "analysis_log_output_roots",
                            "artifact": "analysis_log",
                            "role": "left",
                        }
                    ],
                }
            return passing_comparison(_left, _right, label)

        with tempfile.TemporaryDirectory() as tmp:
            status, report, comparison_count = run_stubbed_comparator(
                Path(tmp),
                comparison_mode="correction",
                comparison_side_effect=baseline_root_failure,
            )
        self.assertEqual(1, status)
        self.assertEqual("fail", report["status"])
        self.assertIn(
            {
                "category": "analysis_log_output_roots",
                "role": "baseline",
            },
            report["failures"],
        )
        self.assertEqual("changed", report["comparisons"][0]["status"])
        self.assertEqual("pass", report["comparisons"][1]["status"])
        self.assertEqual(2, comparison_count)


if __name__ == "__main__":
    unittest.main()
