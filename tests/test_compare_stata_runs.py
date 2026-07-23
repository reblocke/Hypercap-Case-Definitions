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
) -> None:
    path.mkdir()
    (path / "run_manifest.json").write_text(
        json.dumps(
            {
                "legacy_two_argument_mode": legacy,
                "analysis_commit": commit,
                "input": {"sha256": input_sha256},
            }
        ),
        encoding="utf-8",
    )


def passing_comparison(
    _left: Path,
    _right: Path,
    label: str,
) -> dict[str, object]:
    return {"comparison": label, "status": "pass", "failures": []}


def run_stubbed_comparator(
    root: Path,
    *,
    legacy_overrides: dict[str, bool] | None = None,
    commit_overrides: dict[str, str] | None = None,
    input_sha256: str | None = None,
    input_hash_overrides: dict[str, str] | None = None,
    input_file: Path | None = None,
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

    runs = {role: root / role for role in roles}
    for role, path in runs.items():
        write_run_manifest(
            path,
            legacy=legacy_modes[role],
            commit=commits[role],
            input_sha256=input_hashes[role],
        )

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
    ]
    with patch.object(
        compare,
        "compare_pair",
        side_effect=passing_comparison,
    ) as compare_pair:
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

    def test_run_set_failure_details_do_not_expose_commit_values(self) -> None:
        failures = compare.validate_run_set(
            {
                "baseline": Path("/baseline"),
                "candidate_1": Path("/candidate-one"),
                "candidate_2": Path("/candidate-two"),
            },
            {
                "baseline": {
                    "legacy_two_argument_mode": False,
                    "analysis_commit": "baseline-secret",
                },
                "candidate_1": {
                    "legacy_two_argument_mode": False,
                    "analysis_commit": "candidate-secret-one",
                },
                "candidate_2": {
                    "legacy_two_argument_mode": False,
                    "analysis_commit": "candidate-secret-two",
                },
            },
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
            self.assertEqual(
                [
                    "baseline_vs_candidate_1",
                    "candidate_repeatability",
                ],
                [item["comparison"] for item in report["comparisons"]],
            )
            self.assertEqual(2, comparison_count)


if __name__ == "__main__":
    unittest.main()
