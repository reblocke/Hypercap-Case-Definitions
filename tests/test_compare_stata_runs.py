from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

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


if __name__ == "__main__":
    unittest.main()
