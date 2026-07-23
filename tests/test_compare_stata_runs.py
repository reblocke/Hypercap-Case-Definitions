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


if __name__ == "__main__":
    unittest.main()
