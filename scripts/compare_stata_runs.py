#!/usr/bin/env python3
"""Compare one legacy baseline and two guarded Stata runs without exposing values."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import numbers
import re
import sys
from pathlib import Path
from typing import Any, Iterable

from openpyxl import load_workbook
from PIL import Image

import stata_run

ABS_TOL = 1e-12
REL_TOL = 1e-12


def digest_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def load_run(path: Path) -> tuple[dict[str, Any], Path]:
    manifest_path = path / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    relative = manifest.get("artifact_relative_path", ".")
    artifact_root = path if relative in {"", "."} else path / relative
    return manifest, artifact_root


def compare_scalar(
    left: Any,
    right: Any,
    *,
    artifact: str,
    location: str,
    failures: list[dict[str, Any]],
) -> None:
    if (
        isinstance(left, numbers.Real)
        and not isinstance(left, bool)
        and isinstance(right, numbers.Real)
        and not isinstance(right, bool)
    ):
        if math.isclose(float(left), float(right), rel_tol=REL_TOL, abs_tol=ABS_TOL):
            return
    elif left == right:
        return
    failures.append(
        {
            "category": "workbook_value",
            "artifact": artifact,
            "location": location,
        }
    )


def compare_workbook(
    left_path: Path,
    right_path: Path,
    failures: list[dict[str, Any]],
) -> None:
    artifact = left_path.name
    left_book = load_workbook(left_path, data_only=False, read_only=False)
    right_book = load_workbook(right_path, data_only=False, read_only=False)
    try:
        if left_book.sheetnames != right_book.sheetnames:
            failures.append({"category": "workbook_sheets", "artifact": artifact})
            return
        for sheet_name in left_book.sheetnames:
            left_sheet = left_book[sheet_name]
            right_sheet = right_book[sheet_name]
            if (left_sheet.max_row, left_sheet.max_column) != (
                right_sheet.max_row,
                right_sheet.max_column,
            ):
                failures.append(
                    {
                        "category": "workbook_dimensions",
                        "artifact": artifact,
                        "location": sheet_name,
                    }
                )
                continue
            left_merges = sorted(str(item) for item in left_sheet.merged_cells.ranges)
            right_merges = sorted(str(item) for item in right_sheet.merged_cells.ranges)
            if left_merges != right_merges:
                failures.append(
                    {
                        "category": "workbook_merges",
                        "artifact": artifact,
                        "location": sheet_name,
                    }
                )
            for row in left_sheet.iter_rows():
                for left_cell in row:
                    right_cell = right_sheet[left_cell.coordinate]
                    compare_scalar(
                        left_cell.value,
                        right_cell.value,
                        artifact=artifact,
                        location=f"{sheet_name}!{left_cell.coordinate}",
                        failures=failures,
                    )
                    if (
                        left_cell.data_type != right_cell.data_type
                        or left_cell.number_format != right_cell.number_format
                    ):
                        failures.append(
                            {
                                "category": "workbook_cell_metadata",
                                "artifact": artifact,
                                "location": f"{sheet_name}!{left_cell.coordinate}",
                            }
                        )
    finally:
        left_book.close()
        right_book.close()


def compare_png(
    left_path: Path,
    right_path: Path,
    failures: list[dict[str, Any]],
) -> None:
    with Image.open(left_path) as left_image, Image.open(right_path) as right_image:
        left = left_image.convert("RGBA")
        right = right_image.convert("RGBA")
        if left.size != right.size:
            failures.append(
                {"category": "png_dimensions", "artifact": left_path.name}
            )
            return
        if left.tobytes() != right.tobytes():
            failures.append({"category": "png_pixels", "artifact": left_path.name})


def strip_row_listing(lines: list[str]) -> list[str]:
    output: list[str] = []
    skipping = False
    for line in lines:
        if "list hypercap_resp_failure paco2 vbg_co2" in line:
            skipping = True
            continue
        if skipping and "/* Ones I can't do */" in line:
            skipping = False
            output.append(line)
            continue
        if not skipping:
            output.append(line)
    return output


def normalize_log(path: Path, artifact_root: Path, run_root: Path) -> str:
    text = path.read_text(encoding="utf-8", errors="replace")
    start_match = re.search(r"(?m)^.*Pre-processing\s*$", text)
    if not start_match:
        raise ValueError("preprocessing marker unavailable")
    start = text.rfind("/*", 0, start_match.start())
    if start < 0:
        start = start_match.start()

    epilogue = text.find("/* Aggregate-only run metrics", start)
    end = epilogue if epilogue >= 0 else len(text)
    excerpt = text[start:end]
    for root in (artifact_root, run_root, artifact_root.parent):
        excerpt = excerpt.replace(root.as_posix(), "<OUTPUT>")
        excerpt = excerpt.replace(str(root), "<OUTPUT>")
    mac_user_prefix = "/" + "Users" + "/"
    excerpt = re.sub(
        re.escape(mac_user_prefix) + r"[^\n\"]+/full_db[.]dta",
        "<INPUT>/full_db.dta",
        excerpt,
    )
    lines = strip_row_listing(excerpt.splitlines())
    normalized: list[str] = []
    previous_blank = False
    for line in lines:
        clean = line.rstrip()
        blank = not clean
        if blank and previous_blank:
            continue
        normalized.append(clean)
        previous_blank = blank
    return "\n".join(normalized).strip() + "\n"


def compare_logs(
    left_manifest: dict[str, Any],
    left_artifacts: Path,
    left_run: Path,
    right_manifest: dict[str, Any],
    right_artifacts: Path,
    right_run: Path,
    failures: list[dict[str, Any]],
) -> None:
    left_log = left_artifacts / left_manifest["artifacts"]["analysis_log"]
    right_log = right_artifacts / right_manifest["artifacts"]["analysis_log"]
    left_text = normalize_log(left_log, left_artifacts, left_run)
    right_text = normalize_log(right_log, right_artifacts, right_run)
    if left_text == right_text:
        return
    left_lines = left_text.splitlines()
    right_lines = right_text.splitlines()
    differing_line = min(len(left_lines), len(right_lines)) + 1
    for index, (left, right) in enumerate(zip(left_lines, right_lines), start=1):
        if left != right:
            differing_line = index
            break
    failures.append(
        {
            "category": "analysis_log",
            "artifact": "analysis_log",
            "normalized_line": differing_line,
            "left_sha256": digest_text(left_text),
            "right_sha256": digest_text(right_text),
        }
    )


def compare_environment(
    left: dict[str, Any],
    right: dict[str, Any],
    failures: list[dict[str, Any]],
) -> None:
    if left["input"]["sha256"] != right["input"]["sha256"]:
        failures.append({"category": "input_hash", "artifact": "run_manifest"})
    left_status = left.get("driver_status", {})
    right_status = right.get("driver_status", {})
    for key in ("stata_version", "stata_flavor", "operating_system", "machine_type"):
        if left_status.get(key) != right_status.get(key):
            failures.append(
                {
                    "category": "stata_environment",
                    "artifact": "run_manifest",
                    "field": key,
                }
            )
    left_dependencies = {
        row["name"]: (row["found"], row["sha256"], row["version_header"])
        for row in left.get("dependencies", [])
    }
    right_dependencies = {
        row["name"]: (row["found"], row["sha256"], row["version_header"])
        for row in right.get("dependencies", [])
    }
    if left_dependencies != right_dependencies:
        failures.append({"category": "dependencies", "artifact": "dependency_report"})


def compare_pair(
    left_run: Path,
    right_run: Path,
    label: str,
) -> dict[str, Any]:
    left_manifest, left_artifacts = load_run(left_run)
    right_manifest, right_artifacts = load_run(right_run)
    failures: list[dict[str, Any]] = []

    compare_environment(left_manifest, right_manifest, failures)
    for manifest in (left_manifest, right_manifest):
        if manifest.get("status") != "success":
            failures.append({"category": "run_status", "artifact": "run_manifest"})
        if manifest.get("analysis_worktree_dirty"):
            failures.append({"category": "dirty_worktree", "artifact": "run_manifest"})
        if manifest.get("artifacts", {}).get("missing"):
            failures.append({"category": "artifact_inventory", "artifact": "run_manifest"})
        if manifest.get("artifacts", {}).get("control_missing"):
            failures.append({"category": "control_inventory", "artifact": "run_manifest"})

    for name in stata_run.EXPECTED_XLSX:
        compare_workbook(left_artifacts / name, right_artifacts / name, failures)
    for name in stata_run.EXPECTED_PNG:
        compare_png(left_artifacts / name, right_artifacts / name, failures)
    for name in stata_run.EXPECTED_GPH:
        for root in (left_artifacts, right_artifacts):
            path = root / "graph-temp" / name
            if not path.is_file() or path.stat().st_size == 0:
                failures.append({"category": "gph_presence", "artifact": name})

    for manifest, artifact_root in (
        (left_manifest, left_artifacts),
        (right_manifest, right_artifacts),
    ):
        copied = artifact_root / manifest["artifacts"]["copied_do"]
        expected = manifest["code_hashes"][f"analysis:{stata_run.MAIN_DO}"]
        if stata_run.sha256_file(copied) != expected:
            failures.append({"category": "copied_do_hash", "artifact": copied.name})

    compare_logs(
        left_manifest,
        left_artifacts,
        left_run,
        right_manifest,
        right_artifacts,
        right_run,
        failures,
    )
    return {
        "comparison": label,
        "status": "pass" if not failures else "fail",
        "failures": failures,
    }


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--baseline-run", required=True)
    result.add_argument("--candidate-run-1", required=True)
    result.add_argument("--candidate-run-2", required=True)
    result.add_argument("--input-file")
    result.add_argument("--expected-input-sha256")
    result.add_argument("--report", required=True)
    return result


def main(argv: Iterable[str] | None = None) -> int:
    args = parser().parse_args(argv)
    baseline = Path(args.baseline_run).resolve()
    candidate_one = Path(args.candidate_run_1).resolve()
    candidate_two = Path(args.candidate_run_2).resolve()
    comparisons = [
        compare_pair(baseline, candidate_one, "baseline_vs_candidate_1"),
        compare_pair(candidate_one, candidate_two, "candidate_repeatability"),
    ]
    failures: list[dict[str, Any]] = []

    if len({baseline, candidate_one, candidate_two}) != 3:
        failures.append({"category": "run_isolation", "artifact": "run_directories"})

    observed_input_hash = ""
    if args.input_file:
        observed_input_hash = stata_run.sha256_file(Path(args.input_file))
        if (
            args.expected_input_sha256
            and observed_input_hash.lower() != args.expected_input_sha256.lower()
        ):
            failures.append({"category": "post_run_input_hash", "artifact": "input"})

    report = {
        "schema_version": 1,
        "status": "pass"
        if not failures and all(item["status"] == "pass" for item in comparisons)
        else "fail",
        "input_sha256_after_runs": observed_input_hash,
        "comparisons": comparisons,
        "failures": failures,
    }
    report_path = Path(args.report)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(f"Stata run comparison: {report['status'].upper()}")
    return 0 if report["status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
