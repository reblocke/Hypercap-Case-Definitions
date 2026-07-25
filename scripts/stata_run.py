#!/usr/bin/env python3
"""Run the restricted Stata analysis in an isolated, provenance-bearing folder."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import input_manifest

ROOT = Path(__file__).resolve().parents[1]
MAIN_DO = "Hypercapnia Case Definitions.do"
RUN_ID_PATTERN = re.compile(r"^[A-Za-z0-9._-]+$")

EXPECTED_XLSX = (
    "Overall Cohort chars.xlsx",
    *(f"Def{index}-Summary.xlsx" for index in range(1, 11)),
    "Case Definition by Workup.xlsx",
    "Location by Case Definitions.xlsx",
)

EXPECTED_PNG = (
    "Definition Overlap HeatPlot.png",
    "Definition Overlap HeatPlot - Kappa.png",
    "Definition Overlap HeatPlot - Agreement.png",
    "Definition Overlap HeatPlot - PABAK.png",
    "ABG-only Definition Overlap HeatPlot - Kappa.png",
    "VBG-only Definition Overlap HeatPlot - Kappa.png",
    "Both ABG and VBG Definition Overlap HeatPlot - Kappa.png",
    *(f"Loc{index}-Definition Overlap HeatPlot - Kappa.png" for index in range(4)),
    "Unadjusted Prob of Dx Hypercapnia Splines .png",
    "Figure 3 Prob Hypercap ICD.png",
    "South - Unadjusted Prob of Dx Hypercapnia Splines .png",
    "Northeast - Unadjusted Prob of Dx Hypercapnia Splines .png",
    "Midwest - Unadjusted Prob of Dx Hypercapnia Splines .png",
    "West - Unadjusted Prob of Dx Hypercapnia Splines.png",
    "Location - e-Figure 5 Prob Hypercap ICD.png",
)

LEGACY_PNG_ALIASES = {
    "Both ABG and VBG Definition Overlap HeatPlot - Kappa.png": (
        "ABG-VBG Definition Overlap HeatPlot - Kappa.png"
    ),
    "Figure 3 Prob Hypercap ICD.png": "Figure 2 Prob Hypercap ICD.png",
    "Location - e-Figure 5 Prob Hypercap ICD.png": (
        "Location - Figure S3 Prob Hypercap ICD.png"
    ),
}

EXPECTED_GPH = (
    "All_Encounters_Prob_Dx_spline.gph",
    "Emer_Encounters_Prob_Dx_spline.gph",
    "Inp_Encounters_Prob_Dx_spline.gph",
    *(f"Loc{index}_Encounters_Prob_Dx_spline.gph" for index in range(4)),
)

HASHED_RUN_FILES = (
    Path(MAIN_DO),
    Path("data_dictionary.csv"),
    Path("metadata/output_manifest.csv"),
    Path("metadata/stata_dependencies.csv"),
    Path("metadata/upstream_dependency.yml"),
)

HASHED_HARNESS_FILES = (
    Path("scripts/input_manifest.py"),
    Path("scripts/run_stata.sh"),
    Path("scripts/stata_run.py"),
    Path("stata/run_hypercapnia.do"),
    Path("stata/preflight_dependencies.do"),
    Path("stata/validate_input.do"),
)


class RunFailure(RuntimeError):
    """Expected runner failure with a stable machine-readable category."""

    def __init__(self, category: str, message: str):
        super().__init__(message)
        self.category = category


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_json(path: Path, payload: dict[str, Any]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, path)


def git_state(root: Path) -> tuple[str, bool]:
    commit = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=root,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    status = subprocess.run(
        ["git", "status", "--porcelain"],
        cwd=root,
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    return commit, bool(status.strip())


def resolve_stata(value: str | None) -> Path:
    selected = value or os.environ.get("STATA_BIN")
    if selected:
        path = Path(selected).expanduser()
        if path.is_file() and os.access(path, os.X_OK):
            return path.resolve()
        raise RunFailure("missing_stata", "Configured Stata executable is unavailable.")

    for command in ("stata-mp", "stata-se", "stata", "StataMP", "StataSE", "StataBE"):
        located = shutil.which(command)
        if located:
            return Path(located).resolve()

    application_candidates = (
        Path("/Applications/Stata/StataMP.app/Contents/MacOS/StataMP"),
        Path("/Applications/Stata/StataSE.app/Contents/MacOS/StataSE"),
        Path("/Applications/Stata/StataBE.app/Contents/MacOS/StataBE"),
    )
    for path in application_candidates:
        if path.is_file() and os.access(path, os.X_OK):
            return path
    raise RunFailure("missing_stata", "No supported Stata executable was found.")


def stata_mode(path: Path, requested: str) -> tuple[str, str]:
    if requested == "macos-e":
        return "-e", "macos-e"
    if requested == "console-b":
        return "-b", "console-b"
    if ".app/Contents/MacOS/" in path.as_posix():
        return "-e", "macos-e"
    return "-b", "console-b"


def stata_cli_quote(value: str | Path) -> str:
    text = str(value)
    if '"' in text:
        raise RunFailure("unsafe_path", "Stata arguments may not contain double quotes.")
    return f'"{text}"'


def default_run_id(commit: str) -> str:
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    return f"{timestamp}-{commit[:8]}"


def hash_contract_files(analysis_root: Path) -> dict[str, str]:
    hashes: dict[str, str] = {}
    for relative in HASHED_RUN_FILES:
        path = analysis_root / relative
        if path.is_file():
            hashes[f"analysis:{relative.as_posix()}"] = sha256_file(path)
    for relative in HASHED_HARNESS_FILES:
        path = ROOT / relative
        if path.is_file():
            hashes[f"harness:{relative.as_posix()}"] = sha256_file(path)
    return hashes


def parse_key_value_tsv(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    if not path.is_file():
        return values
    with path.open(newline="", encoding="utf-8", errors="replace") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) >= 2:
                values[row[0]] = row[1]
    return values


def dependency_report(raw_path: Path, clean_path: Path) -> list[dict[str, Any]]:
    dependencies: list[dict[str, Any]] = []
    if not raw_path.is_file():
        return dependencies

    with raw_path.open(newline="", encoding="utf-8", errors="replace") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    for row in rows:
        resolved = Path(row.get("resolved_path", "")) if row.get("resolved_path") else None
        file_hash = ""
        version_header = ""
        if resolved and resolved.is_file():
            file_hash = sha256_file(resolved)
            if resolved.suffix.lower() in {".ado", ".scheme"}:
                with resolved.open(encoding="utf-8", errors="replace") as source:
                    for _ in range(5):
                        line = source.readline().strip()
                        if not line:
                            continue
                        if line.startswith("*!"):
                            version_header = line[2:].strip()
                            break
        dependencies.append(
            {
                "name": row.get("name", ""),
                "kind": row.get("kind", ""),
                "required": row.get("required", "") == "true",
                "found": row.get("found", "") == "true",
                "sha256": file_hash,
                "version_header": version_header,
            }
        )

    with clean_path.open("w", newline="", encoding="utf-8") as handle:
        fieldnames = (
            "name",
            "kind",
            "required",
            "found",
            "sha256",
            "version_header",
        )
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(dependencies)
    raw_path.unlink()
    return dependencies


def required_artifact_paths(root: Path, *, legacy: bool = False) -> list[Path]:
    paths = [root / name for name in EXPECTED_XLSX]
    png_names = (
        tuple(LEGACY_PNG_ALIASES.get(name, name) for name in EXPECTED_PNG)
        if legacy
        else EXPECTED_PNG
    )
    paths.extend(root / name for name in png_names)
    paths.extend(root / "graph-temp" / name for name in EXPECTED_GPH)
    return paths


def validate_artifact_inventory(root: Path, *, legacy: bool = False) -> dict[str, Any]:
    expected = required_artifact_paths(root, legacy=legacy)
    missing = [
        path.relative_to(root).as_posix()
        for path in expected
        if not path.is_file() or path.stat().st_size == 0
    ]
    logs = sorted((root / "Logs").glob("*.log")) if (root / "Logs").is_dir() else []
    scripts = sorted((root / "Logs").glob("*.do")) if (root / "Logs").is_dir() else []
    if len(logs) != 1:
        missing.append(f"Logs/*.log (expected 1, found {len(logs)})")
    if len(scripts) != 1:
        missing.append(f"Logs/*.do (expected 1, found {len(scripts)})")
    return {
        "expected_legacy_count": 40,
        "found_legacy_count": len(expected) - len(
            [item for item in missing if not item.startswith("Logs/")]
        )
        + len(logs)
        + len(scripts),
        "missing": missing,
        "analysis_log": logs[0].relative_to(root).as_posix() if len(logs) == 1 else "",
        "copied_do": scripts[0].relative_to(root).as_posix()
        if len(scripts) == 1
        else "",
    }


def validate_control_artifacts(
    run_dir: Path,
    artifact_root: Path,
    inventory: dict[str, Any],
    *,
    legacy: bool,
    include_success: bool,
) -> dict[str, Any]:
    paths = [
        run_dir / "run_manifest.json",
        run_dir / "dependency_report.tsv",
        run_dir / "run_status.tsv",
        artifact_root / inventory.get("analysis_log", ""),
        artifact_root / inventory.get("copied_do", ""),
    ]
    if not legacy:
        paths.extend(
            [
                artifact_root / "input_validation.tsv",
                artifact_root / "run_metrics.tsv",
                artifact_root / "ANALYSIS_COMPLETE",
            ]
        )
    if include_success:
        paths.append(run_dir / "SUCCESS")

    missing = []
    for path in paths:
        if not path.is_file() or path.stat().st_size == 0:
            try:
                label = path.relative_to(run_dir).as_posix()
            except ValueError:
                label = path.name
            missing.append(label)
    return {
        "control_expected_count": len(paths),
        "control_found_count": len(paths) - len(missing),
        "control_missing": missing,
    }


def locate_legacy_artifact_root(run_dir: Path) -> Path:
    legacy_root = run_dir / "legacy-output"
    if not legacy_root.is_dir():
        raise RunFailure(
            "legacy_output_layout",
            "Expected legacy output directory is unavailable.",
        )
    candidates = sorted(path for path in legacy_root.iterdir() if path.is_dir())
    if len(candidates) != 1:
        raise RunFailure(
            "legacy_output_layout",
            f"Expected one legacy date directory, found {len(candidates)}.",
        )
    return candidates[0]


def initial_manifest(
    *,
    run_id: str,
    commit: str,
    dirty: bool,
    approved_input: dict[str, Any],
    mode: str,
    stata: Path,
    legacy: bool,
    code_hashes: dict[str, str],
) -> dict[str, Any]:
    input_record = approved_input["input"]
    approval = approved_input["approval"]
    return {
        "schema_version": 2,
        "run_id": run_id,
        "status": "running",
        "started_at_utc": utc_now(),
        "completed_at_utc": None,
        "analysis_commit": commit,
        "analysis_worktree_dirty": dirty,
        "validation_eligible": not dirty,
        "legacy_two_argument_mode": legacy,
        "input": {
            "logical_name": input_record["logical_name"],
            "size_bytes": input_record["size_bytes"],
            "sha256": input_record["sha256"],
            "approval": {
                "manifest_sha256": approval["manifest_sha256"],
                "upstream_repository": approval["upstream_repository"],
                "producer_commit": approval["producer_commit"],
                "input_schema_version": approval["input_schema_version"],
                "data_dictionary_sha256": approval["data_dictionary_sha256"],
            },
        },
        "stata": {
            "executable_name": stata.name,
            "invocation_mode": mode,
        },
        "code_hashes": code_hashes,
        "dependencies": [],
        "driver_status": {},
        "artifacts": {},
        "metrics": {},
    }


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--input-root", default="data/private")
    result.add_argument("--output-root", default="outputs/stata")
    result.add_argument("--analysis-root", default=str(ROOT))
    result.add_argument("--stata-bin")
    result.add_argument(
        "--stata-mode",
        choices=("auto", "macos-e", "console-b"),
        default="auto",
    )
    result.add_argument("--run-id")
    result.add_argument("--expected-input-sha256")
    result.add_argument("--legacy-two-arg", action="store_true")
    return result


def run(args: argparse.Namespace) -> tuple[int, Path | None]:
    analysis_root = Path(args.analysis_root).expanduser().resolve()
    input_root = Path(args.input_root).expanduser().resolve()
    output_root = Path(args.output_root).expanduser().resolve()
    try:
        same_checkout = analysis_root.samefile(ROOT)
    except OSError:
        same_checkout = False
    if not same_checkout:
        raise RunFailure(
            "analysis_root_mismatch",
            "Analysis root must be the runner's repository checkout.",
        )
    commit, dirty = git_state(analysis_root)
    run_id = args.run_id or default_run_id(commit)
    if not RUN_ID_PATTERN.fullmatch(run_id):
        raise RunFailure("invalid_run_id", "Run ID contains unsafe characters.")

    output_root.mkdir(parents=True, exist_ok=True)
    run_dir = output_root / run_id
    try:
        run_dir.mkdir()
    except FileExistsError as exc:
        raise RunFailure("run_collision", "Run directory already exists.") from exc

    manifest_path = run_dir / "run_manifest.json"
    manifest: dict[str, Any] = {
        "schema_version": 2,
        "run_id": run_id,
        "status": "initializing",
        "started_at_utc": utc_now(),
    }
    atomic_json(manifest_path, manifest)

    try:
        try:
            approved_input = input_manifest.validate_approved_input(
                input_root,
                analysis_root,
            )
        except input_manifest.InputManifestFailure as exc:
            raise RunFailure(exc.category, str(exc)) from exc

        input_hash = approved_input["input"]["sha256"]
        if (
            args.expected_input_sha256
            and input_hash.lower() != args.expected_input_sha256.lower()
        ):
            raise RunFailure("input_hash_mismatch", "Input SHA-256 did not match.")

        stata = resolve_stata(args.stata_bin)
        flag, mode = stata_mode(stata, args.stata_mode)

        analysis_do = analysis_root / MAIN_DO
        if not analysis_do.is_file():
            raise RunFailure("missing_analysis", "Main analysis do-file is unavailable.")

        manifest = initial_manifest(
            run_id=run_id,
            commit=commit,
            dirty=dirty,
            approved_input=approved_input,
            mode=mode,
            stata=stata,
            legacy=args.legacy_two_arg,
            code_hashes=hash_contract_files(analysis_root),
        )
        atomic_json(manifest_path, manifest)

        status_path = run_dir / "run_status.tsv"
        raw_dependencies = run_dir / "dependency_report.raw.tsv"
        if args.legacy_two_arg:
            target_output_root = run_dir / "legacy-output"
            target_output_root.mkdir()
        else:
            target_output_root = output_root

        command = [
            str(stata),
            flag,
            "do",
            stata_cli_quote(ROOT / "stata/run_hypercapnia.do"),
            stata_cli_quote(ROOT),
            stata_cli_quote(analysis_root),
            stata_cli_quote(input_root),
            stata_cli_quote(target_output_root),
            stata_cli_quote(run_id),
            stata_cli_quote(status_path),
            stata_cli_quote(raw_dependencies),
            stata_cli_quote("legacy" if args.legacy_two_arg else "candidate"),
        ]
        process = subprocess.run(command, cwd=run_dir, check=False)
        manifest["stata"]["process_return_code"] = process.returncode
        if process.returncode != 0:
            raise RunFailure(
                "stata_process_exit",
                "The Stata process exited abnormally.",
            )

        status = parse_key_value_tsv(status_path)
        manifest["driver_status"] = status
        manifest["dependencies"] = dependency_report(
            raw_dependencies,
            run_dir / "dependency_report.tsv",
        )
        if status.get("status") != "success":
            raise RunFailure(
                status.get("status", "missing_status"),
                "Stata did not produce a successful fresh status artifact.",
            )

        artifact_root = (
            locate_legacy_artifact_root(run_dir)
            if args.legacy_two_arg
            else run_dir
        )
        inventory = validate_artifact_inventory(artifact_root)
        inventory.update(
            validate_control_artifacts(
                run_dir,
                artifact_root,
                inventory,
                legacy=args.legacy_two_arg,
                include_success=False,
            )
        )
        manifest["artifacts"] = inventory
        manifest["artifact_relative_path"] = artifact_root.relative_to(run_dir).as_posix()
        if inventory["missing"] or inventory["control_missing"]:
            raise RunFailure("incomplete_artifacts", "Required artifacts are incomplete.")

        metrics_path = artifact_root / "run_metrics.tsv"
        if metrics_path.is_file():
            manifest["metrics"] = parse_key_value_tsv(metrics_path)

        try:
            final_approved_input = input_manifest.validate_approved_input(
                input_root,
                analysis_root,
            )
        except input_manifest.InputManifestFailure as exc:
            raise RunFailure(exc.category, str(exc)) from exc
        if (
            final_approved_input["input"] != approved_input["input"]
            or final_approved_input["approval"] != approved_input["approval"]
        ):
            raise RunFailure(
                "input_manifest_run_drift",
                "The approved input changed during the Stata run.",
            )

        success_path = run_dir / "SUCCESS"
        success_path.write_text("status=success\n", encoding="utf-8")
        inventory.update(
            validate_control_artifacts(
                run_dir,
                artifact_root,
                inventory,
                legacy=args.legacy_two_arg,
                include_success=True,
            )
        )
        if inventory["control_missing"]:
            raise RunFailure("incomplete_artifacts", "Run success marker is incomplete.")
        manifest["status"] = "success"
        manifest["completed_at_utc"] = utc_now()
        atomic_json(manifest_path, manifest)
        print(f"Run {run_id}: SUCCESS")
        print(f"Run directory: {run_dir}")
        return 0, run_dir
    except RunFailure as exc:
        (run_dir / "SUCCESS").unlink(missing_ok=True)
        manifest["status"] = "failed"
        manifest["failure_category"] = exc.category
        manifest["completed_at_utc"] = utc_now()
        atomic_json(manifest_path, manifest)
        print(f"Run {run_id}: FAILED ({exc.category})", file=sys.stderr)
        print(str(exc), file=sys.stderr)
        return 1, run_dir
    except (OSError, subprocess.SubprocessError) as exc:
        (run_dir / "SUCCESS").unlink(missing_ok=True)
        manifest["status"] = "failed"
        manifest["failure_category"] = "runner_error"
        manifest["completed_at_utc"] = utc_now()
        atomic_json(manifest_path, manifest)
        print(f"Run {run_id}: FAILED (runner_error)", file=sys.stderr)
        print(type(exc).__name__, file=sys.stderr)
        return 1, run_dir


def main(argv: Iterable[str] | None = None) -> int:
    args = parser().parse_args(argv)
    try:
        return run(args)[0]
    except (OSError, subprocess.SubprocessError, RunFailure) as exc:
        category = exc.category if isinstance(exc, RunFailure) else type(exc).__name__
        print(f"Runner failed before launch: {category}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
