#!/usr/bin/env python3
"""Validate the tracked, public surface of this repository."""

from __future__ import annotations

import csv
import json
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path, PurePosixPath
from typing import Iterable, Mapping, Sequence

ROOT = Path(__file__).resolve().parents[1]

EXPECTED_IDENTITY = {
    "doi": "10.1016/j.chest.2025.08.002",
    "pmid": "40885535",
    "pmcid": "PMC12739763",
    "repository": "https://github.com/reblocke/Hypercap-Case-Definitions",
}

ALLOWED_CSV_PATHS = {
    PurePosixPath("data_dictionary.csv"),
    PurePosixPath("metadata/output_manifest.csv"),
    PurePosixPath("metadata/phenotype_definitions.csv"),
    PurePosixPath("metadata/stata_dependencies.csv"),
}

BANNED_SUFFIXES = {
    ".dta",
    ".gph",
    ".log",
    ".pdf",
    ".png",
    ".sas7bdat",
    ".smcl",
    ".ster",
    ".svg",
    ".tif",
    ".tiff",
    ".xls",
    ".xlsm",
    ".xlsx",
}

BANNED_TOP_LEVEL_DIRECTORIES = {"Data", "Results and Figures", "data", "outputs"}
CONFLICT_MARKER = re.compile(r"(?m)^(?:<{7}(?: .*)?|={7}|>{7}(?: .*)?)$")
VALID_REVIEW_STATUSES = {"blocked", "draft", "needs_review", "verified"}
VALID_WORKFLOW_ROLES = {"context_only", "derived_analysis", "runtime_input"}
VALID_UPSTREAM_STATUSES = {"blocked", "needs_review", "not_applicable", "verified"}


def tracked_files(root: Path = ROOT) -> list[PurePosixPath]:
    """Return paths known to Git, excluding submodules."""

    result = subprocess.run(
        ["git", "ls-files", "-z"],
        cwd=root,
        check=True,
        capture_output=True,
    )
    return [
        PurePosixPath(item.decode("utf-8"))
        for item in result.stdout.split(b"\0")
        if item
    ]


def path_issues(path: PurePosixPath) -> list[str]:
    """Return public-surface violations implied by a tracked path."""

    issues: list[str] = []
    if path.parts and path.parts[0] in BANNED_TOP_LEVEL_DIRECTORIES:
        issues.append(f"{path}: tracked restricted or generated directory")
    if path.suffix.lower() == ".csv" and path not in ALLOWED_CSV_PATHS:
        issues.append(f"{path}: CSV is not on the explicit public allowlist")
    if path.suffix.lower() in BANNED_SUFFIXES:
        issues.append(f"{path}: tracked restricted or generated file type")
    return issues


def content_issues(path: PurePosixPath, text: str) -> list[str]:
    """Return violations found in tracked text content."""

    issues: list[str] = []
    if CONFLICT_MARKER.search(text):
        issues.append(f"{path}: unresolved merge-conflict marker")

    unix_user_prefix = "/" + "Users" + "/"
    unix_home_prefix = "/" + "home" + "/"
    windows_user_prefix = "C:" + chr(92) + "Users" + chr(92)
    if unix_user_prefix in text:
        issues.append(f"{path}: macOS user-home path")
    if re.search(re.escape(unix_home_prefix) + r"[^/\s]+/", text):
        issues.append(f"{path}: Linux user-home path")
    if windows_user_prefix.lower() in text.lower():
        issues.append(f"{path}: Windows user-home path")
    return issues


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def validate_notebook(path: Path, label: PurePosixPath) -> list[str]:
    issues: list[str] = []
    try:
        notebook = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        return [f"{label}: invalid notebook JSON: {exc}"]

    for index, cell in enumerate(notebook.get("cells", [])):
        if cell.get("cell_type") != "code":
            continue
        if cell.get("execution_count") is not None:
            issues.append(f"{label}: code cell {index} has an execution count")
        if cell.get("outputs"):
            issues.append(f"{label}: code cell {index} has embedded outputs")
    return issues


def validate_phenotype_rows(
    rows: Sequence[Mapping[str, str]],
    root: Path = ROOT,
) -> list[str]:
    issues: list[str] = []
    expected = {f"def{index}" for index in range(1, 11)}
    identifiers = [row.get("definition_id", "") for row in rows]
    counts = Counter(identifiers)

    if set(identifiers) != expected or len(identifiers) != 10:
        issues.append("metadata/phenotype_definitions.csv: expected exactly def1-def10")
    duplicates = sorted(key for key, count in counts.items() if key and count > 1)
    if duplicates:
        issues.append(
            "metadata/phenotype_definitions.csv: duplicate IDs: "
            + ", ".join(duplicates)
        )

    required_fields = {
        "definition_id",
        "implemented_rule",
        "missingness_behavior",
        "published_source",
        "published_rule_summary",
        "code_location",
        "implementation_status",
        "source_verification_status",
        "approval_status",
    }
    for index, row in enumerate(rows, start=2):
        missing = sorted(field for field in required_fields if not row.get(field, "").strip())
        if missing:
            issues.append(
                f"metadata/phenotype_definitions.csv:{index}: blank required fields: "
                + ", ".join(missing)
            )
        if row.get("implementation_status") != "observed_from_code":
            issues.append(
                f"metadata/phenotype_definitions.csv:{index}: "
                "implementation status must be observed_from_code"
            )
        if row.get("source_verification_status") not in {
            "blocked",
            "needs_review",
        }:
            issues.append(
                f"metadata/phenotype_definitions.csv:{index}: "
                "source verification must remain blocked or needs_review"
            )
        if row.get("approval_status") != "unapproved":
            issues.append(
                f"metadata/phenotype_definitions.csv:{index}: "
                "approval status must remain unapproved"
            )

        location = row.get("code_location", "").strip()
        location_match = re.fullmatch(r"(.+):([1-9]\d*)", location)
        if location and not location_match:
            issues.append(
                f"metadata/phenotype_definitions.csv:{index}: "
                "code_location must be a repository-relative file:positive-line"
            )
            continue
        if not location_match:
            continue

        relative_text, line_text = location_match.groups()
        relative = PurePosixPath(relative_text)
        if relative.is_absolute() or ".." in relative.parts:
            issues.append(
                f"metadata/phenotype_definitions.csv:{index}: "
                "code_location must remain within the repository"
            )
            continue

        source_path = root.joinpath(*relative.parts)
        if not source_path.is_file():
            issues.append(
                f"metadata/phenotype_definitions.csv:{index}: "
                f"code_location file does not exist: {relative}"
            )
            continue
        try:
            source_lines = source_path.read_text(encoding="utf-8").splitlines()
        except (OSError, UnicodeDecodeError) as exc:
            issues.append(
                f"metadata/phenotype_definitions.csv:{index}: "
                f"cannot read code_location file: {exc}"
            )
            continue

        line_number = int(line_text)
        if line_number > len(source_lines):
            issues.append(
                f"metadata/phenotype_definitions.csv:{index}: "
                f"code_location line {line_number} is out of range"
            )
            continue

        definition_id = row.get("definition_id", "").strip()
        definition_pattern = re.compile(
            rf"^\s*(?:gen|generate)\s+{re.escape(definition_id)}\s*=",
            re.IGNORECASE,
        )
        if not definition_pattern.search(source_lines[line_number - 1]):
            issues.append(
                f"metadata/phenotype_definitions.csv:{index}: "
                f"code_location does not define {definition_id}"
            )
    return issues


def validate_dictionary_rows(rows: Sequence[Mapping[str, str]]) -> list[str]:
    issues: list[str] = []
    roles = Counter(row.get("workflow_role", "") for row in rows)
    required_fields = {
        "allowed_values",
        "definition",
        "readable_name",
        "review_status",
        "type",
        "upstream_derivation_status",
        "variable_name",
        "workflow_role",
    }

    if len(rows) != 90:
        issues.append(f"data_dictionary.csv: expected 90 variable rows, found {len(rows)}")
    expected_roles = {
        "runtime_input": 56,
        "context_only": 7,
        "derived_analysis": 27,
    }
    if roles != Counter(expected_roles):
        issues.append(
            "data_dictionary.csv: workflow-role counts must be "
            "56 runtime_input, 7 context_only, and 27 derived_analysis"
        )

    names = [row.get("variable_name", "") for row in rows]
    duplicates = sorted(name for name, count in Counter(names).items() if name and count > 1)
    if duplicates:
        issues.append("data_dictionary.csv: duplicate variables: " + ", ".join(duplicates))
    if "full_db.dta" in names:
        issues.append("data_dictionary.csv: input artifact must not be a variable row")

    for index, row in enumerate(rows, start=2):
        missing = sorted(field for field in required_fields if not row.get(field, "").strip())
        if missing:
            issues.append(
                f"data_dictionary.csv:{index}: blank required fields: "
                + ", ".join(missing)
            )
        if row.get("source_or_origin") == "generated output":
            issues.append(
                f"data_dictionary.csv:{index}: output artifact must be in output manifest"
            )
        if row.get("workflow_role") not in VALID_WORKFLOW_ROLES:
            issues.append(f"data_dictionary.csv:{index}: invalid workflow_role")
        if row.get("upstream_derivation_status") not in VALID_UPSTREAM_STATUSES:
            issues.append(
                f"data_dictionary.csv:{index}: invalid upstream_derivation_status"
            )
        if row.get("review_status") not in VALID_REVIEW_STATUSES:
            issues.append(f"data_dictionary.csv:{index}: invalid review_status")
        if row.get("review_status") == "reviewed_from_code":
            issues.append(
                f"data_dictionary.csv:{index}: reviewed_from_code is not a verification status"
            )
    return issues


def _parse_flat_yaml(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        if ":" not in line:
            raise ValueError(f"expected key: value line, found {raw_line!r}")
        key, value = line.split(":", 1)
        values[key.strip()] = value.strip()
    return values


def validate_metadata(root: Path = ROOT) -> list[str]:
    issues: list[str] = []

    phenotype_path = root / "metadata/phenotype_definitions.csv"
    dictionary_path = root / "data_dictionary.csv"
    outputs_path = root / "metadata/output_manifest.csv"
    dependencies_path = root / "metadata/stata_dependencies.csv"
    upstream_path = root / "metadata/upstream_dependency.yml"

    required_paths = [
        phenotype_path,
        dictionary_path,
        outputs_path,
        dependencies_path,
        upstream_path,
    ]
    for path in required_paths:
        if not path.is_file():
            issues.append(f"{path.relative_to(root)}: required metadata file is missing")
    if issues:
        return issues

    issues.extend(validate_phenotype_rows(_read_csv(phenotype_path), root))
    issues.extend(validate_dictionary_rows(_read_csv(dictionary_path)))

    output_rows = _read_csv(outputs_path)
    output_ids = [row.get("output_id", "") for row in output_rows]
    if len(output_rows) != 12 or len(set(output_ids)) != 12:
        issues.append("metadata/output_manifest.csv: expected 12 unique outputs")
    expected_counts = {
        "run_provenance": 9,
        "overall_cohort": 1,
        "definition_overlap_heatmaps": 4,
        "testing_strategy_kappa": 3,
        "definition_summary": 10,
        "workup_summary": 1,
        "location_summary": 1,
        "location_kappa": 4,
        "encounter_spline_intermediates": 4,
        "paco2_spline": 1,
        "location_splines": 9,
        "consort_diagram": 1,
    }
    for index, row in enumerate(output_rows, start=2):
        if row.get("tracked") != "false":
            issues.append(
                f"metadata/output_manifest.csv:{index}: generated output must be untracked"
            )
        if not row.get("path_pattern", "").startswith("outputs/"):
            issues.append(
                f"metadata/output_manifest.csv:{index}: output must be under outputs/"
            )
        output_id = row.get("output_id", "")
        expected_required = "false" if output_id == "consort_diagram" else "true"
        if row.get("required_for_stata_success") != expected_required:
            issues.append(
                f"metadata/output_manifest.csv:{index}: invalid Stata success requirement"
            )
        try:
            observed_count = int(row.get("expected_count", ""))
        except ValueError:
            observed_count = -1
        if observed_count != expected_counts.get(output_id):
            issues.append(
                f"metadata/output_manifest.csv:{index}: unexpected artifact count"
            )
        if not row.get("validation_rule", "").strip():
            issues.append(
                f"metadata/output_manifest.csv:{index}: validation rule is required"
            )

    dependency_rows = _read_csv(dependencies_path)
    expected_dependencies = {
        "cleanplots",
        "colorpalette",
        "colrspace",
        "diagt",
        "gtools",
        "heatplot",
        "kappaetc",
        "missings",
        "mkspline2",
        "moremata",
        "table1_mc",
        "xblc",
    }
    dependency_names = {row.get("name", "") for row in dependency_rows}
    if dependency_names != expected_dependencies or len(dependency_rows) != 12:
        issues.append(
            "metadata/stata_dependencies.csv: expected 12 direct and transitive dependencies"
        )
    for index, row in enumerate(dependency_rows, start=2):
        if row.get("publication_required_version") != "UNRESOLVED":
            issues.append(
                f"metadata/stata_dependencies.csv:{index}: publication version "
                "must remain UNRESOLVED"
            )
        expected_required = "false" if row.get("name") == "gtools" else "true"
        if row.get("required_for_current_analysis") != expected_required:
            issues.append(
                f"metadata/stata_dependencies.csv:{index}: invalid required flag"
            )
        if not row.get("preflight_check", "").strip():
            issues.append(
                f"metadata/stata_dependencies.csv:{index}: preflight check is required"
            )

    try:
        upstream = _parse_flat_yaml(upstream_path)
    except (OSError, ValueError) as exc:
        issues.append(f"metadata/upstream_dependency.yml: {exc}")
    else:
        expected_upstream = {
            "producer_commit": "UNRESOLVED",
            "input_schema_version": "UNRESOLVED",
            "expected_artifact": "full_db.dta",
            "observed_validation_checkout_commit": (
                "1185a6bc9957a02cb24be5f1f7fa10c48d8a4c13"
            ),
            "observed_validation_artifact_path": (
                "Data/derived/hypercapnia/preprocessing/full_db.dta"
            ),
            "access_classification": "restricted",
            "redistributable": "false",
            "verification_status": "blocked",
        }
        for key, expected in expected_upstream.items():
            if upstream.get(key) != expected:
                issues.append(
                    f"metadata/upstream_dependency.yml: {key} must be {expected}"
                )

    return issues


def validate_identity(root: Path = ROOT) -> list[str]:
    issues: list[str] = []
    requirements = {
        "README.md": EXPECTED_IDENTITY.values(),
        "llms.txt": EXPECTED_IDENTITY.values(),
        "CITATION.cff": (
            EXPECTED_IDENTITY["doi"],
            EXPECTED_IDENTITY["repository"],
        ),
    }
    for relative, tokens in requirements.items():
        path = root / relative
        if not path.is_file():
            issues.append(f"{relative}: required identity file is missing")
            continue
        text = path.read_text(encoding="utf-8")
        for token in tokens:
            if token not in text:
                issues.append(f"{relative}: missing canonical identity token {token}")

    command_files = ("README.md", "llms.txt", "AGENTS.md")
    for relative in command_files:
        path = root / relative
        if not path.is_file():
            issues.append(f"{relative}: required workflow file is missing")
            continue
        text = path.read_text(encoding="utf-8")
        for command in (
            "make check",
            "make diagram-smoke",
            "make stata-run",
            "make stata-compare",
        ):
            if command not in text:
                issues.append(f"{relative}: missing canonical command {command}")
    return issues


def find_issues(
    root: Path = ROOT,
    paths: Iterable[PurePosixPath] | None = None,
) -> list[str]:
    issues: list[str] = []
    selected = list(paths) if paths is not None else tracked_files(root)

    for relative in selected:
        issues.extend(path_issues(relative))
        path = root / relative
        if not path.is_file():
            continue
        if relative.suffix.lower() == ".ipynb":
            issues.extend(validate_notebook(path, relative))
        try:
            text = path.read_text(encoding="utf-8")
        except (UnicodeDecodeError, OSError):
            continue
        issues.extend(content_issues(relative, text))

    issues.extend(validate_identity(root))
    issues.extend(validate_metadata(root))
    return sorted(set(issues))


def main() -> int:
    issues = find_issues()
    if issues:
        print("Public-surface audit failed:", file=sys.stderr)
        for issue in issues:
            print(f"- {issue}", file=sys.stderr)
        return 1
    print("Public-surface audit passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
