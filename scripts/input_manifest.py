#!/usr/bin/env python3
"""Approve and validate the restricted Stata input without exposing its contents."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence

INPUT_FILENAME = "full_db.dta"
MANIFEST_FILENAME = "full_db.manifest.json"
MANIFEST_SCHEMA_VERSION = 1
APPROVED_VERIFICATION_BASIS = "scientific_owner_approved_historical_evidence"
APPROVED_VERIFICATION_SCOPE = "artifact_producer_and_observed_schema_only"
MANIFEST_KEYS = frozenset(
    {
        "schema_version",
        "logical_name",
        "size_bytes",
        "sha256",
        "approved_at_utc",
        "upstream_repository",
        "producer_commit",
        "input_schema_version",
        "data_dictionary_sha256",
    }
)
SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")
COMMIT_PATTERN = re.compile(r"^[0-9a-f]{40}$")
SCHEMA_NAME_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
UTC_TIMESTAMP_PATTERN = re.compile(r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z$")


class InputManifestFailure(RuntimeError):
    """Expected approval or validation failure with a stable category."""

    def __init__(self, category: str, message: str):
        super().__init__(message)
        self.category = category


def _utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _hash_and_size(path: Path) -> tuple[str, int]:
    digest = hashlib.sha256()
    size = 0
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
            size += len(chunk)
    return digest.hexdigest(), size


def _sha256_file(path: Path) -> str:
    return _hash_and_size(path)[0]


def _parse_flat_yaml(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except (OSError, UnicodeError) as error:
        raise InputManifestFailure(
            "input_manifest_unapproved",
            "The repository input-provenance authority is unavailable.",
        ) from error

    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        if ":" not in line:
            raise InputManifestFailure(
                "input_manifest_unapproved",
                "The repository input-provenance authority is invalid.",
            )
        key, raw_value = line.split(":", 1)
        key = key.strip()
        value = raw_value.strip()
        if not key or key in values:
            raise InputManifestFailure(
                "input_manifest_unapproved",
                "The repository input-provenance authority is invalid.",
            )
        if (
            len(value) >= 2
            and value[0] == value[-1]
            and value[0] in {'"', "'"}
        ):
            value = value[1:-1]
        values[key] = value
    return values


def _authority(analysis_root: Path) -> dict[str, str]:
    metadata_path = analysis_root / "metadata" / "upstream_dependency.yml"
    dictionary_path = analysis_root / "data_dictionary.csv"
    if not metadata_path.is_file() or not dictionary_path.is_file():
        raise InputManifestFailure(
            "input_manifest_unapproved",
            "The repository input-provenance authority is incomplete.",
        )

    metadata = _parse_flat_yaml(metadata_path)
    repository = metadata.get("upstream_repository", "")
    producer_commit = metadata.get("producer_commit", "")
    input_schema_version = metadata.get("input_schema_version", "")
    expected_artifact = metadata.get("expected_artifact", "")
    verification_status = metadata.get("verification_status", "")
    verification_basis = metadata.get("verification_basis", "")
    verification_scope = metadata.get("verification_scope", "")

    if (
        verification_status != "verified"
        or verification_basis != APPROVED_VERIFICATION_BASIS
        or verification_scope != APPROVED_VERIFICATION_SCOPE
        or expected_artifact != INPUT_FILENAME
        or not repository
        or repository == "UNRESOLVED"
        or not COMMIT_PATTERN.fullmatch(producer_commit)
        or not SCHEMA_NAME_PATTERN.fullmatch(input_schema_version)
        or input_schema_version == "UNRESOLVED"
    ):
        raise InputManifestFailure(
            "input_manifest_unapproved",
            "The repository has not approved the restricted input provenance.",
        )

    try:
        dictionary_hash = _sha256_file(dictionary_path)
    except OSError as error:
        raise InputManifestFailure(
            "input_manifest_unapproved",
            "The repository data contract is unavailable.",
        ) from error

    return {
        "upstream_repository": repository,
        "producer_commit": producer_commit,
        "input_schema_version": input_schema_version,
        "data_dictionary_sha256": dictionary_hash,
    }


def _atomic_json(path: Path, payload: dict[str, Any]) -> None:
    temporary_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary_path = Path(handle.name)
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_path, path)
    except OSError:
        if temporary_path is not None:
            try:
                temporary_path.unlink(missing_ok=True)
            except OSError:
                pass
        raise


def approve_input(
    input_root: Path,
    analysis_root: Path,
    approval: str,
    replace: str = "",
) -> Path:
    """Create an adjacent approval manifest after explicit confirmation."""

    if approval != "YES":
        raise InputManifestFailure(
            "input_manifest_confirmation",
            "Input approval requires the literal confirmation YES.",
        )
    if replace not in {"", "YES"}:
        raise InputManifestFailure(
            "input_manifest_confirmation",
            "Manifest replacement requires the literal confirmation YES.",
        )

    input_file = input_root / INPUT_FILENAME
    manifest_path = input_root / MANIFEST_FILENAME
    if not input_file.is_file():
        raise InputManifestFailure(
            "missing_input",
            "The required restricted input file is unavailable.",
        )
    if manifest_path.exists() and replace != "YES":
        raise InputManifestFailure(
            "input_manifest_exists",
            "An input approval manifest already exists; replacement was not confirmed.",
        )

    authority = _authority(analysis_root)
    try:
        input_hash, input_size = _hash_and_size(input_file)
    except OSError as error:
        raise InputManifestFailure(
            "missing_input",
            "The required restricted input file could not be read.",
        ) from error

    payload: dict[str, Any] = {
        "schema_version": MANIFEST_SCHEMA_VERSION,
        "logical_name": INPUT_FILENAME,
        "size_bytes": input_size,
        "sha256": input_hash,
        "approved_at_utc": _utc_now(),
        **authority,
    }
    try:
        _atomic_json(manifest_path, payload)
    except OSError as error:
        raise InputManifestFailure(
            "input_manifest_unapproved",
            "The input approval manifest could not be written.",
        ) from error
    return manifest_path


def _read_manifest(path: Path) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise InputManifestFailure(
            "input_manifest_malformed",
            "The input approval manifest is malformed.",
        ) from error
    if not isinstance(payload, dict) or set(payload) != MANIFEST_KEYS:
        raise InputManifestFailure(
            "input_manifest_malformed",
            "The input approval manifest does not have the required schema.",
        )

    if (
        type(payload["schema_version"]) is not int
        or payload["schema_version"] != MANIFEST_SCHEMA_VERSION
        or type(payload["logical_name"]) is not str
        or payload["logical_name"] != INPUT_FILENAME
        or type(payload["size_bytes"]) is not int
        or payload["size_bytes"] < 0
        or type(payload["sha256"]) is not str
        or not SHA256_PATTERN.fullmatch(payload["sha256"])
        or type(payload["approved_at_utc"]) is not str
        or not UTC_TIMESTAMP_PATTERN.fullmatch(payload["approved_at_utc"])
        or type(payload["upstream_repository"]) is not str
        or not payload["upstream_repository"]
        or type(payload["producer_commit"]) is not str
        or not COMMIT_PATTERN.fullmatch(payload["producer_commit"])
        or type(payload["input_schema_version"]) is not str
        or not SCHEMA_NAME_PATTERN.fullmatch(payload["input_schema_version"])
        or type(payload["data_dictionary_sha256"]) is not str
        or not SHA256_PATTERN.fullmatch(payload["data_dictionary_sha256"])
    ):
        raise InputManifestFailure(
            "input_manifest_malformed",
            "The input approval manifest contains invalid field values.",
        )

    try:
        datetime.strptime(payload["approved_at_utc"], "%Y-%m-%dT%H:%M:%SZ")
    except ValueError as error:
        raise InputManifestFailure(
            "input_manifest_malformed",
            "The input approval manifest contains an invalid approval time.",
        ) from error
    return payload


def validate_approved_input(input_root: Path, analysis_root: Path) -> dict[str, Any]:
    """Validate the adjacent approval, repository provenance, and input bytes."""

    input_file = input_root / INPUT_FILENAME
    manifest_path = input_root / MANIFEST_FILENAME
    if not input_file.is_file():
        raise InputManifestFailure(
            "missing_input",
            "The required restricted input file is unavailable.",
        )
    if not manifest_path.is_file():
        raise InputManifestFailure(
            "input_manifest_missing",
            "The required input approval manifest is unavailable.",
        )

    manifest = _read_manifest(manifest_path)
    authority = _authority(analysis_root)

    provenance_fields = (
        "upstream_repository",
        "producer_commit",
        "input_schema_version",
    )
    if any(manifest[field] != authority[field] for field in provenance_fields):
        raise InputManifestFailure(
            "input_manifest_provenance_mismatch",
            "The approved input provenance does not match the repository authority.",
        )
    if manifest["data_dictionary_sha256"] != authority["data_dictionary_sha256"]:
        raise InputManifestFailure(
            "input_manifest_contract_mismatch",
            "The approved input data contract does not match the repository contract.",
        )

    try:
        input_hash, input_size = _hash_and_size(input_file)
    except OSError as error:
        raise InputManifestFailure(
            "missing_input",
            "The required restricted input file could not be read.",
        ) from error
    if manifest["sha256"] != input_hash or manifest["size_bytes"] != input_size:
        raise InputManifestFailure(
            "input_manifest_input_mismatch",
            "The restricted input does not match its approval manifest.",
        )

    try:
        manifest_hash = _sha256_file(manifest_path)
    except OSError as error:
        raise InputManifestFailure(
            "input_manifest_malformed",
            "The input approval manifest could not be read.",
        ) from error

    return {
        "input_file": input_file,
        "input": {
            "logical_name": INPUT_FILENAME,
            "size_bytes": input_size,
            "sha256": input_hash,
        },
        "approval": {
            "manifest_sha256": manifest_hash,
            **authority,
        },
    }


def parser() -> argparse.ArgumentParser:
    command_parser = argparse.ArgumentParser(
        description="Approve a restricted full_db.dta without printing sensitive details."
    )
    subparsers = command_parser.add_subparsers(dest="command", required=True)
    approve_parser = subparsers.add_parser(
        "approve",
        help="Create or explicitly replace the adjacent approval manifest.",
    )
    approve_parser.add_argument("--input-root", type=Path, required=True)
    approve_parser.add_argument("--analysis-root", type=Path, required=True)
    approve_parser.add_argument("--approve", required=True)
    approve_parser.add_argument("--replace", default="")
    return command_parser


def main(argv: Sequence[str] | None = None) -> int:
    args = parser().parse_args(argv)
    try:
        if args.command == "approve":
            approve_input(
                input_root=args.input_root,
                analysis_root=args.analysis_root,
                approval=args.approve,
                replace=args.replace,
            )
    except InputManifestFailure as error:
        print(f"ERROR [{error.category}]: {error}", file=sys.stderr)
        return 1
    print("Approved restricted input manifest.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
