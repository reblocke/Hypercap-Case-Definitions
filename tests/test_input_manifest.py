from __future__ import annotations

import contextlib
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import input_manifest  # noqa: E402


PRODUCER_COMMIT = "44f49748d415e92b7d50b50d86b8fdea29f6cb07"
UPSTREAM_REPOSITORY = "https://github.com/reblocke/trinetx-hypercapnia-code"
INPUT_SCHEMA_VERSION = "hypercapnia-full-db-v1"


def write_authority(root: Path) -> None:
    metadata = root / "metadata"
    metadata.mkdir()
    (metadata / "upstream_dependency.yml").write_text(
        "\n".join(
            (
                "schema_version: 1",
                f"upstream_repository: {UPSTREAM_REPOSITORY}",
                f"producer_commit: {PRODUCER_COMMIT}",
                f"input_schema_version: {INPUT_SCHEMA_VERSION}",
                "expected_artifact: full_db.dta",
                "verification_status: verified",
                "verification_basis: scientific_owner_approved_historical_evidence",
                "verification_scope: artifact_producer_and_observed_schema_only",
                "",
            )
        ),
        encoding="utf-8",
    )
    (root / "data_dictionary.csv").write_text(
        "variable_name,description\nexample,public contract\n",
        encoding="utf-8",
    )


class InputManifestTests(unittest.TestCase):
    def make_fixture(self, base: Path) -> tuple[Path, Path]:
        analysis_root = base / "analysis"
        input_root = base / "restricted-input"
        analysis_root.mkdir()
        input_root.mkdir()
        write_authority(analysis_root)
        (input_root / "full_db.dta").write_bytes(b"restricted fixture bytes")
        return input_root, analysis_root

    def assert_category(
        self,
        expected: str,
        callable_object: object,
        *args: object,
        **kwargs: object,
    ) -> input_manifest.InputManifestFailure:
        with self.assertRaises(input_manifest.InputManifestFailure) as raised:
            callable_object(*args, **kwargs)  # type: ignore[operator]
        self.assertEqual(expected, raised.exception.category)
        return raised.exception

    def test_approval_requires_literal_confirmation(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            input_root, analysis_root = self.make_fixture(Path(tmp))
            error = self.assert_category(
                "input_manifest_confirmation",
                input_manifest.approve_input,
                input_root,
                analysis_root,
                "yes",
            )
        self.assertNotIn(str(input_root), str(error))

    def test_existing_manifest_requires_explicit_replacement(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            input_root, analysis_root = self.make_fixture(Path(tmp))
            manifest_path = input_manifest.approve_input(
                input_root, analysis_root, "YES"
            )
            original = manifest_path.read_bytes()
            self.assert_category(
                "input_manifest_exists",
                input_manifest.approve_input,
                input_root,
                analysis_root,
                "YES",
            )
            self.assertEqual(original, manifest_path.read_bytes())
            input_manifest.approve_input(
                input_root, analysis_root, "YES", replace="YES"
            )
        self.assertTrue(original)

    def test_approval_is_atomic_path_free_and_exact_schema(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            input_root, analysis_root = self.make_fixture(base)
            manifest_path = input_manifest.approve_input(
                input_root, analysis_root, "YES"
            )
            payload = json.loads(manifest_path.read_text(encoding="utf-8"))
            temporary_files = list(input_root.glob("*.tmp"))
        self.assertEqual(input_manifest.MANIFEST_KEYS, set(payload))
        self.assertEqual(1, payload["schema_version"])
        self.assertEqual("full_db.dta", payload["logical_name"])
        self.assertEqual(PRODUCER_COMMIT, payload["producer_commit"])
        self.assertEqual(INPUT_SCHEMA_VERSION, payload["input_schema_version"])
        self.assertNotIn(str(base), json.dumps(payload))
        self.assertNotIn("row", json.dumps(payload).lower())
        self.assertNotIn("value", json.dumps(payload).lower())
        self.assertEqual([], temporary_files)

    def test_validate_returns_only_bounded_approval_details(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            input_root, analysis_root = self.make_fixture(Path(tmp))
            input_manifest.approve_input(input_root, analysis_root, "YES")
            result = input_manifest.validate_approved_input(
                input_root, analysis_root
            )
        self.assertEqual(input_root / "full_db.dta", result["input_file"])
        self.assertEqual(
            {"logical_name", "size_bytes", "sha256"}, set(result["input"])
        )
        self.assertEqual(
            {
                "manifest_sha256",
                "upstream_repository",
                "producer_commit",
                "input_schema_version",
                "data_dictionary_sha256",
            },
            set(result["approval"]),
        )

    def test_missing_input_and_manifest_have_stable_categories(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            input_root, analysis_root = self.make_fixture(base)
            (input_root / "full_db.dta").unlink()
            self.assert_category(
                "missing_input",
                input_manifest.validate_approved_input,
                input_root,
                analysis_root,
            )
            (input_root / "full_db.dta").write_bytes(b"restored")
            self.assert_category(
                "input_manifest_missing",
                input_manifest.validate_approved_input,
                input_root,
                analysis_root,
            )

    def test_malformed_json_extra_keys_and_invalid_fields_fail(self) -> None:
        invalid_payloads: list[object] = [
            "{not json",
            {"unexpected": True},
        ]
        for invalid in invalid_payloads:
            with self.subTest(invalid=invalid):
                with tempfile.TemporaryDirectory() as tmp:
                    input_root, analysis_root = self.make_fixture(Path(tmp))
                    manifest_path = input_root / input_manifest.MANIFEST_FILENAME
                    if isinstance(invalid, str):
                        manifest_path.write_text(invalid, encoding="utf-8")
                    else:
                        manifest_path.write_text(
                            json.dumps(invalid), encoding="utf-8"
                        )
                    self.assert_category(
                        "input_manifest_malformed",
                        input_manifest.validate_approved_input,
                        input_root,
                        analysis_root,
                    )

        invalid_fields = {
            "schema_version": 2,
            "logical_name": "other.dta",
            "approved_at_utc": "2026-99-99T99:99:99Z",
            "sha256": "A" * 64,
            "size_bytes": True,
        }
        for field, value in invalid_fields.items():
            with self.subTest(field=field):
                with tempfile.TemporaryDirectory() as tmp:
                    input_root, analysis_root = self.make_fixture(Path(tmp))
                    manifest_path = input_manifest.approve_input(
                        input_root, analysis_root, "YES"
                    )
                    payload = json.loads(manifest_path.read_text())
                    payload[field] = value
                    manifest_path.write_text(json.dumps(payload), encoding="utf-8")
                    self.assert_category(
                        "input_manifest_malformed",
                        input_manifest.validate_approved_input,
                        input_root,
                        analysis_root,
                    )

    def test_unapproved_repository_authority_fails_safely(self) -> None:
        mutations = (
            ("verification_status: verified", "verification_status: blocked"),
            (
                "verification_basis: scientific_owner_approved_historical_evidence",
                "verification_basis: unsupported",
            ),
            (
                "verification_scope: artifact_producer_and_observed_schema_only",
                "verification_scope: overly_broad",
            ),
        )
        for original, replacement in mutations:
            with self.subTest(replacement=replacement):
                with tempfile.TemporaryDirectory() as tmp:
                    input_root, analysis_root = self.make_fixture(Path(tmp))
                    metadata_path = (
                        analysis_root / "metadata" / "upstream_dependency.yml"
                    )
                    metadata_path.write_text(
                        metadata_path.read_text().replace(original, replacement),
                        encoding="utf-8",
                    )
                    error = self.assert_category(
                        "input_manifest_unapproved",
                        input_manifest.approve_input,
                        input_root,
                        analysis_root,
                        "YES",
                    )
                self.assertNotIn(PRODUCER_COMMIT, str(error))
                self.assertNotIn(str(input_root), str(error))

    def test_input_drift_fails(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            input_root, analysis_root = self.make_fixture(Path(tmp))
            input_manifest.approve_input(input_root, analysis_root, "YES")
            (input_root / "full_db.dta").write_bytes(b"changed restricted bytes")
            self.assert_category(
                "input_manifest_input_mismatch",
                input_manifest.validate_approved_input,
                input_root,
                analysis_root,
            )

    def test_provenance_drift_fails(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            input_root, analysis_root = self.make_fixture(Path(tmp))
            manifest_path = input_manifest.approve_input(
                input_root, analysis_root, "YES"
            )
            metadata_path = analysis_root / "metadata" / "upstream_dependency.yml"
            metadata_path.write_text(
                metadata_path.read_text().replace(
                    INPUT_SCHEMA_VERSION,
                    "hypercapnia-full-db-v2",
                ),
                encoding="utf-8",
            )
            error = self.assert_category(
                "input_manifest_provenance_mismatch",
                input_manifest.validate_approved_input,
                input_root,
                analysis_root,
            )
        self.assertNotIn(str(manifest_path), str(error))

    def test_data_contract_drift_fails(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            input_root, analysis_root = self.make_fixture(Path(tmp))
            input_manifest.approve_input(input_root, analysis_root, "YES")
            (analysis_root / "data_dictionary.csv").write_text(
                "variable_name,description\nchanged,new contract\n",
                encoding="utf-8",
            )
            self.assert_category(
                "input_manifest_contract_mismatch",
                input_manifest.validate_approved_input,
                input_root,
                analysis_root,
            )

    def test_cli_output_never_prints_paths_or_hashes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            input_root, analysis_root = self.make_fixture(Path(tmp))
            stdout = io.StringIO()
            stderr = io.StringIO()
            with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
                code = input_manifest.main(
                    [
                        "approve",
                        "--input-root",
                        str(input_root),
                        "--analysis-root",
                        str(analysis_root),
                        "--approve",
                        "YES",
                    ]
                )
            output = stdout.getvalue() + stderr.getvalue()
            manifest = json.loads(
                (input_root / input_manifest.MANIFEST_FILENAME).read_text()
            )
        self.assertEqual(0, code)
        self.assertNotIn(str(input_root), output)
        self.assertNotIn(manifest["sha256"], output)
        self.assertEqual("Approved restricted input manifest.\n", stdout.getvalue())


if __name__ == "__main__":
    unittest.main()
