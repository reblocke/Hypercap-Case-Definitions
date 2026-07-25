from __future__ import annotations

import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path, PurePosixPath

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import check_public_surface as audit  # noqa: E402


def valid_upstream_record() -> dict[str, str]:
    return {
        "schema_version": "2",
        "upstream_repository": (
            "https://github.com/reblocke/trinetx-hypercapnia-code"
        ),
        "producer_commit": audit.APPROVED_PRODUCER_COMMIT,
        "input_schema_version": audit.APPROVED_INPUT_SCHEMA,
        "expected_artifact": "full_db.dta",
        "observed_validation_checkout_commit": (
            "1185a6bc9957a02cb24be5f1f7fa10c48d8a4c13"
        ),
        "observed_validation_artifact_path": (
            "Data/derived/hypercapnia/preprocessing/full_db.dta"
        ),
        "access_classification": "restricted",
        "redistributable": "false",
        "verification_status": "verified",
        "verification_basis": audit.APPROVED_VERIFICATION_BASIS,
        "verification_scope": audit.APPROVED_VERIFICATION_SCOPE,
        "verification_date": "2026-07-23",
        "verified_derivations": "hypercap_on_abg,hypercap_resp_failure",
        "derivation_source_file": "stata/do/10_preprocessing.do",
        "derivation_source_blob_git_oid": (
            "5842c635bc37a0477eed351295a8d2109b8037e3"
        ),
        "derivation_source_sha256": (
            "be64c9e09c32fa3615c523ebcfb75dbe631ecb8b87c05e7710fe940134e1e451"
        ),
        "derivation_verification_date": "2026-07-24",
        "historical_clean_worktree_recorded": "false",
        "historical_source_hashes_recorded": "false",
        "notes": (
            "The historical build did not record clean-worktree state or "
            "source-file hashes. This public repository does not reproduce "
            "upstream construction."
        ),
    }


def copy_identity_surface(destination: Path) -> None:
    for relative in (
        "README.md",
        "llms.txt",
        "CITATION.cff",
        "AGENTS.md",
        "data_dictionary.md",
        "docs/REPRODUCIBILITY.md",
    ):
        target = destination / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text((ROOT / relative).read_text(encoding="utf-8"), encoding="utf-8")


class PublicSurfaceUnitTests(unittest.TestCase):
    def test_conflict_marker_is_detected(self) -> None:
        marker = "<" * 7 + " ours\ncontent\n" + "=" * 7 + "\n"
        issues = audit.content_issues(PurePosixPath("example.md"), marker)
        self.assertTrue(any("merge-conflict" in issue for issue in issues))

    def test_local_home_path_is_detected(self) -> None:
        local_path = "/" + "Users" + "/person/project/file.txt"
        issues = audit.content_issues(PurePosixPath("example.md"), local_path)
        self.assertTrue(any("user-home" in issue for issue in issues))

    def test_restricted_path_is_detected(self) -> None:
        issues = audit.path_issues(PurePosixPath("data/private/full_db.dta"))
        self.assertTrue(any("restricted" in issue for issue in issues))

    def test_generated_output_path_is_detected(self) -> None:
        issues = audit.path_issues(PurePosixPath("outputs/figures/result.png"))
        self.assertTrue(any("generated" in issue for issue in issues))

    def test_unapproved_csv_is_detected(self) -> None:
        issues = audit.path_issues(PurePosixPath("private_extract.csv"))
        self.assertTrue(any("allowlist" in issue for issue in issues))

    def test_duplicate_phenotype_is_detected(self) -> None:
        rows = []
        for index in range(1, 11):
            rows.append(
                {
                    "definition_id": f"def{index}",
                    "implemented_rule": "rule",
                    "missingness_behavior": "missing becomes zero",
                    "published_source": "source",
                    "published_rule_summary": "summary",
                    "code_location": "file",
                    "implementation_status": "observed_from_code",
                    "source_verification_status": "needs_review",
                    "approval_status": "unapproved",
                }
            )
        rows[-1]["definition_id"] = "def1"
        issues = audit.validate_phenotype_rows(rows)
        self.assertTrue(any("duplicate IDs" in issue for issue in issues))

    def test_phenotype_cannot_be_silently_approved(self) -> None:
        rows = []
        for index in range(1, 11):
            rows.append(
                {
                    "definition_id": f"def{index}",
                    "implemented_rule": "rule",
                    "missingness_behavior": "missing becomes zero",
                    "published_source": "source",
                    "published_rule_summary": "summary",
                    "code_location": "file",
                    "implementation_status": "observed_from_code",
                    "source_verification_status": "needs_review",
                    "approval_status": "unapproved",
                }
            )
        rows[0]["approval_status"] = "approved"
        issues = audit.validate_phenotype_rows(rows)
        self.assertTrue(any("not in the approved allowlist" in issue for issue in issues))

    def test_approved_phenotype_must_remain_verified_and_approved(self) -> None:
        rows = audit._read_csv(ROOT / "metadata/phenotype_definitions.csv")
        target = next(row for row in rows if row["definition_id"] == "def4")
        target["source_verification_status"] = "needs_review"
        target["approval_status"] = "unapproved"
        issues = audit.validate_phenotype_rows(rows)
        self.assertTrue(
            any("approved definition source must be verified" in issue for issue in issues)
        )
        self.assertTrue(
            any("approved definition must carry approved status" in issue for issue in issues)
        )

    def test_selected_upstream_derivations_are_allowlisted(self) -> None:
        rows = audit._read_csv(ROOT / "data_dictionary.csv")
        aggregate = next(
            row for row in rows if row["variable_name"] == "hypercap_on_abg"
        )
        aggregate["upstream_derivation_status"] = "blocked"
        unselected = next(row for row in rows if row["variable_name"] == "paco2")
        unselected["upstream_derivation_status"] = "verified"
        issues = audit.validate_dictionary_rows(rows)
        self.assertTrue(any("must be verified" in issue for issue in issues))
        self.assertTrue(any("must be blocked" in issue for issue in issues))

    def test_malformed_phenotype_code_location_is_detected(self) -> None:
        rows = [
            {
                "definition_id": "def1",
                "code_location": "analysis.do:not-a-line",
            }
        ]
        issues = audit.validate_phenotype_rows(rows)
        self.assertTrue(
            any("file:positive-line" in issue for issue in issues),
            "\n".join(issues),
        )

    def test_out_of_range_phenotype_code_location_is_detected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "analysis.do").write_text("gen def1 = 1\n", encoding="utf-8")
            rows = [
                {
                    "definition_id": "def1",
                    "code_location": "analysis.do:2",
                }
            ]
            issues = audit.validate_phenotype_rows(rows, root)
        self.assertTrue(
            any("line 2 is out of range" in issue for issue in issues),
            "\n".join(issues),
        )

    def test_wrong_phenotype_code_location_is_detected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "analysis.do").write_text("gen def2 = 1\n", encoding="utf-8")
            rows = [
                {
                    "definition_id": "def1",
                    "code_location": "analysis.do:1",
                }
            ]
            issues = audit.validate_phenotype_rows(rows, root)
        self.assertTrue(
            any("does not define def1" in issue for issue in issues),
            "\n".join(issues),
        )

    def test_identity_drift_is_detected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for name in ("README.md", "llms.txt", "CITATION.cff"):
                (root / name).write_text("incorrect metadata\n", encoding="utf-8")
            issues = audit.validate_identity(root)
        self.assertTrue(any("canonical identity token" in issue for issue in issues))

    def test_upstream_producer_schema_and_verification_are_enforced(self) -> None:
        mutations = {
            "producer_commit": "0" * 40,
            "input_schema_version": "different-schema",
            "verification_status": "blocked",
            "verification_basis": "unsupported_basis",
            "verification_scope": "overly_broad",
        }
        for field, replacement in mutations.items():
            with self.subTest(field=field):
                record = valid_upstream_record()
                record[field] = replacement
                issues = audit.validate_upstream_record(record)
                self.assertTrue(
                    any(field in issue for issue in issues),
                    "\n".join(issues),
                )

    def test_upstream_historical_limitations_are_enforced(self) -> None:
        for omitted in (
            "did not record clean-worktree state",
            "source-file hashes",
            "does not reproduce upstream construction",
        ):
            with self.subTest(omitted=omitted):
                record = valid_upstream_record()
                record["notes"] = record["notes"].replace(omitted, "omitted")
                issues = audit.validate_upstream_record(record)
                self.assertTrue(
                    any(omitted in issue for issue in issues),
                    "\n".join(issues),
                )

    def test_direct_stata_batch_invocation_is_detected(self) -> None:
        for relative in ("README.md", "AGENTS.md"):
            with self.subTest(relative=relative), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                copy_identity_surface(root)
                document = root / relative
                document.write_text(
                    document.read_text(encoding="utf-8")
                    + '\nstata-mp -b do "Hypercapnia Case Definitions.do"\n',
                    encoding="utf-8",
                )
                issues = audit.validate_identity(root)
                self.assertTrue(
                    any(
                        issue.startswith(f"{relative}: direct Stata batch command")
                        for issue in issues
                    ),
                    "\n".join(issues),
                )

    def test_provenance_tokens_are_required_in_public_docs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            copy_identity_surface(root)
            reproducibility = root / "docs/REPRODUCIBILITY.md"
            reproducibility.write_text(
                reproducibility.read_text(encoding="utf-8").replace(
                    audit.APPROVED_INPUT_SCHEMA,
                    "removed-schema-token",
                ),
                encoding="utf-8",
            )
            issues = audit.validate_identity(root)
        self.assertTrue(
            any(
                "docs/REPRODUCIBILITY.md: missing approved input token" in issue
                for issue in issues
            ),
            "\n".join(issues),
        )

    def test_dictionary_output_row_is_detected(self) -> None:
        fieldnames = [
            "variable_name",
            "source_or_origin",
            "workflow_role",
            "upstream_derivation_status",
            "review_status",
        ]
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "dictionary.csv"
            with path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerow(
                    {
                        "variable_name": "result.png",
                        "source_or_origin": "generated output",
                        "workflow_role": "derived_analysis",
                        "upstream_derivation_status": "not_applicable",
                        "review_status": "draft",
                    }
                )
            with path.open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
        issues = audit.validate_dictionary_rows(rows)
        self.assertTrue(any("output artifact" in issue for issue in issues))

    def test_unselected_full_db_derivation_cannot_be_silently_verified(self) -> None:
        with (ROOT / "data_dictionary.csv").open(
            newline="",
            encoding="utf-8",
        ) as handle:
            rows = list(csv.DictReader(handle))
        target = next(row for row in rows if row["variable_name"] == "paco2")
        target["upstream_derivation_status"] = "verified"
        issues = audit.validate_dictionary_rows(rows)
        self.assertTrue(
            any("input derivation must be blocked" in issue for issue in issues),
            "\n".join(issues),
        )


class RepositoryContractTests(unittest.TestCase):
    def test_consort_diagram_uses_portable_tiff_conversion(self) -> None:
        notebook = json.loads(
            (ROOT / "Case Definitions Consort.ipynb").read_text(encoding="utf-8")
        )
        source = "\n".join(
            "".join(cell.get("source", []))
            for cell in notebook["cells"]
            if cell.get("cell_type") == "code"
        )

        self.assertIn("from PIL import Image", source)
        self.assertIn("format='png'", source)
        self.assertNotIn("format='tiff'", source)
        self.assertIn("format='TIFF'", source)
        self.assertIn("dpi=(300, 300)", source)
        self.assertIn("compression='tiff_lzw'", source)
        self.assertIn("temporary_tiff.replace(tiff_path)", source)
        self.assertIn("rendered_png.unlink(missing_ok=True)", source)

    def test_repository_contract(self) -> None:
        issues = audit.find_issues(ROOT)
        self.assertEqual([], issues, "\n".join(issues))


if __name__ == "__main__":
    unittest.main()
