from __future__ import annotations

import csv
import sys
import tempfile
import unittest
from pathlib import Path, PurePosixPath

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import check_public_surface as audit  # noqa: E402


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
        self.assertTrue(any("must remain unapproved" in issue for issue in issues))

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


class RepositoryContractTests(unittest.TestCase):
    def test_repository_contract(self) -> None:
        issues = audit.find_issues(ROOT)
        self.assertEqual([], issues, "\n".join(issues))


if __name__ == "__main__":
    unittest.main()
