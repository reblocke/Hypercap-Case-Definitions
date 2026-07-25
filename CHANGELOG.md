# Changelog

## Unreleased

- Applied the owner-approved HCD-001 decisions for SA-001, SA-003 through
  SA-006, SA-010, and SA-011 while preserving SA-002, SA-007, SA-008, and
  SA-009 as unresolved.
- Executed the exact HCD-001 scientific-correction commit twice against the
  same approved input, recorded correction impact, and confirmed complete
  guarded runs and candidate repeatability.
- Made preserved-run log comparison relocation-safe by normalizing only Stata
  file notifications that share one root for the known expected artifact
  inventory; inconsistent roots, unknown paths, and substantive transcript
  differences remain comparison failures.
- Made CONSORT TIFF export portable by rendering an intermediate PNG,
  converting it with the locked Pillow dependency, validating the TIFF, and
  removing intermediate renderer files.
- Added explicit one-time approval for the adjacent restricted input, bound new
  runs and comparisons to that manifest, and made `make stata-run` the sole
  supported scientific execution interface.
- Recorded the owner-approved historical producer commit and observed input
  schema while retaining the missing upstream source-hash, clean-worktree, and
  variable-derivation limitations.
- Executed the exact approved-input implementation twice from a clean isolated
  checkout and confirmed approval-gated completion, baseline equivalence,
  candidate repeatability, and post-comparison input identity.
- Added a guarded restricted-data Stata runner with unique run identifiers,
  preflight dependency checks, data-dictionary-driven input validation,
  provenance manifests, explicit completion status, and a 40-artifact success
  gate.
- Added a value-suppressing comparator for one legacy baseline and two
  independent candidate runs, including semantic workbook, decoded-pixel,
  normalized-log, environment, dependency, and repeatability checks.
- Hardened comparison adjudication so current input drift, incorrect run roles,
  missing or mismatched candidate commits, and omission of the current input
  fail before equivalence or repeatability labels are assigned.
- Required contained and distinct artifact roots, distinct candidate run IDs,
  live completion controls, and a zero Stata process return code before a
  validation run or comparison can pass.
- Executed the exact evidence-hardening implementation twice from a clean
  isolated checkout and confirmed zero process exits, complete live controls,
  baseline equivalence, and independent candidate repeatability.
- Preserved Stata contract-import failures through cleanup, separately locked
  the diagnostic-performance program, and made phenotype source locations
  executable-code checks rather than unchecked line-number annotations.
- Removed the only row-level `list` output from the analysis log while locking
  the remaining scientific analysis body against unintentional changes.
- Expanded the locked Python environment for workbook and image comparison and
  documented the restricted validation protocol.
- Completed a clean legacy-baseline run and two clean candidate runs against
  the available restricted input; baseline equivalence, candidate
  repeatability, output completeness, and post-run input identity passed.
- Re-adjudicated those preserved runs with the strengthened comparator and
  confirmed the corrected contract-import failure path plus the current
  restricted-input preflight without repeating unchanged scientific analyses.
- Executed the exact review-remediation commit twice from a clean isolated
  checkout and confirmed complete guarded runs, baseline equivalence,
  repeatability, and post-run input identity with the strengthened comparator.
- Defined public, restricted downstream, and upstream reproducibility boundaries.
- Added code-derived phenotype, upstream provenance, Stata dependency, and generated-output metadata without approving unresolved scientific rules.
- Normalized the data dictionary to one row per input or derived variable and separated generated artifacts into an output manifest.
- Added a locked Python 3.11 environment, data-free repository checks, explicit notebook execution, and hosted CI.
- Reframed repository documentation around the final CHEST article metadata, PubMed record, and PMC/NLM full-text pointer.
- Added `llms.txt`, repository-specific `AGENTS.md`, data dictionary files, `.gitignore`, and notebook requirements.
- Made the Stata workflow repo-root runnable with argument-driven restricted-input and output roots.
- Moved generated Stata logs, copied do-files, figures, tables, and temporary graph files under ignored `outputs/` paths.
- Corrected the Cavalot value-label assignment from `def8` to `def9`.
- Corrected the Midwest regional spline block to filter `location == 2` instead of duplicating the West filter.
