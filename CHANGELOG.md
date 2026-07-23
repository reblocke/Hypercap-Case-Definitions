# Changelog

## Unreleased

- Added a guarded restricted-data Stata runner with unique run identifiers,
  preflight dependency checks, data-dictionary-driven input validation,
  provenance manifests, explicit completion status, and a 40-artifact success
  gate.
- Added a value-suppressing comparator for one legacy baseline and two
  independent candidate runs, including semantic workbook, decoded-pixel,
  normalized-log, environment, dependency, and repeatability checks.
- Removed the only row-level `list` output from the analysis log while locking
  the remaining scientific analysis body against unintentional changes.
- Expanded the locked Python environment for workbook and image comparison and
  documented the restricted validation protocol.
- Completed a clean legacy-baseline run and two clean candidate runs against
  the available restricted input; baseline equivalence, candidate
  repeatability, output completeness, and post-run input identity passed.
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
