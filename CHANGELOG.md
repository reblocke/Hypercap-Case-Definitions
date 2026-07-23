# Changelog

## Unreleased

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
