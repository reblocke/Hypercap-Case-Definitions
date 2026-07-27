# Repository Instructions for Coding Agents

## Purpose

This public repository contains Stata analysis code and a Python/Graphviz cohort
diagram notebook for the *CHEST* article “The Consistency of Hypercapnic
Respiratory Failure Case Definitions in Electronic Health Record Data.”

## Public and Data-Safety Rules

- Do not add PHI, restricted TriNetX datasets, row-level derivatives,
  credentials, machine-specific paths, private drafts, or
  publisher-formatted article files.
- Treat `data/private/full_db.dta`, any `Data/` directory, and generated `.dta`
  files as local-only restricted artifacts.
- Keep generated outputs under ignored `outputs/` paths.
- Do not attach restricted data, generated results, logs, or local manifests to
  public releases.
- Link to DOI, PubMed, and PMC records instead of copying article text.

## Orientation

- Start with `README.md` for the public overview and supported commands.
- Use `llms.txt` for the concise machine-readable project index.
- Read `data_dictionary.md` and `data_dictionary.csv` before changing variables
  or case-definition logic.
- Read `docs/SCIENTIFIC_ALIGNMENT.md` before changing scientific behavior.
- Use `CITATION.cff` for structured citation and release metadata.

## Releases

- `v1.0.0` is the paper-associated historical implementation.
- `v2.0.0` is the current LLM-assisted reproducibility and
  scientific-alignment update.
- Do not move, recreate, or delete release tags.
- Do not describe the public checks as reproduction of the restricted analysis
  or article estimates.

## Supported Workflow

Approve an eligible restricted input and run the Stata analysis from the
repository root:

```bash
make input-approve INPUT_ROOT="/approved/restricted/hypercapnia" APPROVE_RESTRICTED_INPUT=YES
make stata-run STATA_BIN="/path/to/stata" INPUT_ROOT="/approved/restricted/hypercapnia"
```

`INPUT_ROOT` contains `full_db.dta` and its ignored adjacent
`full_db.manifest.json`. `OUTPUT_ROOT` defaults to `outputs/stata`.
`make stata-run` is the sole supported scientific execution interface.

For controlled validation, reuse the preserved legacy baseline and run the
candidate twice in clean isolated checkouts:

```bash
make stata-compare BASELINE_RUN="..." CANDIDATE_RUN_1="..." CANDIDATE_RUN_2="..." INPUT_ROOT="..."
```

Keep manifests, hashes, comparison reports, and generated results in ignored
`outputs/` paths. Public validation summaries may report commits, environments,
checks, and pass/fail outcomes, but not result values or row-level content.

Public data-free verification:

```bash
python3 -m pip install --require-hashes -r requirements.txt
make check
make diagram-smoke
```

The diagram workflow requires the Python `graphviz` package and the system
Graphviz `dot` executable.

## Change Discipline

- Preserve scientific logic unless a change is explicitly authorized and
  validated.
- Preserve unresolved scientific-alignment items rather than inferring a
  resolution.
- Preserve machine-readable approval enums, input contracts, and dependency
  requirements.
- Keep path handling repository-relative and argument-driven.
- Update the data dictionary when adding or renaming variables.
- Avoid unrelated refactors, formatting churn, and speculative abstractions.

## Verification Before Publication

- Run `git diff --check`.
- Run `make check` with the locked Python 3.11 environment.
- Run `make diagram-smoke` when Python 3.11 and Graphviz are available.
- Validate `CITATION.cff` with `cffconvert`.
- Confirm no user-home paths, Windows local paths, restricted artifacts, or
  generated root-level files are tracked.
- Confirm the guarded runner cannot report success without an approved,
  unchanged input, a zero Stata process exit, fresh status and completion
  controls, and the complete required artifact inventory.
- Confirm the comparator requires independent candidate evidence, contained
  artifact roots, live completion controls, and current input identity.
- Confirm README, `llms.txt`, `CITATION.cff`, and the data dictionary agree on
  DOI `10.1016/j.chest.2025.08.002`, PMID `40885535`, PMCID `PMC12739763`, and
  the restricted TriNetX boundary.
- If Stata or the restricted input is unavailable, document the skipped
  restricted run rather than creating synthetic patient-like data.
