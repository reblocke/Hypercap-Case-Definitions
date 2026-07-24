# AGENTS

## Project Purpose

This public repository contains Stata analysis code and a Python/Graphviz diagram notebook for the CHEST article "The Consistency of Hypercapnic Respiratory Failure Case Definitions in Electronic Health Record Data."

## Public and Data-Safety Rules

- Do not add PHI, restricted TriNetX datasets, row-level derived data, credentials, local paths, private drafts, or publisher-formatted article files.
- Treat `data/private/full_db.dta`, any `Data/` directory, and generated `.dta` files as local-only restricted artifacts.
- Link DOI, PubMed, and PMC/NLM records instead of copying article text into repository docs.
- Keep generated outputs under ignored `outputs/` paths unless an explicit release workflow says otherwise.

## How to Orient Quickly

- Start with `README.md` for the human-facing overview and run commands.
- Use `llms.txt` for the concise machine-readable project index.
- Use `data_dictionary.md` and `data_dictionary.csv` before changing variable names, labels, or case-definition logic.
- Use `CITATION.cff` for structured citation metadata.

## Workflow

Canonical Stata command from the repository root:

```bash
make input-approve INPUT_ROOT="/approved/restricted/hypercapnia" APPROVE_RESTRICTED_INPUT=YES
make stata-run STATA_BIN="/path/to/stata" INPUT_ROOT="/approved/restricted/hypercapnia"
```

`INPUT_ROOT` contains `full_db.dta` and its ignored adjacent
`full_db.manifest.json`; `OUTPUT_ROOT` defaults to `outputs/stata`.
`make stata-run` is the sole supported scientific execution interface. The
guarded runner must fail clearly for a missing or unapproved input, approval
drift, dependency, unsafe or colliding run ID, failed input contract, failed
Stata status, or incomplete legacy artifact inventory.

For equivalence validation, reuse the preserved legacy baseline and run the
guarded candidate twice in clean isolated checkouts, then invoke:

```bash
make stata-compare BASELINE_RUN="..." CANDIDATE_RUN_1="..." CANDIDATE_RUN_2="..." INPUT_ROOT="..."
```

Keep manifests, dependency hashes, input hashes, comparison details, and all
generated results under ignored `outputs/` paths. A public validation summary
may name commits, environment, checks, and pass/fail status only.

Optional CONSORT notebook workflow:

```bash
python3 -m pip install --require-hashes -r requirements.txt
make diagram-smoke
```

The notebook requires the Python `graphviz` package and the system Graphviz `dot` binary.

Public data-free verification:

```bash
make check
make diagram-smoke
```

Passing public checks does not establish that the restricted Stata analysis or
article estimates were reproduced.

## Verification Before Publishing Changes

- Run `git diff --check`.
- Run `make check`.
- Run `make diagram-smoke` when Python 3.11 and Graphviz are available.
- Validate `CITATION.cff` as YAML and, when available, with `cffconvert`.
- Confirm no hard-coded user-home paths, Windows local paths, legacy generated-output roots, root `.gph`, root `.log`, or restricted-data paths are introduced.
- Confirm the guarded runner is collision-safe and that `SUCCESS` cannot be written without fresh Stata status, completion, and all 40 legacy artifacts.
- Confirm new runs cannot resolve or launch Stata without a valid adjacent input approval matching the tracked producer/schema and data dictionary.
- Confirm the runner revalidates that approval before `SUCCESS`, the comparator revalidates it after comparisons, and filesystem aliases cannot satisfy distinct run roles.
- Confirm the runner rejects a split runner/analysis checkout before creating a run directory.
- Confirm comparison reports never expose workbook values or row-level log content.
- Confirm README, `llms.txt`, `CITATION.cff`, and the data dictionary agree on DOI `10.1016/j.chest.2025.08.002`, PMID `40885535`, PMCID `PMC12739763`, and the restricted TriNetX data boundary.
- If Stata is unavailable or the restricted data are absent, document the skipped smoke check rather than creating synthetic patient-like data.
