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
stata-mp -b do "Hypercapnia Case Definitions.do" "data/private" "outputs/stata"
```

The first argument is the input root containing `full_db.dta`; the second argument is the output root. The script should fail clearly if the restricted input is absent.

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
- Confirm README, `llms.txt`, `CITATION.cff`, and the data dictionary agree on DOI `10.1016/j.chest.2025.08.002`, PMID `40885535`, PMCID `PMC12739763`, and the restricted TriNetX data boundary.
- If Stata is unavailable or the restricted data are absent, document the skipped smoke check rather than creating synthetic patient-like data.
