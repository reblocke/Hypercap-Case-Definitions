# Restricted Validation Summary

## Scope

This summary records the HCD-000B downstream validation completed on
2026-07-22 Mountain Time (2026-07-23 UTC). It intentionally omits restricted
input locations and hashes, dataset dimensions, result values, logs, and
detailed comparison evidence.

The validated analysis commits were:

- legacy baseline: `3520a95663d69beeebe48ad8948300506f4f10de`;
- guarded candidate: `ba7337bdb24494859e668a26383d91202eb9f76d`.

The restricted artifact was observed in upstream checkout
`1185a6bc9957a02cb24be5f1f7fa10c48d8a4c13`. That checkout is context only:
the exact commit that produced the article input remains `UNRESOLVED`.

## Environment

- Stata 18 IC on macOS, Apple Silicon;
- Python 3.11.15 for the runner, comparison, and public checks;
- Graphviz 14.1.5 for the data-free diagram smoke test.

Publication-time Stata package versions remain unresolved. All dependencies
required by the current analysis passed the runtime preflight, and their local
hashes and observed version headers remain in ignored run manifests.

## Validation Results

| Check | Outcome |
| --- | --- |
| Public repository tests, safety audit, and CFF validation | Pass |
| Data-free diagram notebook smoke test | Pass |
| Restricted input contract preflight | Pass |
| Clean legacy baseline execution and complete output gate | Pass |
| Clean candidate execution 1 and complete output gate | Pass |
| Clean candidate execution 2 and complete output gate | Pass |
| Input identity before, across, and after the three runs | Pass |
| Baseline versus candidate semantic workbook comparison | Pass |
| Baseline versus candidate decoded PNG-pixel comparison | Pass |
| Required Stata graph presence and copied-do integrity | Pass |
| Path- and wrapper-normalized analysis transcript comparison | Pass |
| Candidate run-to-run repeatability | Pass |

The transcript comparison removes the legacy row listing, generated output
paths, Stata line wrapping caused by path length, prompt-only blank lines, and
the neutral wrapper tail. It does not suppress substantive analysis commands
or results. Workbook cells and structure and decoded image pixels are compared
directly rather than inferred from transcript equality.

## Interpretation and Boundary

These checks support that HCD-000B preserves the baseline downstream results
for the restricted `full_db.dta` available during validation and that two
independent candidate runs are repeatable in the observed environment. The
scientific analysis body is test-locked to HCD-000A except for removal of the
row-level listing command.

This is not proof that the available file is the exact publication input, does
not reproduce upstream TriNetX construction, and does not resolve the
scientific-alignment items documented in `SCIENTIFIC_ALIGNMENT.md`. Detailed
run manifests, dependency records, comparison reports, and generated outputs
remain local under ignored `outputs/` paths.
