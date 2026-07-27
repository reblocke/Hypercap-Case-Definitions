# Restricted Validation Evidence

## Purpose and Boundary

This document summarizes value-suppressed validation of the guarded Stata
workflow and the scientific-alignment update. It is a technical appendix, not a
substitute for the article or a public reproduction of its results.

The summary omits restricted input locations and hashes, dataset dimensions,
result values, logs, workbook contents, images, and detailed comparison
reports. Those artifacts remain local under ignored `outputs/` paths.

The validation does not establish that the available file was the exact
publication input, reproduce upstream TriNetX construction, or resolve the
scientific questions that remain open in
[`SCIENTIFIC_ALIGNMENT.md`](SCIENTIFIC_ALIGNMENT.md).

## Releases Represented

| Release | Validation interpretation |
| --- | --- |
| `v1.0.0` | Paper-associated historical implementation. The preserved pre-correction scientific body was compared with a legacy baseline, and two guarded candidate runs were repeatable. |
| `v2.0.0` | Current LLM-assisted update. The corrected scientific implementation was run twice and compared in correction mode; subsequent repository changes did not alter the test-locked scientific body. |

The historical implementation is also preserved as annotated tag
`hcd-000b-historical`. Both that tag and `v1.0.0` resolve to
`e5671dc35a30ba433588f0793de48d0e97a6d3df`.

## Validated Commits

| Role | Commit |
| --- | --- |
| Legacy baseline | `3520a95663d69beeebe48ad8948300506f4f10de` |
| Initial guarded candidate | `ba7337bdb24494859e668a26383d91202eb9f76d` |
| Comparator and input-import remediation candidate | `17ffcf9c16ff11254717fe0d8b234bfa9b79a7c1` |
| Input-approval contract candidate | `85076dc76c4de246c7f390c71d91ce54e8a16750` |
| Evidence-integrity candidate | `ccbea198694b2e6c97c1e233686883c7f6742feb` |
| Scientific-alignment candidate | `9f8fecd2384edd3ae3fd4b0c8ecbc8536769bc19` |

The restricted artifact was observed in upstream checkout
`1185a6bc9957a02cb24be5f1f7fa10c48d8a4c13`. The public provenance record
associates the approved input with producer commit
`44f49748d415e92b7d50b50d86b8fdea29f6cb07` and observed schema
`hypercapnia-full-db-v1`.

This historical association did not independently rebuild the upstream
dataset. The historical build did not preserve source-file hashes or
clean-worktree state.

## Observed Environment

- Stata 18 IC on macOS, Apple Silicon;
- Python 3.11.15 for the runner, comparator, and public checks; and
- Graphviz 14.1.5 for the data-free diagram smoke test.

Publication-time versions of community Stata packages remain unresolved.
Dependencies required by the current analysis passed the runtime preflight;
their observed details remain in local run manifests.

## Evidence Summary

| Check | Outcome |
| --- | --- |
| Public tests, public-surface audit, and citation validation | Pass |
| Data-free diagram notebook smoke test | Pass |
| Restricted input contract and adjacent approval checks | Pass |
| Legacy baseline and guarded candidate completion | Pass |
| Two independent candidate executions for each guarded-workflow milestone | Pass |
| Zero Stata process exits and fresh completion controls | Pass |
| Complete expected artifact inventory | Pass |
| Current input identity before, across, and after comparison | Pass |
| Contained, distinct artifact roots and distinct candidate run identifiers | Pass |
| Semantic workbook comparison | Pass |
| Decoded image-pixel comparison | Pass |
| Required Stata graph and copied-analysis-file integrity | Pass |
| Normalized substantive analysis-log comparison | Pass |
| Candidate run-to-run repeatability | Pass |
| Scientific-correction impact classification | Pass with documented changes |
| Relocated preserved-run comparison | Pass |

The comparator normalizes only known wrapper and path differences. It does not
suppress substantive analysis commands or results. Workbook cells and decoded
image pixels are compared directly rather than inferred from log equality.

## Historical Implementation

The historical validation used a clean legacy baseline and two clean guarded
candidate runs against the same restricted input. The scientific analysis body
was test-locked to the baseline implementation except for removal of a
row-level listing command.

The preserved baseline matched the first candidate. The candidates matched each
other, completed with the expected artifacts, and retained the same input
identity through comparison.

Subsequent guarded-workflow milestones were each executed twice from clean,
isolated checkouts. Those runs confirmed:

- failed input-contract imports retain their nonzero status;
- a run cannot pass after a nonzero Stata process exit;
- the adjacent input approval is validated before and after execution;
- candidate run identifiers and artifact roots are independent;
- completion controls are checked from disk rather than trusted from stale
  manifests; and
- comparison reports suppress result values and row-level log content.

## Current Scientific-Alignment Update

The scientific-alignment candidate
`9f8fecd2384edd3ae3fd4b0c8ecbc8536769bc19` was executed twice from separate
clean checkouts against the same approved restricted input. Both runs:

- exited zero;
- passed input, approval, status, completion, and artifact checks;
- confirmed that the transparent Bülbül and Meservey rules matched their
  verified upstream aggregate flags; and
- matched for every compared workbook cell, image pixel, graph-presence check,
  normalized log, dependency record, input reference, and copied analysis
  file.

Correction mode classified the historical comparison as changed while
requiring exact candidate repeatability. The changed artifact groups included:

- Calvo and Wilson definition summaries;
- overall, testing-strategy, and regional agreement heatmaps;
- workup and location summaries; and
- the normalized analysis transcript.

The Bülbül and Meservey summaries were unchanged because the transparent rules
matched their verified upstream aggregates for this input. Vonderbank and Chung
boundary corrections did not change their summaries for this input.

All seven spline images changed relative to the historical baseline even though
the corresponding fitted-model and prediction transcript segments were
unchanged. The two corrected runs produced identical spline pixels. This
supports repeatability of the observed corrected output but does not resolve
SA-007, the known row-selection sensitivity in the plotting code.

## Interpretation

The evidence supports two narrow conclusions:

1. the paper-associated historical implementation preserved the observed
   baseline downstream results for the restricted input available during
   validation; and
2. the current corrected scientific implementation produced repeatable
   candidate outputs in the observed environment.

A correction-mode pass records historical impact and candidate repeatability.
It is not an equivalence claim between the current and historical releases.

The remaining SA-002, SA-007, SA-008, and SA-009 questions are unresolved.
Detailed manifests, dependency records, comparison reports, and generated
results remain local and are not release assets.
