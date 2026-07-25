# Restricted Validation Summary

## Scope

This summary records the HCD-000B downstream validation completed on
2026-07-22 Mountain Time (2026-07-23 UTC), the review-remediation validation
completed on 2026-07-23 Mountain Time, and the approved-input and
evidence-hardening validations completed on 2026-07-24 Mountain Time. It also
records the HCD-001 scientific-correction validation completed on 2026-07-24
Mountain Time. It intentionally omits restricted input locations and hashes,
dataset dimensions, result values, logs, and detailed comparison evidence.

The validated analysis commits were:

- legacy baseline: `3520a95663d69beeebe48ad8948300506f4f10de`;
- initial guarded candidate: `ba7337bdb24494859e668a26383d91202eb9f76d`;
- review-remediation candidate:
  `17ffcf9c16ff11254717fe0d8b234bfa9b79a7c1`;
- approved-input candidate:
  `85076dc76c4de246c7f390c71d91ce54e8a16750`;
- evidence-hardening candidate:
  `ccbea198694b2e6c97c1e233686883c7f6742feb`;
- HCD-001 scientific-correction candidate:
  `9f8fecd2384edd3ae3fd4b0c8ecbc8536769bc19`.

The restricted artifact was observed in upstream checkout
`1185a6bc9957a02cb24be5f1f7fa10c48d8a4c13`. Subsequent owner adjudication
bound the approved restricted input to producer commit
`44f49748d415e92b7d50b50d86b8fdea29f6cb07` and observed schema
`hypercapnia-full-db-v1`. This historical assignment did not independently
rebuild the upstream dataset; source-file hashes and clean-worktree state were
not preserved.

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
| Two clean review-remediation candidate executions and complete output gates | Pass |
| Adjacent input approval before and after each new run | Pass |
| Two clean approved-input candidate executions and complete output gates | Pass |
| Two clean evidence-hardening candidate executions and zero process exits | Pass |
| Two clean HCD-001 candidate executions and complete output gates | Pass |
| HCD-001 aggregate-definition cross-checks in both runs | Pass |
| Live artifact-root, run-ID, and completion-control comparison gates | Pass |
| Input identity before, across, and after the three runs | Pass |
| HCD-000B baseline versus candidate semantic workbook comparison | Pass |
| HCD-000B baseline versus candidate decoded PNG-pixel comparison | Pass |
| Required Stata graph presence and copied-do integrity | Pass |
| Path- and wrapper-normalized analysis transcript comparison | Pass |
| Candidate run-to-run repeatability | Pass |
| HCD-001 correction-impact classification | Pass with documented changes |

The transcript comparison removes the legacy row listing, generated output
paths, Stata line wrapping caused by path length, prompt-only blank lines, and
the neutral wrapper tail. It does not suppress substantive analysis commands
or results. Workbook cells and structure and decoded image pixels are compared
directly rather than inferred from transcript equality.

## Review-Finding Re-adjudication

On 2026-07-23 Mountain Time, the preserved three-run evidence was re-compared
after hardening the comparator. The re-adjudication intentionally omitted the
optional expected-input-hash argument so that the default documented command
exercised the corrected behavior. It confirmed:

- the supplied baseline was a legacy-interface run;
- both candidates were guarded runs from the same nonempty candidate commit;
- the current restricted input still matched the input hash in every run
  manifest;
- baseline equivalence passed; and
- candidate repeatability passed.

A targeted Stata negative-path smoke also confirmed that a failed contract
import returns its original nonzero status after frame cleanup. A separate
positive-path preflight against the current restricted input passed.

The exact review-remediation commit
`17ffcf9c16ff11254717fe0d8b234bfa9b79a7c1` was then executed twice from a
clean isolated checkout against the same restricted input. Both guarded runs
passed fresh Stata status, completion, and 40-artifact inventory gates. The
strengthened comparator confirmed that the current input matched all three run
manifests, the preserved legacy baseline matched the first remediation run, and
the two remediation runs were repeatable. Detailed run and comparison evidence
remains local and ignored.

## Approved-Input Validation

The exact approved-input commit
`85076dc76c4de246c7f390c71d91ce54e8a16750` was executed twice from a clean
isolated checkout against the same restricted input. The adjacent local
approval was created through `make input-approve` and bound the input to the
public producer/schema record and current data-dictionary contract without
recording a local path or row-level value.

Both schema-v2 guarded runs passed pre-launch and post-analysis approval
validation, fresh Stata status and completion checks, and the complete
40-artifact inventory. Their recorded approval references matched the current
manifest. The schema-v2 comparator validated the current approval before
loading the runs and again after comparison, accepted the preserved historical
baseline without a newer approval reference, and confirmed:

- current input identity across the baseline and both candidates;
- baseline equivalence for the first approved-input run; and
- candidate run-to-run repeatability.

The comparison report contained no failures. Detailed manifests, hashes,
artifacts, and comparison evidence remain ignored and local.

## Evidence-Hardening Validation

The exact evidence-hardening commit
`ccbea198694b2e6c97c1e233686883c7f6742feb` was executed twice from a clean
isolated checkout against the same approved restricted input. The guarded runs
used distinct run identifiers and directories. Both Stata processes exited
zero, and both runs passed fresh Stata status, approval revalidation, all live
completion controls, and the complete 40-artifact inventory.

The strengthened schema-v2 comparator rejected neither run during its new
precomparison gates: every artifact root was contained within its run
directory and distinct from the other artifact roots, the candidate run
identifiers were distinct, and the required controls remained present on disk.
The preserved legacy baseline matched the first evidence-hardening run, the
two evidence-hardening runs were repeatable, and the current approved input
remained unchanged. The final comparison report contained no failures;
detailed evidence remains ignored and local.

## HCD-001 Scientific-Correction Validation

The pre-correction public implementation was preserved as annotated tag
`hcd-000b-historical`. The exact HCD-001 candidate commit
`9f8fecd2384edd3ae3fd4b0c8ecbc8536769bc19` was then executed twice from
separate clean isolated checkouts against the same approved restricted input,
using distinct run identifiers and directories.

Both Stata processes exited zero. Both runs passed the fresh status,
completion, approval-revalidation, input-contract, and complete 40-artifact
gates. In both runs, the guarded contract confirmed that the transparent
Bülbül and Meservey simulated definitions matched their verified upstream
aggregate flags. The two candidate runs matched for every semantic workbook
cell, decoded PNG pixel, required graph presence check, normalized analysis
log, dependency record, input reference, and copied analysis file.

Correction mode classified the historical comparison as changed while
requiring exact candidate repeatability and all evidence-integrity checks. The
changed artifact set included:

- the Calvo and Wilson definition-summary workbooks;
- the four overall definition-overlap heatmaps;
- the three mutually exclusive testing-strategy heatmaps;
- the workup and location summary workbooks;
- the four regional definition-overlap heatmaps; and
- the normalized analysis transcript.

The Bülbül and Meservey summary workbooks remained unchanged, consistent with
their new transparent rules matching the verified upstream aggregates. The
Vonderbank and Chung boundary corrections did not change their summary
workbooks for this input.

All seven spline PNGs also changed relative to the historical baseline. The
corresponding fitted-model and prediction-output transcript segments were
identical, and the two corrected runs produced identical spline pixels. This
is consistent with the pre-existing SA-007 concern: the plot filters on the
original `paco2_rounded` field after `xblc` creates the prediction grid, so
earlier row-order changes can affect which predicted points are drawn. SA-007
therefore remains unresolved and no unapproved plotting change was made.
Detailed cell locations, hashes, images, logs, and run evidence remain ignored
and local.

## Interpretation and Boundary

These checks support that HCD-000B preserves the baseline downstream results
for the restricted `full_db.dta` available during validation and that two
independent candidate runs are repeatable in the observed environment. The
scientific analysis body is test-locked to HCD-000A except for removal of the
row-level listing command. The approved-input validation additionally supports
that the exact `85076dc` infrastructure commit enforces the recorded local
input approval before and after execution and comparison. The
evidence-hardening validation supports that the exact `ccbea19` implementation
also enforces zero Stata process exits, distinct run evidence, contained
artifact roots, and live completion controls.

The HCD-001 validation supports that the exact `9f8fecd` implementation applies
the owner-approved SA-001, SA-003, SA-004, SA-005, SA-006, SA-010, and SA-011
decisions and is repeatable in the observed environment. Correction-mode
`PASS` means that the historical impact was recorded and the corrected
candidates were repeatable; it does not mean that the corrected artifacts are
equivalent to the historical baseline.

This is not proof that the available file is the exact publication input, does
not reproduce upstream TriNetX construction, and does not resolve the
remaining SA-002, SA-007, SA-008, or SA-009 scientific-alignment items
documented in `SCIENTIFIC_ALIGNMENT.md`. Detailed run manifests, dependency
records, comparison reports, and generated outputs remain local under ignored
`outputs/` paths.
