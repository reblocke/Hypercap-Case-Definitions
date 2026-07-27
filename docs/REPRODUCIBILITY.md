# Reproducibility Guide

This repository accompanies:

> Locke BW, et al. The Consistency of Hypercapnic Respiratory Failure Case
> Definitions in Electronic Health Record Data. *CHEST*. 2026;169(1):230-243.
> doi: [10.1016/j.chest.2025.08.002](https://doi.org/10.1016/j.chest.2025.08.002).

The repository supports three distinct reproducibility layers. A successful
check at one layer does not establish success at another.

## Public, Data-Free Checks

A public clone can:

- validate repository metadata and data-safety rules;
- check the documented phenotype, dependency, input, and output contracts;
- confirm that the tracked notebook contains no execution state or embedded
  output; and
- render the cohort diagram from fixed published aggregate counts.

Use Python 3.11:

```bash
python3 -m pip install --require-hashes -r requirements.txt
make check
make diagram-smoke
```

The diagram smoke test also requires the system Graphviz `dot` executable.
Generated notebooks and figures are written under ignored `outputs/` paths.
The notebook does not calculate its counts from patient-level data or Stata
results.

These checks establish that the public repository is internally consistent and
data-free. They do not reproduce the restricted analysis or validate article
estimates.

## Restricted Downstream Analysis

### Required Input

The Stata workflow starts with an analysis-ready, restricted encounter-level
file named `full_db.dta`.

The approved provenance record identifies:

- upstream repository:
  `https://github.com/reblocke/trinetx-hypercapnia-code`;
- producer commit:
  `44f49748d415e92b7d50b50d86b8fdea29f6cb07`; and
- observed schema: `hypercapnia-full-db-v1`.

This assignment was adjudicated from historical repository evidence. The
historical build did not preserve source-file hashes or clean-worktree state,
and this repository does not reproduce that build. The selected upstream
derivations of `hypercap_on_abg` and `hypercap_resp_failure` were separately
verified from the producer-commit Git object. Other upstream derivations remain
outside the verified scope.

The complete public provenance record is
[`metadata/upstream_dependency.yml`](../metadata/upstream_dependency.yml).

### Input Approval

Create an adjacent local approval manifest before the first run:

```bash
make input-approve \
  INPUT_ROOT="/approved/restricted/hypercapnia" \
  APPROVE_RESTRICTED_INPUT=YES
```

The command writes ignored `full_db.manifest.json` beside `full_db.dta`. The
manifest binds the input bytes to the approved producer, schema, and current
data-dictionary contract. It contains no local path, row count, or patient-level
value. Replacing an existing approval also requires
`REPLACE_APPROVED_INPUT=YES`.

### Guarded Stata Run

Run the analysis from the repository root:

```bash
make stata-run \
  STATA_BIN="/path/to/stata" \
  INPUT_ROOT="/approved/restricted/hypercapnia" \
  OUTPUT_ROOT="outputs/stata"
```

`make stata-run` is the sole supported scientific execution interface.

The guarded runner:

- rejects missing, unapproved, changed, or contract-incompatible input;
- rejects unsafe or colliding run identifiers and split analysis checkouts;
- records the analysis commit, dirty state, input and code hashes, environment,
  and dependency hashes;
- checks required community Stata commands before analysis;
- requires a zero Stata process exit, fresh successful status, explicit
  completion, and the complete expected artifact inventory;
- revalidates the input and approval after analysis; and
- writes `SUCCESS` only after every required condition passes.

Every run receives a unique ignored directory under `outputs/stata/`. Run
manifests, logs, copied analysis files, tables, figures, and temporary Stata
graphs remain local.

## Restricted Validation Protocol

Use this protocol when assessing a change to the scientific analysis or guarded
execution workflow:

1. Use clean, isolated checkouts for the preserved baseline and candidate.
2. Validate the same adjacent restricted-input approval before every run.
3. Reuse the preserved historical baseline.
4. Run the candidate twice with distinct run identifiers and directories.
5. Compare the baseline and both candidates:

   ```bash
   make stata-compare \
     BASELINE_RUN="outputs/validation/baseline" \
     CANDIDATE_RUN_1="outputs/validation/candidate-1" \
     CANDIDATE_RUN_2="outputs/validation/candidate-2" \
     INPUT_ROOT="/approved/restricted/hypercapnia"
   ```

For an explicitly adjudicated scientific correction, add
`COMPARISON_MODE=correction`:

```bash
make stata-compare \
  COMPARISON_MODE=correction \
  BASELINE_RUN="outputs/validation/baseline" \
  CANDIDATE_RUN_1="outputs/validation/candidate-1" \
  CANDIDATE_RUN_2="outputs/validation/candidate-2" \
  INPUT_ROOT="/approved/restricted/hypercapnia"
```

Before reporting a result, the comparator verifies:

- one legacy baseline and two guarded candidates;
- matching nonempty candidate commits and distinct candidate run identifiers;
- contained, distinct artifact roots and live completion controls;
- current input and approval identity across all runs;
- environment and dependency compatibility;
- semantic workbook content, decoded image pixels, required Stata graph files,
  copied analysis files, and normalized analysis logs; and
- candidate repeatability.

The comparator reports discrepancy categories, locations, and hashes. It does
not expose workbook values or row-level log content.

In equivalence mode, baseline equivalence and candidate repeatability must both
pass. In correction mode, baseline differences are classified as correction
impact, while candidate repeatability and all evidence-integrity checks remain
mandatory. A correction-mode pass does not mean that corrected artifacts are
equivalent to the historical baseline.

Sanitized results are in [`VALIDATION.md`](VALIDATION.md). Detailed reports,
hashes, manifests, logs, and generated artifacts remain ignored and local.

## Upstream TriNetX Construction

This repository does not reproduce:

- the TriNetX query or export;
- diagnosis and procedure code lists;
- laboratory extraction and calendar-day windowing;
- most upstream variable derivations; or
- assembly and validation of `full_db.dta`.

Reproducing those steps requires the upstream preprocessing repository,
appropriate TriNetX access, and historical source evidence that is not
available here.

## Scientific Alignment

The phenotype inventory distinguishes six adjudicated simulated rules from four
definitions whose source alignment remains unresolved. The adjudication applies
only to the documented simulated rules, not to unavailable source-study
exclusions, settings, or repeat-measurement criteria.

See [`SCIENTIFIC_ALIGNMENT.md`](SCIENTIFIC_ALIGNMENT.md) before changing any
case definition, model, time window, missingness rule, or figure logic.

## Data Safety

Do not add PHI, restricted TriNetX data, row-level derivatives, credentials,
machine-specific paths, private drafts, or publisher-formatted article files.
Do not attach restricted or generated artifacts to public releases.
