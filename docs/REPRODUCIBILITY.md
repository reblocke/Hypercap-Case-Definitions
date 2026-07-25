# Reproducibility

This repository is a public, article-specific analysis companion for:

> Locke BW, et al. The Consistency of Hypercapnic Respiratory Failure Case
> Definitions in Electronic Health Record Data. *CHEST*. 2026;169(1):230-243.
> doi: [10.1016/j.chest.2025.08.002](https://doi.org/10.1016/j.chest.2025.08.002).

The repository supports several different levels of reproducibility. They should
not be treated as equivalent.

## Level 1: Public, Data-Free Checks

A clean public clone can:

- validate repository metadata and the public data-safety boundary;
- validate the code-derived phenotype, dependency, input, and output metadata;
- confirm that the tracked notebook is unexecuted and contains no embedded
  outputs; and
- execute the CONSORT-style diagram notebook using fixed published aggregate
  counts.

Install the locked Python environment with Python 3.11:

```bash
python3 -m pip install --require-hashes -r requirements.txt
```

Run all data-free checks:

```bash
make check
```

Render the diagram:

```bash
make diagram-smoke
```

The rendered notebook and TIFF are written under ignored `outputs/` paths. The
smoke command reports the installed Graphviz version because Graphviz is a
system dependency and image bytes may differ across versions or platforms.

The notebook contains fixed aggregate counts from the published analysis. It
does not calculate those counts from patient-level data or Stata outputs.

## Level 2: Restricted Downstream Analysis

The Stata workflow starts with an analysis-ready, restricted encounter-level
file named `full_db.dta`. The approved restricted input is bound to upstream
producer commit `44f49748d415e92b7d50b50d86b8fdea29f6cb07` and the
repository-defined observed schema `hypercapnia-full-db-v1`. This
owner-approved assignment is based on historical evidence; the upstream build
did not preserve source-file hashes or a clean-worktree attestation and is not
reproduced here. The selected derivations of `hypercap_on_abg` and
`hypercap_resp_failure` were separately verified from the approved
producer-commit Git object and are rechecked by the guarded input contract.

With approved access to that file, create its adjacent local approval manifest
once:

```bash
make input-approve \
  INPUT_ROOT="/approved/restricted/hypercapnia" \
  APPROVE_RESTRICTED_INPUT=YES
```

This writes ignored `full_db.manifest.json` beside `full_db.dta`. Its exact
fields are `schema_version`, `logical_name`, `size_bytes`, `sha256`,
`approved_at_utc`, `upstream_repository`, `producer_commit`,
`input_schema_version`, and `data_dictionary_sha256`. It contains no local path,
row count, or patient-level value. Replacement requires the separate explicit
setting `REPLACE_APPROVED_INPUT=YES`.

With licensed Stata, the canonical guarded command is:

```bash
make stata-run \
  STATA_BIN="/path/to/stata" \
  INPUT_ROOT="/approved/restricted/hypercapnia" \
  OUTPUT_ROOT="outputs/stata"
```

The runner:

- refuses a missing or unapproved input, approval drift, unsafe run identifier,
  or existing run directory;
- records the analysis commit, dirty-worktree state, input and code hashes,
  approved producer/schema reference, Stata environment, and dependency hashes
  in an ignored manifest;
- preflights required community Stata dependencies;
- validates all documented runtime inputs before analysis;
- requires both a zero Stata process return code and explicit successful status
  and completion artifacts;
- revalidates the input bytes, adjacent approval, and tracked authority after
  analysis; and
- writes `SUCCESS` only when all 40 expected artifacts are present and
  nonempty.

The input and output arguments may point to other approved local directories.
Neither the restricted input, the local manifest, detailed comparison report,
nor row-level derivatives may be committed.

`make stata-run` is the sole supported scientific execution interface. The
internal legacy argument mode is retained only to interpret preserved
historical validation evidence and must not be used to create new scientific
runs. The runner and analysis root must resolve to the same checkout so the
recorded commit, dirty state, harness hashes, and executed files describe one
tree.

## Restricted Validation Protocol

For a code change that is intended to preserve results:

1. Use clean, isolated checkouts of the named baseline and candidate commits.
2. Validate the adjacent restricted-input approval before each run.
3. Reuse the preserved legacy baseline; do not create a new scientific
   baseline through direct do-file invocation.
4. Run the candidate twice using different run identifiers.
5. Compare the three isolated outputs with:

   ```bash
   make stata-compare \
     BASELINE_RUN="outputs/validation/baseline" \
     CANDIDATE_RUN_1="outputs/validation/candidate-1" \
     CANDIDATE_RUN_2="outputs/validation/candidate-2" \
     INPUT_ROOT="/approved/restricted/hypercapnia"
   ```

For an owner-approved scientific correction, use the same protocol with
`COMPARISON_MODE=correction`:

   ```bash
   make stata-compare \
     COMPARISON_MODE=correction \
     BASELINE_RUN="outputs/validation/baseline" \
     CANDIDATE_RUN_1="outputs/validation/candidate-1" \
     CANDIDATE_RUN_2="outputs/validation/candidate-2" \
     INPUT_ROOT="/approved/restricted/hypercapnia"
   ```

Before assigning equivalence or repeatability labels, the comparator validates
the current adjacent approval, requires a legacy baseline, requires two guarded
candidates whose approval references match that manifest, and requires matching
nonempty candidate commit identifiers plus distinct nonempty candidate run
identifiers. A preserved legacy baseline may lack the new approval reference;
if it records one, that reference must match. Every artifact root must resolve
inside its own run directory and remain distinct from the other artifact
roots. The comparator recomputes the artifact and completion-control inventory
from disk, rehashes the current input, and compares that hash with every run
manifest even when no optional expected-hash pin is supplied. It then checks
environment and input identity, dependency hashes, semantic workbook content,
decoded PNG pixels, required nonempty Stata graph files, copied-do hashes, and
normalized logs. It revalidates the input and approval immediately before
reporting and atomically replaces any prior comparison report, so a stale pass
cannot survive an incomplete comparison. It reports only discrepancy
categories, locations, and hashes—not cell values or row-level content.
In equivalence mode, baseline equivalence and candidate repeatability are
separate required conditions. In correction mode, historical differences are
reported as correction impact and do not themselves fail the report, but
candidate repeatability and all evidence-integrity gates remain required. The
three intentionally renamed Figure 3, e-Figure 5, and both-test heatmap files
are mapped explicitly to their historical names. A sanitized record of the
completed HCD-000B validation is in [`VALIDATION.md`](VALIDATION.md); detailed
evidence remains ignored and local.

The public checks do not execute Stata and do not establish that article
estimates were reproduced. A full reproduction claim requires a controlled run
against the approved input, the required community-contributed Stata commands,
and a documented software environment.

## Level 3: Upstream TriNetX Construction

This repository does not reproduce:

- the TriNetX query or export;
- diagnosis and procedure code-list construction;
- laboratory extraction and calendar-day windowing;
- derivation of upstream flags other than the two selected aggregates verified
  below; or
- assembly and validation of `full_db.dta`.

The upstream repository and the owner-approved historical producer/schema
assignment are recorded in `metadata/upstream_dependency.yml`. This assignment
does not independently reproduce upstream construction. The historical build
did not record source-file hashes or clean-worktree state. The
`hypercap_on_abg` and `hypercap_resp_failure` derivations were verified
separately from the producer-commit Git object; complete derivations for other
variables remain unavailable for review in this repository.

## Scientific Alignment Boundary

The ten definitions in `metadata/phenotype_definitions.csv` distinguish six
owner-approved simulated rules from four definitions that remain unapproved.
Approval applies only to the documented simulated rules, not to unavailable
source-study exclusions, settings, or repeat-measurement criteria. Resolved and
unresolved decisions are recorded in `docs/SCIENTIFIC_ALIGNMENT.md`. Public
validation must preserve unresolved items rather than infer an additional
scientific resolution.

## Data-Safety Boundary

Do not add PHI, restricted TriNetX data, row-level derived data, credentials,
machine-specific paths, private drafts, or publisher-formatted article files.
Generated notebooks, figures, tables, logs, and Stata files belong under ignored
`outputs/` paths.
