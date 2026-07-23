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
file named `full_db.dta`. With approved access to that file and licensed Stata,
the canonical guarded command is:

```bash
make stata-run \
  STATA_BIN="/path/to/stata" \
  INPUT_ROOT="/approved/restricted/hypercapnia" \
  OUTPUT_ROOT="outputs/stata"
```

The runner:

- refuses a missing input, unsafe run identifier, or existing run directory;
- records the analysis commit, dirty-worktree state, input and code hashes,
  Stata environment, and dependency hashes in an ignored manifest;
- preflights required community Stata dependencies;
- validates all documented runtime inputs before analysis;
- treats Stata's explicit status and completion artifacts—not its shell return
  code alone—as authoritative; and
- writes `SUCCESS` only when all 40 expected legacy artifacts are present and
  nonempty.

The input and output arguments may point to other approved local directories.
Neither the restricted input, the local manifest, detailed comparison report,
nor row-level derivatives may be committed.

Direct two-argument execution remains available as a compatibility path:

```bash
stata-mp -b do "Hypercapnia Case Definitions.do" "data/private" "outputs/stata"
```

That path does not provide the runner's collision protection, complete
provenance manifest, or output-inventory gate.

## Restricted Validation Protocol

For a code change that is intended to preserve results:

1. Use clean, isolated checkouts of the named baseline and candidate commits.
2. Verify the restricted input hash before each run.
3. Run the baseline once using the legacy two-argument interface.
4. Run the candidate twice using different run identifiers.
5. Compare the three isolated outputs with:

   ```bash
   make stata-compare \
     BASELINE_RUN="outputs/validation/baseline" \
     CANDIDATE_RUN_1="outputs/validation/candidate-1" \
     CANDIDATE_RUN_2="outputs/validation/candidate-2" \
     INPUT_ROOT="/approved/restricted/hypercapnia"
   ```

The comparator checks environment and input identity, dependency hashes,
semantic workbook content, decoded PNG pixels, required nonempty Stata graph
files, copied-do hashes, and normalized logs. It reports only discrepancy
categories, locations, and hashes—not cell values or row-level content.
Passing baseline equivalence and candidate repeatability are separate required
conditions. A sanitized record of the completed HCD-000B validation is in
[`VALIDATION.md`](VALIDATION.md); detailed evidence remains ignored and local.

The public checks do not execute Stata and do not establish that article
estimates were reproduced. A full reproduction claim requires a controlled run
against the approved input, the required community-contributed Stata commands,
and a documented software environment.

## Level 3: Upstream TriNetX Construction

This repository does not reproduce:

- the TriNetX query or export;
- diagnosis and procedure code-list construction;
- laboratory extraction and calendar-day windowing;
- derivation of upstream flags; or
- assembly and validation of `full_db.dta`.

The upstream repository is identified in
`metadata/upstream_dependency.yml`, but the exact producer commit and schema
version for the article dataset remain `UNRESOLVED`. A checkout commit and
artifact location observed during local validation are recorded only as
context; they do not establish which commit produced the restricted file.

## Scientific Alignment Boundary

The ten definitions in `metadata/phenotype_definitions.csv` are an inventory of
the current Stata implementation, not clinical approval. Known differences
between the final article and the current code are recorded in
`docs/SCIENTIFIC_ALIGNMENT.md`. Public validation must preserve unresolved
items rather than infer or approve a scientific resolution.

## Data-Safety Boundary

Do not add PHI, restricted TriNetX data, row-level derived data, credentials,
machine-specific paths, private drafts, or publisher-formatted article files.
Generated notebooks, figures, tables, logs, and Stata files belong under ignored
`outputs/` paths.
