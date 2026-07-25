# Hypercapnic Respiratory Failure Case Definitions

[![DOI](https://img.shields.io/badge/DOI-10.1016%2Fj.chest.2025.08.002-blue)](https://doi.org/10.1016/j.chest.2025.08.002)
[![PubMed](https://img.shields.io/badge/PubMed-40885535-green)](https://pubmed.ncbi.nlm.nih.gov/40885535/)
[![PMC](https://img.shields.io/badge/PMC-PMC12739763-green)](https://pmc.ncbi.nlm.nih.gov/articles/PMC12739763/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Stata analysis code (and a small CONSORT-style diagram notebook) for the CHEST
article **"The Consistency of Hypercapnic Respiratory Failure Case Definitions
in Electronic Health Record Data."**

## Article Links

- Final article DOI: <https://doi.org/10.1016/j.chest.2025.08.002>
- PubMed: <https://pubmed.ncbi.nlm.nih.gov/40885535/>
- PubMed Central / NLM full-text record: <https://pmc.ncbi.nlm.nih.gov/articles/PMC12739763/>
- Journal: *CHEST*. 2026;169(1):230-243.

## Project Summary

This repository contains code supporting a study asking whether common
electronic-health-record case definitions for hypercapnic respiratory failure
identify the same patients. The analysis applies study-specific
operationalizations based on 10 published definitions to 2022 adult
emergency-department and inpatient encounters from the TriNetX Research
Network, then compares agreement, cohort characteristics, mortality, and
diagnosis-code performance against laboratory-based hypercapnia measures.
The implemented rules and their source-alignment status are documented in
[`docs/SCIENTIFIC_ALIGNMENT.md`](docs/SCIENTIFIC_ALIGNMENT.md).

This repository contains code and documentation only; it does not include
patient-level data. The restricted workflow begins with a locally supplied,
analysis-ready `full_db.dta`. It does not reproduce the upstream TriNetX query,
export, code-list construction, laboratory windowing, or analytic-dataset
assembly.

## Authors, Funding, and Disclosures

Article authors: Brian W. Locke, W. Wayne Richards, Ramkiran Gouripeddi, Jeanette P. Brown, Dustin Anderson-Bell, Joseph Finkelstein, Krishna M. Sundar, Ithan D. Peltan, and Samuel M. Brown.

Repository maintainer: Brian W. Locke, ORCID `0000-0002-3588-5238`, GitHub `@reblocke`.

Support for this research listed in the article includes the American Thoracic
Society ASPIRE Fellowship and grant, NIH Ruth L. Kirschstein National Research
Service Award `5T32HL105321`, the National Center for Advancing Translational
Sciences, and the National Institute of General Medical Sciences. Use the
article record for the authoritative funding and disclosure statement.

## Repository Contents

| Path | Purpose |
| --- | --- |
| `Hypercapnia Case Definitions.do` | Main Stata workflow for cohort filtering, case-definition operationalization, agreement analyses, descriptive tables, Cox models, diagnosis-code performance, and figures. |
| `scripts/input_manifest.py` | One-time restricted-input approval and validation for the adjacent local manifest. |
| `scripts/run_stata.sh`, `scripts/stata_run.py` | Guarded restricted-data runner with unique run folders, provenance capture, and artifact-completeness checks. |
| `scripts/compare_stata_runs.py` | Value-suppressing equivalence or correction-impact comparator for one legacy baseline and two candidate runs. |
| `stata/` | Dependency preflight, input-contract validation, and neutral Stata driver files. |
| `Case Definitions Consort.ipynb` | Python/Graphviz notebook for the CONSORT-style case-definition diagram. |
| `data_dictionary.md`, `data_dictionary.csv` | Human- and machine-readable documentation for expected input and derived variables. |
| `docs/REPRODUCIBILITY.md` | Public, restricted downstream, and upstream reproducibility boundaries. |
| `docs/VALIDATION.md` | Sanitized equivalence, correction-impact, and repeatability results for the guarded Stata workflow. |
| `docs/SCIENTIFIC_ALIGNMENT.md` | Resolved and unresolved final-article alignment decisions. |
| `metadata/` | Phenotype approval inventory, upstream dependency record, Stata dependency inventory, and generated-output manifest. |
| `CITATION.cff` | Structured citation metadata for the repository and the preferred CHEST article citation. |
| `llms.txt` | Concise machine-readable project index for search and retrieval. |
| `AGENTS.md` | Repository-specific instructions for automated coding tools. |

## Data Requirements

The Stata workflow requires a restricted TriNetX encounter-level file named
`full_db.dta`. Its default location is:

```text
data/private/full_db.dta
```

Set `INPUT_ROOT` in the supported `make` commands to use another approved
directory.

The dataset is one row per emergency-department or inpatient encounter and must contain first-calendar-day laboratory, diagnosis, procedure, demographic, comorbidity, location, and mortality variables described in the data dictionary. The analysis uses TriNetX calendar-day lab resolution; first-day windows follow that convention.

TriNetX data must be re-requested under an investigator's institutional TriNetX agreement. Do not commit source data, derived row-level data, local exports, or other patient-level files to this repository.

The approved input record identifies upstream producer commit
`44f49748d415e92b7d50b50d86b8fdea29f6cb07` and schema
`hypercapnia-full-db-v1`. This provenance is based on historical repository
evidence rather than a reproduced upstream build; source-file hashes and
clean-worktree status were not preserved. The derivations of `hypercap_on_abg`
and `hypercap_resp_failure` were separately verified from the producer commit;
other upstream derivations remain outside the verified scope. See
[`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md) for the complete boundary.

## Workflow

### Public, Data-Free Checks

Use Python 3.11. Install the fully locked environment:

```bash
python3 -m pip install --require-hashes -r requirements.txt
```

Validate the public repository without Stata or restricted data:

```bash
make check
```

Render the fixed-count CONSORT-style diagram into ignored `outputs/` paths:

```bash
make diagram-smoke
```

The diagram notebook uses fixed published aggregate counts. Rendering it does
not reproduce those counts from the restricted analysis data and does not
validate the article's numerical results. See
[`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md) for the complete boundary.
`make diagram-smoke` uses the active locked Python 3.11 environment and requires
the system Graphviz `dot` executable.

### Restricted Stata Analysis

Install the required community Stata packages before running the full workflow.
The runner checks direct and transitive dependencies before loading the
restricted input. The complete code-derived inventory is in
`metadata/stata_dependencies.csv`.

After placing the approved `full_db.dta` under `INPUT_ROOT`, create its local
approval manifest once:

```bash
make input-approve \
  INPUT_ROOT="/approved/restricted/hypercapnia" \
  APPROVE_RESTRICTED_INPUT=YES
```

This writes an ignored `full_db.manifest.json` beside the restricted data. The
manifest records the input hash and approved provenance and schema, but no local
path or row-level value. Replacing an existing manifest also requires
`REPLACE_APPROVED_INPUT=YES`.

Then run the canonical guarded workflow from the repository root:

```bash
make stata-run \
  STATA_BIN="/path/to/stata" \
  INPUT_ROOT="/approved/restricted/hypercapnia" \
  OUTPUT_ROOT="outputs/stata"
```

`INPUT_ROOT` is the directory containing `full_db.dta`. The runner creates a
unique ignored run directory, records provenance, verifies the approval and
input contract before and after analysis, and writes `SUCCESS` only after Stata
completes successfully and the required artifact inventory is complete.
Existing run directories are never reused.

`make stata-run` is the sole supported scientific execution interface. Run it
from the intended checkout.

To assess a change, reuse the preserved legacy baseline, run the guarded
candidate twice in clean isolated checkouts, then compare them:

```bash
make stata-compare \
  BASELINE_RUN="outputs/validation/baseline" \
  CANDIDATE_RUN_1="outputs/validation/candidate-1" \
  CANDIDATE_RUN_2="outputs/validation/candidate-2" \
  INPUT_ROOT="/approved/restricted/hypercapnia"
```

The default `COMPARISON_MODE=equivalence` requires the baseline and candidates
to match. For a documented, owner-approved scientific correction, add
`COMPARISON_MODE=correction` to the same command. Correction mode reports every
baseline-to-candidate difference as correction impact rather than treating the
difference itself as a failure; a passing report is not evidence that individual
differences were scientifically approved. The candidate runs must still match
each other, and all evidence-integrity checks must pass.

The comparator requires one legacy baseline and two distinct guarded candidates
from the same commit. It revalidates the approved input and run evidence, then
compares workbook contents, decoded image pixels, required Stata graphs,
normalized logs, dependencies, and candidate repeatability. The complete
protocol is described in
[`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md). Detailed reports remain
under ignored `outputs/`; sanitized validation outcomes are recorded in
[`docs/VALIDATION.md`](docs/VALIDATION.md).

The repository does not establish exact publication-time Stata package
versions. The code-derived dependency inventory marks them `UNRESOLVED` and
lists the packages required by the current analysis.

## Outputs

Generated outputs are written under ignored `outputs/` folders and should not
be committed as source files. The guarded Stata workflow produces unique
`outputs/stata/<run_id>/` folders containing a local provenance manifest, status
and dependency records, logs, copied do-files, tables, heatmaps, spline figures,
and temporary Stata graph files. The notebook writes the CONSORT diagram under
`outputs/figures/`.

The machine-readable inventory is `metadata/output_manifest.csv`.

Key paper-facing artifacts include:

- cohort characteristics table;
- case-definition relative-sensitivity, raw-agreement, kappa, and PABAK heatmaps;
- case-definition-by-workup and location summary tables;
- case-definition-specific descriptive summaries and survival analyses;
- diagnosis-code performance summaries against ABG and any-blood-gas reference standards;
- probability of hypercapnic respiratory failure code by day-1 PaCO2 overall, by encounter type, and by region.

## Citation

If using this repository, cite both the article and the specific repository commit or release.

**Article**

Locke BW, Richards WW, Gouripeddi R, Brown JP, Anderson-Bell D, Finkelstein J, Sundar KM, Peltan ID, Brown SM. The Consistency of Hypercapnic Respiratory Failure Case Definitions in Electronic Health Record Data. *CHEST*. 2026;169(1):230-243. doi:10.1016/j.chest.2025.08.002.

**Repository**

Locke BW, Richards WW, Gouripeddi R, Brown JP, Anderson-Bell D, Finkelstein J, Sundar KM, Peltan ID, Brown SM. Hypercapnic Respiratory Failure Case Definitions: repository materials. GitHub: <https://github.com/reblocke/Hypercap-Case-Definitions>.

Structured metadata are available in `CITATION.cff`.

## License

Code and repository-authored documentation are released under the MIT License. TriNetX data, patient-level derivatives, third-party software, and publisher-hosted article content are not covered by this repository license.

## Contact

Open a GitHub issue or pull request for repository-specific questions. For other correspondence, contact Brian W. Locke at `brian.locke@imail.org`.
