# Hypercapnic Respiratory Failure Case Definitions

[![DOI](https://img.shields.io/badge/DOI-10.1016%2Fj.chest.2025.08.002-blue)](https://doi.org/10.1016/j.chest.2025.08.002)
[![PubMed](https://img.shields.io/badge/PubMed-40885535-green)](https://pubmed.ncbi.nlm.nih.gov/40885535/)
[![PMC](https://img.shields.io/badge/PMC-PMC12739763-green)](https://pmc.ncbi.nlm.nih.gov/articles/PMC12739763/)
[![Release](https://img.shields.io/github/v/release/reblocke/Hypercap-Case-Definitions)](https://github.com/reblocke/Hypercap-Case-Definitions/releases/latest)
[![Public reproducibility](https://github.com/reblocke/Hypercap-Case-Definitions/actions/workflows/public-reproducibility.yml/badge.svg)](https://github.com/reblocke/Hypercap-Case-Definitions/actions/workflows/public-reproducibility.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

> Code and reproducibility materials for the *CHEST* article
> **“The Consistency of Hypercapnic Respiratory Failure Case Definitions in
> Electronic Health Record Data.”**

This repository contains the Stata analysis and a Python/Graphviz cohort
diagram notebook. It does not contain patient-level TriNetX data.

## Status and Releases

| Version | Purpose |
| --- | --- |
| [`v2.0.0`](https://github.com/reblocke/Hypercap-Case-Definitions/releases/tag/v2.0.0) | Current LLM-assisted reproducibility and scientific-alignment update. Scientific decisions were adjudicated by the repository maintainer and study author. |
| [`v1.0.0`](https://github.com/reblocke/Hypercap-Case-Definitions/releases/tag/v1.0.0) | Paper-associated historical implementation before the current scientific-alignment corrections. |

The article was [published online on August 28,
2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12739763/) and appears in the
January 2026 issue of *CHEST* (169[1]:230-243).

The paper-associated release preserves the validated historical implementation,
but it does not establish the exact publication input or reproduce the upstream
TriNetX build. The current release includes documented scientific corrections
and should not be treated as identical to the paper-associated implementation.

## Study at a Glance

The study asks whether published electronic-health-record case definitions for
hypercapnic respiratory failure identify the same patients.

The analysis applies study-specific operationalizations of 10 published
definitions to 2022 adult emergency-department and inpatient encounters. It
examines agreement, cohort characteristics, mortality, and diagnosis-code
performance against laboratory-based hypercapnia measures.

Some source-study criteria could not be simulated from the available data.
Resolved and unresolved differences between the article and the current code
are listed in the [scientific alignment
register](docs/SCIENTIFIC_ALIGNMENT.md).

## Links and Identifiers

| Resource | Link |
| --- | --- |
| Article DOI | [10.1016/j.chest.2025.08.002](https://doi.org/10.1016/j.chest.2025.08.002) |
| PubMed | [PMID 40885535](https://pubmed.ncbi.nlm.nih.gov/40885535/) |
| PubMed Central | [PMCID PMC12739763](https://pmc.ncbi.nlm.nih.gov/articles/PMC12739763/) |
| Repository | [reblocke/Hypercap-Case-Definitions](https://github.com/reblocke/Hypercap-Case-Definitions) |
| Releases | [Versioned code snapshots](https://github.com/reblocke/Hypercap-Case-Definitions/releases) |

## Data Availability and Scope

The Stata workflow starts with a restricted, analysis-ready TriNetX file named
`full_db.dta`. The default local location is
`data/private/full_db.dta`.

The repository does not reproduce:

- the TriNetX query or export;
- diagnosis and procedure code-list construction;
- laboratory extraction or calendar-day windowing; or
- assembly of the analysis-ready dataset.

Investigators must obtain the data under their own institutional TriNetX
agreement. Do not commit source data, row-level derivatives, local exports, or
generated analysis results.

The [reproducibility guide](docs/REPRODUCIBILITY.md) documents the approved
input contract, upstream provenance, supported execution path, and limitations.
Passing public checks does not reproduce the restricted Stata analysis or
validate the article’s numerical results.

## Quick Start

### Public, Data-Free Checks

Use Python 3.11 and install the locked dependencies:

```bash
python3 -m pip install --require-hashes -r requirements.txt
```

Run the public repository checks:

```bash
make check
make diagram-smoke
```

`make diagram-smoke` also requires the system Graphviz `dot` executable. The
notebook renders fixed published aggregate counts into ignored `outputs/`
paths; it does not derive those counts from patient-level data.

### Restricted Stata Analysis

With approved access to `full_db.dta`, create its adjacent local approval
record:

```bash
make input-approve \
  INPUT_ROOT="/approved/restricted/hypercapnia" \
  APPROVE_RESTRICTED_INPUT=YES
```

Then run the guarded analysis:

```bash
make stata-run \
  STATA_BIN="/path/to/stata" \
  INPUT_ROOT="/approved/restricted/hypercapnia"
```

`make stata-run` is the sole supported scientific execution interface. Required
community Stata packages are listed in
[`metadata/stata_dependencies.csv`](metadata/stata_dependencies.csv).

Controlled comparisons use `make stata-compare` with one preserved baseline and
two independent candidate runs. See the [restricted validation
protocol](docs/REPRODUCIBILITY.md#restricted-validation-protocol) before using
that command.

## Paper-to-Code Map

| Paper component | Repository source |
| --- | --- |
| Cohort preparation, case definitions, agreement analyses, models, tables, and figures | [`Hypercapnia Case Definitions.do`](Hypercapnia%20Case%20Definitions.do) |
| Implemented phenotype rules and source alignment | [`metadata/phenotype_definitions.csv`](metadata/phenotype_definitions.csv) |
| Input and derived-variable documentation | [`data_dictionary.md`](data_dictionary.md) and [`data_dictionary.csv`](data_dictionary.csv) |
| Cohort diagram | [`Case Definitions Consort.ipynb`](Case%20Definitions%20Consort.ipynb) |
| Expected generated artifacts | [`metadata/output_manifest.csv`](metadata/output_manifest.csv) |
| Reproducibility boundaries and execution details | [`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md) |
| Restricted validation evidence | [`docs/VALIDATION.md`](docs/VALIDATION.md) |
| Article-to-code scientific decisions | [`docs/SCIENTIFIC_ALIGNMENT.md`](docs/SCIENTIFIC_ALIGNMENT.md) |

## Citation

Please cite both the article and the repository release or commit used.

Locke BW, Richards WW, Gouripeddi R, Brown JP, Anderson-Bell D, Finkelstein J,
Sundar KM, Peltan ID, Brown SM. The Consistency of Hypercapnic Respiratory
Failure Case Definitions in Electronic Health Record Data. *CHEST*.
2026;169(1):230-243. doi:
[10.1016/j.chest.2025.08.002](https://doi.org/10.1016/j.chest.2025.08.002).

Structured citation metadata are available in [`CITATION.cff`](CITATION.cff).

## Funding and Disclosures

The article record is authoritative for authorship, contributions, funding, and
disclosures. Support listed there includes the American Thoracic Society ASPIRE
Fellowship and grant and awards from the National Institutes of Health.

## License

Repository-authored code and documentation are available under the
[MIT License](LICENSE). The license does not cover TriNetX data, patient-level
derivatives, third-party software, or publisher-hosted article content.

## Contact

Open a GitHub issue or pull request for repository-specific questions.
The repository maintainer is Brian W. Locke
([ORCID 0000-0002-3588-5238](https://orcid.org/0000-0002-3588-5238)).
