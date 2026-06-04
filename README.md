# Hypercapnic Respiratory Failure Case Definitions

[![DOI](https://img.shields.io/badge/DOI-10.1016%2Fj.chest.2025.08.002-blue)](https://doi.org/10.1016/j.chest.2025.08.002)
[![PubMed](https://img.shields.io/badge/PubMed-40885535-green)](https://pubmed.ncbi.nlm.nih.gov/40885535/)
[![PMC](https://img.shields.io/badge/PMC-PMC12739763-green)](https://pmc.ncbi.nlm.nih.gov/articles/PMC12739763/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Stata analysis code and a small CONSORT-style diagram notebook for the CHEST article **"The Consistency of Hypercapnic Respiratory Failure Case Definitions in Electronic Health Record Data."**

## Article Links

- Final article DOI: <https://doi.org/10.1016/j.chest.2025.08.002>
- PubMed: <https://pubmed.ncbi.nlm.nih.gov/40885535/>
- PubMed Central / NLM full-text record: <https://pmc.ncbi.nlm.nih.gov/articles/PMC12739763/>
- Journal: *CHEST*. 2026;169(1):230-243.

## Project Summary

This repository supports a study asking whether common electronic-health-record case definitions for hypercapnic respiratory failure identify the same patients. The analysis emulates 10 published definitions in 2022 adult emergency-department and inpatient encounters from the TriNetX Research Network, then compares agreement, cohort characteristics, mortality, and diagnosis-code performance against laboratory-based hypercapnia measures.

The repository intentionally contains code and documentation only. TriNetX-derived patient-level data are restricted and cannot be redistributed.

## Authors, Funding, and Disclosures

Article authors: Brian W. Locke, W. Wayne Richards, Ramkiran Gouripeddi, Jeanette P. Brown, Dustin Anderson-Bell, Joseph Finkelstein, Krishna M. Sundar, Ithan D. Peltan, and Samuel M. Brown.

Repository maintainer: Brian W. Locke, ORCID `0000-0002-3588-5238`, GitHub `@reblocke`.

Support listed in the article includes the American Thoracic Society ASPIRE Fellowship and grant, NIH Ruth L. Kirschstein National Research Service Award `5T32HL105321`, the National Center for Advancing Translational Sciences, and the National Institute of General Medical Sciences. Use the article record for the authoritative funding and disclosure statement.

## Repository Contents

| Path | Purpose |
| --- | --- |
| `Hypercapnia Case Definitions.do` | Main Stata workflow for cohort filtering, case-definition emulation, agreement analyses, descriptive tables, Cox models, diagnosis-code performance, and figures. |
| `Case Definitions Consort.ipynb` | Python/Graphviz notebook for the CONSORT-style case-definition diagram. |
| `data_dictionary.md`, `data_dictionary.csv` | Human- and machine-readable documentation for expected inputs, derived variables, case-definition flags, and outputs. |
| `CITATION.cff` | Structured citation metadata for the repository and the preferred CHEST article citation. |
| `llms.txt` | Concise machine-readable project index for search, retrieval, and future coding agents. |
| `AGENTS.md` | Repository-specific working rules for future coding agents. |

## Data Requirements

The Stata workflow expects a restricted TriNetX encounter-level Stata dataset named:

```text
data/private/full_db.dta
```

The dataset is one row per emergency-department or inpatient encounter and must contain first-calendar-day laboratory, diagnosis, procedure, demographic, comorbidity, location, and mortality variables described in the data dictionary. The analysis uses TriNetX calendar-day lab resolution; first-day windows follow that convention.

TriNetX data must be re-requested under an investigator's institutional TriNetX agreement. Do not commit source data, derived row-level data, local exports, or other patient-level files to this repository.

## Workflow

Install the required community Stata packages before running the full workflow. Observed dependencies include `missings`, `table1_mc`, `heatplot`, `kappaetc`, `diagt`, `mkspline2`, `xblc`, `cleanplots`, and related graphics/table dependencies.

Canonical Stata run from the repository root:

```bash
stata-mp -b do "Hypercapnia Case Definitions.do" "data/private" "outputs/stata"
```

The first argument is the directory containing `full_db.dta`; the second argument is the output root. If arguments are omitted, the script defaults to `data/private` and `outputs/stata`.

Optional notebook workflow:

```bash
python -m pip install -r requirements.txt
jupyter nbconvert --execute "Case Definitions Consort.ipynb"
```

The notebook requires both the Python `graphviz` package and the system Graphviz `dot` executable.

## Outputs

Generated outputs are written under ignored `outputs/` folders and should not be committed as source files. The Stata workflow produces dated run folders containing logs, copied do-files, tables, heatmaps, spline figures, and temporary Stata graph files. The notebook writes the CONSORT diagram under `outputs/figures/`.

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
