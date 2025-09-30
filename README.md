# Hypercapnic Respiratory Failure Case Definitions — Stata Analysis Code

**Manuscript:** *The Consistency of Hypercapnic Respiratory Failure Case Definitions in Electronic Health Record Data (CHEST)*. 
Online ahead of print August 28, 2025. doi:10.1016/j.chest.2025.08.002.

This repository contains Stata code to reproduce the analyses comparing commonly used case definitions for **hypercapnic respiratory failure** (HRF) in electronic health record (EHR) data.

---

## Quick links

- **Manuscript (CHEST):** doi:10.1016/j.chest.2025.08.002  
- **Contact (maintainer):** Brian W. Locke (brian.locke@imail.org)  
- **License:** MIT (see `LICENSE`)

---

## Summary (what this code reproduces)

- **Question.** Do different Hypercapnic Respiratory Failure case definitions used in EHR studies identify the same patients?  
- **Data.** 2022 adult ED and inpatient encounters from the **TriNetX Research Network**, which aggregates de‑identified EHR data from **76 U.S. healthcare organizations** (≈ **115M** patients at time of pull). Data extracted **2023‑06‑04**.  
- **Cohort.** 515,286 eligible encounters after applying completeness checks.  
- **Methods.** Emulated 10 published Hypercapnic Respiratory Failure case definitions on first calendar day of ED/hospital admission; compared agreement and cohort characteristics. Agreement metrics included **Cohen’s κ**, **relative sensitivity**, **positive percent agreement**; additional **PABAK**. Mortality via Cox models; diagnosis‑code accuracy vs ABG via logistic GEE with **restricted cubic splines**.  
- **Key findings.** Limited agreement across definitions (**median κ = 0.35**). Hypercapnic Respiratory Failure **diagnosis codes are insensitive** for first‑day ABG‑confirmed hypercapnia (**sensitivity 23.5%**), though PPV among coded patients with ABG is moderate. Cohorts differ materially in ventilatory support rates and mortality.

---

## What’s in this repository

- `code/`  
  - `hypercap_case_definitions.do` — main analysis script (cohort emulation, agreement metrics, figures/tables).
  - `Case Definitions Consort.ipynb` - jupyter notebook to create the consort diagram
- `output/` (created by the do‑files)  
  - Derived tables (case‑definition cross‑tabs, κ matrix, baseline characteristics).  
  - Figures (κ heatmap; PaCO₂ vs probability of HRF code; PRISMA‑style diagram).  
- `LICENSE` — MIT.  
- `README.md` — this file.

---

## Data requirements (replication‑ready)

Analyses assume one row **per encounter** with encounter‑day laboratory and administrative signals resolvable to **calendar day** of ED/hospital admission. TriNetX only resolves labs to calendar day; our first‑day windows use that convention. TriNetX data must be independently acquired as per their licensing agreements. 

### Minimum variables

**Identifiers and timing**
- `patient_id`, `encounter_id`
- `admission_date` (date), `encounter_type` (ED vs inpatient)

**Demography/comorbidity (used for cohort description)**
- `age`, `sex`, `race`, `ethnicity`, `BMI` (kg/m²)

**Laboratory (first calendar day of encounter)**
- **ABG PaCO₂ (mmHg)** — LOINC: **2019‑8**, **2026‑3**, **32771‑8**  
- **VBG pCO₂ (mmHg)** — LOINC: **115577‑6**, **2021‑4**  
- **Arterial pH**, **serum bicarbonate (mEq/L)** as available  
- **BMI** — LOINC: **39156‑5**  

**Diagnoses (ICD‑10‑CM)**
- HRF‑specific: **J96.02, J96.12, J96.22, J96.92**  
- Obesity hypoventilation: **E66.2**  

**Procedures (respiratory support)**
- ICD‑10‑PCS: **5A09459, 5A0945B, 5A09559, 5A0955B** (assistance with resp. ventilation: negative/continuous)  
- ICD‑10‑PCS: **5A09358, 5A09458, 5A09558** (intermittent CPAP: 24h, 24–96h, 96+h)  
- ICD‑10‑PCS: **5A0935Z, 5A0945Z, 5A0955Z** (assistance with resp. ventilation: 24h, 24–96h, 96+h)  
- ICD‑10‑PCS: **5A1935Z, 5A1945Z, 5A1955Z** (respiratory ventilation)  
- CPT: **101509, 1014859, 94002, 94003, 94660** (ventilation management)

**Administrative**
- Billing for critical care (if available)
- Death date or in‑network 60‑day mortality indicator (for outcomes)

---

## Case definitions emulated (high‑level)

We emulated **10** Hypercapnic Respiratory Failure case definitions from published EHR‑based studies. Types of criteria utilized: 
- **ABG‑based hypercapnia** at presentation (e.g., PaCO₂ ≥ 45 mmHg; with/without acidosis criteria).  
- **VBG‑based hypercapnia** at presentation.(e.g., PvCO2 ≥ 50 mmHg)
- **Diagnosis‑code** (ICD‑10‑CM).  
- **Procedure‑based** (ventilatory support codes), with or without laboratory corroboration.  

Exact simulated criteria and study sources are documented in the manuscript Tables/Figures and mirrored in code comments.

---

## Data availability and governance

- Data come from **TriNetX Research Network** and cannot be redistributed here. Investigators can **re‑request** equivalent datasets directly from TriNetX under their institutional agreements and reproduce the analysis using this code.  
- Because only **de‑identified** records were used, the study was determined **exempt** by the University of Utah IRB

---

## Support

This research was supported by: 
- the American Thoracic Society Academic Sleep Pulmonary Integrated Research/Clinical Fellowship (ASPIRE) Fellowship and grant 
- the National Institutes of Health under Ruth L. Kirschstein National Research Service Award 5T32HL105321 
- the National Center for Advancing Translational Sciences 
- the National Institute of General Medical Sciences 

---

## Citation

If you use this code or reproduce the analysis, cite the paper and this repository:

**Paper**  
Locke BW, Richards WW, Gouripeddi R, Brown JP, Anderson-Bell D, Finkelstein J, Sundar KM, Peltan ID, Brown SM. *The Consistency of Hypercapnic Respiratory Failure Case Definitions in Electronic Health Record Data.* **Chest**. 2025. doi:10.1016/j.chest.2025.08.002.

**Software (this repository)**  
Locke BW. **Hypercapnic Respiratory Failure Case Definitions — Stata Analysis Code.** MIT License. URL: GitHub repository (this page).

---

## License

MIT — see `LICENSE`. You may reuse or adapt the code with attribution.
