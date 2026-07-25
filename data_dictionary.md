# Data Dictionary

This dictionary documents the restricted input variables expected by
`Hypercapnia Case Definitions.do` and the main variables derived by the current
Stata implementation. It is documentation only; no TriNetX row-level data are
included in this repository.

The machine-readable companion file is `data_dictionary.csv`.

The CSV contains one row per variable: 56 runtime inputs, seven context-only
input fields that are not referenced by the current do-file, and 27 variables
derived in Stata. Input and output artifacts are documented separately in
`metadata/upstream_dependency.yml` and `metadata/output_manifest.csv`.

`verified` means that the documented rule was checked against the named source
and approved for this repository. `draft` means that a description or
derivation was observed from the downstream code but has not been verified
against an authoritative source. `needs_review` marks a known question requiring
human confirmation. `blocked` is reserved for questions that cannot be resolved
without missing upstream documentation.

## Restricted Input

| File | Unit of observation | Required status | Notes |
| --- | --- | --- | --- |
| `data/private/full_db.dta` | One adult emergency-department or inpatient encounter | Required, local only | TriNetX-derived encounter-level analytic dataset bound to the approved producer/schema record and an adjacent local approval manifest; see `metadata/upstream_dependency.yml`. |

The approved restricted input is bound to upstream producer commit
`44f49748d415e92b7d50b50d86b8fdea29f6cb07` and the repository-defined
observed schema `hypercapnia-full-db-v1`. This owner-approved assignment is
based on historical evidence; the upstream build did not preserve source-file
hashes or a clean-worktree attestation and is not reproduced here. The
derivations of `hypercap_on_abg` and `hypercap_resp_failure` were separately
verified from the producer-commit Git object and are rechecked by the guarded
input contract. This does not verify any other source-variable derivation marked
`blocked` or `needs_review` in `data_dictionary.csv`.

## Source Variable Groups

| Group | Variables | Notes |
| --- | --- | --- |
| Identifiers and timing | `patient_id`, `encounter_id`, `admission_date`, `first_encounter`, `encounter_type` | Encounter-level identifiers and timing; identifiers remain restricted. Only `first_encounter` is referenced by the current do-file. |
| Encounter setting | `is_emer`, `is_inp`, `location` | Analysis keeps emergency or inpatient encounters and uses `location` for regional sensitivity analyses. |
| Demographics | `age_at_encounter`, `female`, `sex`, `black_race`, `race`, `hisp_eth`, `ethnicity`, `bmi` | Public documentation avoids row-level values. `sex`, `race`, and `ethnicity` are context-only fields; the analysis uses the corresponding derived flags. |
| First-day blood gas and chemistry | `paco2`, `paco2_flag`, `highest_paco2_flag`, `vbg_co2`, `vbg_ph`, `vbg_po2`, `highest_vbg_co2_flag`, `vbg_or_abg_co2_flag`, `has_abg`, `has_vbg`, `has_vbg_and_cat`, `abg_ph`, `serum_hco3`, `acidemia`, `hypercap_on_abg` | First-calendar-day lab resolution follows TriNetX constraints. |
| Comorbidities | `chf`, `ckd`, `copd`, `nmd`, `osa`, `prim_met_alk`, `combo_met_alk` | Used for cohort description and reviewer-response summaries. |
| Diagnoses and procedures | `hypercap_resp_failure`, `ohs_code`, `has_j9602`, `has_j9612`, `has_j9622`, `has_j9692`, `has_j9600`, `has_j9601`, `has_j9610`, `has_j9611`, `has_j9620`, `has_j9621`, `has_j9690`, `has_j9691`, `resp_acid_dx`, `sleep_hypovent_dx`, `cchs_dx`, `other_sleep_hypovent_dx`, `other_abn_of_br`, `vent_proc`, `niv_proc`, `imv_proc`, `cc_time` | Binary flags derived upstream from diagnosis/procedure code lists. |
| Mortality and follow-up | `died`, `months_death_or_cens` | Used for survival analyses and short-term mortality indicators. |

## Emulated Case Definitions

| Variable | Label | Operational rule in this repository | Review status |
| --- | --- | --- | --- |
| `def1` | Adler | `paco2 >= 47.25` and `vent_proc == 1` | draft |
| `def2` | Thille | `paco2 >= 45`, `abg_ph < 7.35`, and `vent_proc == 1` | draft |
| `def3` | Ouanes-Besbes | `paco2 >= 45` and `abg_ph < 7.35` | draft |
| `def4` | Calvo | `paco2 > 45`, `abg_ph < 7.35`, and `niv_proc == 1` | verified |
| `def5` | Bülbül | nonmissing `paco2 >= 45`, asserted against `hypercap_on_abg` | verified |
| `def6` | Meservey | any of `ohs_code`, `has_j9602`, `has_j9612`, `has_j9622`, or `has_j9692`, asserted against `hypercap_resp_failure` | verified |
| `def7` | Vonderbank | `paco2 >= 45` or qualifying VBG criteria `vbg_co2 >= 45` and `vbg_ph > 7.35` | verified |
| `def8` | Wilson | `paco2 >= 50` and `7.35 <= abg_ph <= 7.45` | verified |
| `def9` | Cavalot | `paco2 >= 45` and `abg_ph <= 7.35`, or `vbg_co2 >= 50` and `vbg_ph <= 7.34` | draft |
| `def10` | Chung | `paco2 >= 45` and `abg_ph <= 7.45` | verified |

The complete code-versus-source inventory, missingness behavior, unavailable
criteria, and approval state are recorded in
`metadata/phenotype_definitions.csv`. The simulated rules for `def4`, `def5`,
`def6`, `def7`, `def8`, and `def10` are approved. This approval does not imply
that non-simulated source-study exclusions, settings, or repeat-measurement
criteria were implemented.

## Main Derived Variables

| Variable | Definition |
| --- | --- |
| `abg_vbg_confusion_matrix` | Testing-strategy category: neither ABG/VBG, ABG only, VBG only, or both on the first day. |
| `highest_any_flag` | Either `highest_paco2_flag` or `highest_vbg_co2_flag` is positive. |
| `died_1mo` | Death within the first month among encounters with death observed. |
| `died_2mo` | Death within two months among encounters with death observed. |
| `death_time` | Month of death for survival summaries. |
| `paco2_rounded` | PaCO2 rounded to 0.1 mmHg for spline models. |
| `rc*` | Restricted-cubic-spline basis terms generated by `mkspline2`. |
| `pa`, `odds`, `lb`, `ub` | `xblc` spline prediction grid and odds-scale intervals. |
| `prob_hypercap`, `pr_lb`, `pr_ub` | Probability-scale transformation of spline odds and intervals. |
| `log_odds`, `log_lb`, `log_ub` | Log-scale transformation of spline odds and intervals. |

## Related Metadata

- `metadata/phenotype_definitions.csv` records the ten current implementation
  rules and distinguishes approved simulated rules from unresolved definitions.
- `metadata/output_manifest.csv` records the 12 generated artifact families.
- `metadata/stata_dependencies.csv` records directly invoked
  community-contributed commands and the graphics scheme.
- `metadata/upstream_dependency.yml` records the restricted input boundary and
  the owner-approved historical producer/schema assignment and its limitations.
