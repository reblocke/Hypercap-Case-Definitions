# Scientific Alignment Register

This register compares the final article with the current public implementation
and records the scientific owner's adjudications. HCD-001 items are marked
resolved only after the controlled restricted-data validation described below.

Authoritative article source:
[PMC12739763](https://pmc.ncbi.nlm.nih.gov/articles/PMC12739763/).

| ID | Published method or definition | Current implementation | Status | Required owner |
| --- | --- | --- | --- | --- |
| SA-001 | The article states that Stata 18 was used. | The scientific owner confirmed that the historical analysis used Stata 17; the do-file retains `version 17.0`. | RESOLVED | Methods owner |
| SA-002 | The PaCO2 spline analysis used logistic GEE clustered by patient with an independent correlation structure. | The do-file uses ordinary `logistic` and does not use `patient_id` for clustering. | UNRESOLVED | Methods owner |
| SA-003 | The published Wilson simulated criterion is PaCO2 at least 50 mmHg with pH 7.35-7.45. | `def8` now uses PaCO2 at least 50 mmHg with pH 7.35-7.45. | RESOLVED | Scientific owner |
| SA-004 | The published Vonderbank simulated VBG criterion uses venous pH above 7.35. The published Calvo PaCO2 criterion is above 45 mmHg, and Chung uses pH at most 7.45. | The three simulated criteria now use those published strict or inclusive boundaries. | RESOLVED | Scientific owner |
| SA-005 | The sensitivity analysis excludes four definition pairs with intentionally conflicting target populations. | The code now names Wilson-Thille, Wilson-Ouanes-Besbes, Wilson-Calvo, and Wilson-Cavalot explicitly rather than selecting pairs by their observed kappa. | RESOLVED | Methods owner |
| SA-006 | Testing-strategy analyses describe ABG-only, VBG-only, and both-test subgroups. | The three analyses now filter on mutually exclusive testing categories 1, 2, and 3. Pairwise kappa remains missing when either definition is constant within a subgroup. | RESOLVED | Scientific owner |
| SA-007 | The spline figures represent diagnosis-code probability over PaCO2. | The plotting filters use the original `paco2_rounded` field after `xblc` generates prediction-grid variable `pa`. HCD-001 changed all seven spline PNGs despite identical model-output transcript segments, while two corrected runs produced identical pixels; this is direct evidence of the known row-selection sensitivity. | UNRESOLVED | Methods owner |
| SA-008 | The article discusses informed presence and incomplete emulation of some source criteria. | Every implemented definition ultimately converts missing defining evidence to zero, and several published exclusions/settings were intentionally not simulated. | UNRESOLVED | Scientific owner |
| SA-009 | The article reports two-month mortality. | `died_2mo` uses `months_death_or_cens <= 1`; the upstream time encoding is not documented here. | UNRESOLVED | Data owner |
| SA-010 | Bülbül and Meservey depend on laboratory and diagnosis-code constructs described in the article. | The selected upstream aggregates were verified at the approved producer commit; `def5` and `def6` are now generated transparently from their components and asserted against those aggregates. The Bülbül persistence criterion and published Meservey exclusions remain non-simulated under SA-008. | RESOLVED | Data owner |
| SA-011 | The final article uses Figure 3 for the PaCO2 spline and e-Figure 5 for the regional spline. | The generated filenames now use Figure 3 and e-Figure 5; the comparison workflow maps the historical filenames explicitly. | RESOLVED | Scientific owner |

## Approved HCD-001 Decisions

- Preserve the pre-correction implementation as annotated tag
  `hcd-000b-historical`.
- Retain Stata 17 compatibility behavior for SA-001.
- Apply the published operational boundaries for SA-003 and SA-004.
- Exclude exactly the four owner-approved definition pairs listed in SA-005.
- Use mutually exclusive testing-strategy subgroups for SA-006.
- Make the two verified aggregate-dependent definitions transparent and
  contract-checked for SA-010.
- Align generated figure filenames with the final article for SA-011.
- Defer SA-002, SA-007, SA-008, and SA-009 without changing their scientific
  behavior.

The exact corrected commit `9f8fecd2384edd3ae3fd4b0c8ecbc8536769bc19`
was executed twice from clean isolated checkouts against the same approved
restricted input. Both runs passed the guarded completion contract and were
artifact-for-artifact repeatable. The correction-impact comparison and the
remaining SA-007 finding are summarized in
[`VALIDATION.md`](VALIDATION.md).

## Rules for Resolving an Item

An approved item may become `RESOLVED` only when:

1. the designated owner approves the intended rule or method;
2. the upstream data contract needed to implement it is available;
3. the change is made in a separate scientific-correction release; and
4. the revised analysis is compared with a controlled baseline using the same
   approved restricted input, and two independent candidate runs agree.

Do not update an expected result merely to make a regression check pass.
