# Scientific Alignment Register

This register compares the final article with the current public
implementation. It records decisions made by the repository maintainer and
study author and identifies questions that remain unresolved.

Authoritative article source:
[PMC12739763](https://pmc.ncbi.nlm.nih.gov/articles/PMC12739763/).

## Release Context

- `v1.0.0` preserves the paper-associated historical implementation.
- `v2.0.0` contains the adjudicated corrections listed below.
- A resolved item has a documented decision, an implementable data contract,
  controlled restricted-data validation, and repeatable candidate runs.
- An unresolved item must remain unchanged until the appropriate scientific,
  methods, or data authority adjudicates it.

## Article-to-Code Decisions

| ID | Published method or definition | Current implementation | Status | Resolution authority |
| --- | --- | --- | --- | --- |
| SA-001 | The article states that Stata 18 was used. | The study author confirmed that the historical analysis targeted Stata 17; the do-file retains `version 17.0`. | RESOLVED | Study methods |
| SA-002 | The PaCO2 spline analysis used logistic GEE clustered by patient with an independent correlation structure. | The do-file uses ordinary `logistic` and does not use `patient_id` for clustering. | UNRESOLVED | Study methods |
| SA-003 | The published Wilson simulated criterion is PaCO2 at least 50 mmHg with pH 7.35-7.45. | `def8` uses PaCO2 at least 50 mmHg with pH 7.35-7.45. | RESOLVED | Scientific definition |
| SA-004 | The published Vonderbank simulated VBG criterion uses venous pH above 7.35. Calvo uses PaCO2 above 45 mmHg, and Chung uses pH at most 7.45. | The three simulated criteria use the published strict or inclusive boundaries. | RESOLVED | Scientific definition |
| SA-005 | The sensitivity analysis excludes four definition pairs with intentionally conflicting target populations. | The code explicitly excludes Wilson-Thille, Wilson-Ouanes-Besbes, Wilson-Calvo, and Wilson-Cavalot rather than selecting pairs by observed kappa. | RESOLVED | Study methods |
| SA-006 | Testing-strategy analyses describe ABG-only, VBG-only, and both-test subgroups. | The analyses use mutually exclusive categories 1, 2, and 3. Pairwise kappa remains missing when either definition is constant within a subgroup. | RESOLVED | Scientific definition |
| SA-007 | The spline figures represent diagnosis-code probability over PaCO2. | Plotting filters use the original `paco2_rounded` field after `xblc` creates prediction-grid variable `pa`. The corrected implementation changed all seven spline images despite unchanged fitted-model transcript segments; two corrected runs produced identical pixels. | UNRESOLVED | Study methods |
| SA-008 | The article discusses informed presence and incomplete emulation of some source criteria. | Every implemented definition ultimately converts missing defining evidence to zero, and several published exclusions or settings were not simulated. | UNRESOLVED | Scientific definition |
| SA-009 | The article reports two-month mortality. | `died_2mo` uses `months_death_or_cens <= 1`; the upstream time encoding is not documented here. | UNRESOLVED | Data provenance |
| SA-010 | Bülbül and Meservey depend on laboratory and diagnosis-code constructs described in the article. | The selected upstream aggregates were verified at the approved producer commit. `def5` and `def6` are generated transparently from their components and asserted against those aggregates. The Bülbül persistence criterion and published Meservey exclusions remain non-simulated under SA-008. | RESOLVED | Data provenance |
| SA-011 | The final article uses Figure 3 for the PaCO2 spline and e-Figure 5 for the regional spline. | Generated filenames use Figure 3 and e-Figure 5; the comparison workflow maps historical filenames explicitly. | RESOLVED | Scientific definition |

## Decisions Included in v2.0.0

- Retain Stata 17 compatibility behavior for SA-001.
- Apply the published operational boundaries for SA-003 and SA-004.
- Exclude exactly the four definition pairs listed in SA-005.
- Use mutually exclusive testing-strategy subgroups for SA-006.
- Make the two verified aggregate-dependent definitions transparent and
  contract-checked for SA-010.
- Align generated figure filenames with the final article for SA-011.
- Preserve SA-002, SA-007, SA-008, and SA-009 without changing their scientific
  behavior.

The exact corrected commit
`9f8fecd2384edd3ae3fd4b0c8ecbc8536769bc19` was executed twice from clean,
isolated checkouts against the same approved restricted input. Both runs passed
the guarded completion contract and were artifact-for-artifact repeatable.
Historical correction impact and the remaining SA-007 concern are summarized
in [`VALIDATION.md`](VALIDATION.md).

## Resolution Standard

Before changing an unresolved item:

1. obtain an explicit decision from the relevant scientific, methods, or data
   authority;
2. confirm that the required upstream data contract is available;
3. implement the change in a separate scientific-correction release; and
4. compare the revised analysis with a controlled baseline using the same
   approved restricted input and two independent candidate runs.

Do not update an expected result merely to make a regression check pass.
