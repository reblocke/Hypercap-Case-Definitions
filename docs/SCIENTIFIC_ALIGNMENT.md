# Scientific Alignment Register

This register compares the final article with the current public implementation.
It records review needs; it does not authorize changes to the analysis.

Authoritative article source:
[PMC12739763](https://pmc.ncbi.nlm.nih.gov/articles/PMC12739763/).

| ID | Published method or definition | Current implementation | Status | Required owner |
| --- | --- | --- | --- | --- |
| SA-001 | Stata 18 was used for the published analysis. | The do-file declares `version 17.0`, which requests Stata 17 compatibility behavior. | UNRESOLVED | Methods owner |
| SA-002 | The PaCO2 spline analysis used logistic GEE clustered by patient with an independent correlation structure. | The do-file uses ordinary `logistic` and does not use `patient_id` for clustering. | UNRESOLVED | Methods owner |
| SA-003 | The published Wilson simulated criterion is PaCO2 at least 50 mmHg with pH 7.35-7.45. | `def8` uses PaCO2 at least 45 mmHg with pH 7.35-7.45. | UNRESOLVED | Scientific owner |
| SA-004 | The published Vonderbank simulated VBG criterion uses venous pH above 7.35. The published Calvo PaCO2 criterion is above 45 mmHg, and Chung uses pH at most 7.45. | The code uses VBG pH at least 7.35, Calvo through `def3` with PaCO2 at least 45, and Chung with pH below 7.45. | UNRESOLVED | Scientific owner |
| SA-005 | The sensitivity analysis excludes four definition pairs with intentionally conflicting target populations. | The code selects observed pairwise kappa values greater than zero rather than naming the four pairs. | UNRESOLVED | Methods owner |
| SA-006 | Testing-strategy analyses describe ABG, VBG, and both-test subgroups. | Outputs labeled ABG-only and VBG-only filter on any ABG or any VBG, so encounters with both tests can enter both groups. | UNRESOLVED | Scientific owner |
| SA-007 | The spline figures represent diagnosis-code probability over PaCO2. | The plotting filters use the original `paco2_rounded` field after `xblc` generates prediction-grid variable `pa`; behavior may depend on row order. | UNRESOLVED | Methods owner |
| SA-008 | The article discusses informed presence and incomplete emulation of some source criteria. | Every implemented definition ultimately converts missing defining evidence to zero, and several published exclusions/settings were intentionally not simulated. | UNRESOLVED | Scientific owner |
| SA-009 | The article reports two-month mortality. | `died_2mo` uses `months_death_or_cens <= 1`; the upstream time encoding is not documented here. | UNRESOLVED | Data owner |
| SA-010 | Bülbül and Meservey depend on laboratory and diagnosis-code constructs described in the article. | `def5` and `def6` copy upstream flags whose complete derivations and producer version are unavailable here. | UNRESOLVED | Data owner |
| SA-011 | The final article uses Figure 3 for the PaCO2 spline and e-Figure 5 for the regional spline. | Current output names refer to Figure 2 and Figure S3. | UNRESOLVED | Scientific owner |

## Rules for Resolving an Item

An item may leave `UNRESOLVED` only when:

1. the designated owner approves the intended rule or method;
2. the upstream data contract needed to implement it is available;
3. the change is made in a separate scientific-correction ticket; and
4. the revised analysis is compared with a controlled baseline using the same
   approved restricted input.

Do not update an expected result merely to make a regression check pass.
