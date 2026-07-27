# Changelog

This changelog describes public, citable releases. Detailed scientific
decisions and validation evidence are maintained in
[`docs/SCIENTIFIC_ALIGNMENT.md`](docs/SCIENTIFIC_ALIGNMENT.md) and
[`docs/VALIDATION.md`](docs/VALIDATION.md).

## Unreleased

No changes yet.

## [v2.0.0] - 2026-07-27

### Current LLM-Assisted Update

- Aligned the implemented case-definition boundaries, sensitivity-analysis
  exclusions, testing-strategy subgroups, selected upstream aggregate checks,
  and figure filenames with decisions made by the repository maintainer and
  study author.
- Added a guarded restricted-data runner, input approval contract, provenance
  capture, complete artifact checks, and a value-suppressing comparison
  workflow.
- Added locked public checks, portable cohort-diagram rendering, structured
  metadata, and release-aware reader and machine documentation.
- Recorded the remaining scientific-alignment questions without changing their
  unresolved behavior.
- Used LLM tools to assist code review, reproducibility hardening, and
  documentation. Scientific decisions remained the responsibility of the
  repository maintainer and study author.

This release does not include restricted data or generated analysis results.
Its public checks do not reproduce the article estimates.

## [v1.0.0] - 2026-07-27

### Paper-Associated Historical Implementation

- Preserves the validated implementation before the current
  scientific-alignment corrections.
- Includes public reproducibility packaging around the historical scientific
  analysis.
- Does not establish the exact publication input, reproduce upstream TriNetX
  construction, or include restricted data or generated results.

[Unreleased]: https://github.com/reblocke/Hypercap-Case-Definitions/compare/v2.0.0...HEAD
[v2.0.0]: https://github.com/reblocke/Hypercap-Case-Definitions/releases/tag/v2.0.0
[v1.0.0]: https://github.com/reblocke/Hypercap-Case-Definitions/releases/tag/v1.0.0
