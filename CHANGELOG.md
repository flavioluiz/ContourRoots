# Changelog

All notable changes to ContourRoots are documented here. The project
follows [semantic versioning](https://semver.org/).

## [0.1.0] — 2026-09-26

First release.

### Added
- `croots`, `cpoles`, `czeros`, `cpzmap` and `ndpair`: a MATLAB-style
  interface to the contour-based solver, accepting function handles,
  coefficient vectors, symbolic expressions and continuous SISO LTI models.
- Pole-zero cancellation for numerator/denominator pairs.
- Time-delay tools: `critical_delays`, `unstable_root_count`,
  `rhp_root_bound`, `delay_roots`, `delay_root_count` and the Padé
  utilities `pade_delay`, `pade_characteristic`, `pade_critical_delays`,
  `match_pade_roots`.
- Examples: quick starts, time-delay and Padé studies, distributed-parameter
  models (heat, string, duct, beam) and a beam coupled to an oscillator,
  with finite-element validation.
- Documentation: README, getting started, eight tutorials, function
  reference, diagnostics, troubleshooting, validation, related software,
  and a PDF manual.
- Tests: unit, regression (including a baseline of the original research
  code) and documentation tests; `buildfile.m` tasks.

### Notes
- The numerical core (`complex_spectrum`) is unchanged from the research
  code it was extracted from; `characteristic_roots` and `transfer_poles`
  remain available under their original names.
