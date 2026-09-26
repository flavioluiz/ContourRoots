# Changelog

All notable changes to ContourRoots are documented here. The project
follows [semantic versioning](https://semver.org/).

## [0.3.0] — 2026-09-26

### Added
- Aeroelasticity example: flutter of a two-degree-of-freedom wing section
  with exact Theodorsen aerodynamics (`examples/aeroelasticity`,
  `examples/models/aeroelastic_*.m`, `theodorsen_laplace.m`), validated
  against NASA/TP-2015-218765 and Kaiser and Quero (2022); comparisons with
  Jones, quasi-steady and p-k aerodynamics; Tutorial 9 and a new manual
  chapter. `flutter_quickstart.m` shows a stability analysis in a few lines:
  unstable-root count in the right half-plane, modes, root locus in the
  airspeed, and flutter speed with `fzero`.

### Fixed
- Newton steps and finite-difference stencils never evaluate the function
  outside the current search cell. Before, a search whose rectangle ended
  next to a branch cut (e.g. the right half-plane for Theodorsen's
  function) could step onto the cut and stop with an error.

## [0.2.0] — 2026-09-26

Correctness release, following an independent review of 0.1.0 (carried out
with Codex). Several
situations in which 0.1.0 could return an incomplete or wrong result
without warning are fixed. Results that were correct are unchanged: the
regression baseline recorded from the original research code still passes.

### Fixed
- **Non-vectorized function handles.** A handle that returns a scalar for
  vector input (e.g. `@(s) prod([s-1, s+1])`) was treated as a constant, so
  `croots` could return no roots with status `numerically_complete`. Such
  handles are now evaluated point by point.
- **`unstable_root_count` and `delay_root_count`.** The count used a fixed
  number of samples and replaced negative or unresolved winding numbers by
  zero; for large delays it could report 0 unstable roots where there were
  thousands. Sampling is now adaptive, with the acceptance rule of the core
  solver, and an unresolved count is returned as `NaN` with the warning
  `unstable_root_count:Inconclusive`. `delay_roots` no longer treats an
  inconclusive count as a match.
- **`critical_delays` scale dependence.** Coefficients are normalized and all
  tolerances are relative, so multiplying `N` and `D` by a constant or
  rescaling time no longer loses crossings or small frequencies. The
  `(omega, T)` refinement uses scaled variables, and `zeroRootForAllDelays`
  uses a relative test.
- **`pade_characteristic`** removed legitimate small leading coefficients and
  could return a polynomial of lower degree than the Padé model (degree 8
  instead of 13 for `T = 0.1`, `n = 12`). Only exact zeros are removed now.
- **`pade_critical_delays`** used a fixed sampling grid, so enlarging `Tmax`
  could remove crossings. Each admissible branch is now solved directly from
  the continuous, monotone Padé phase; the option `GridSize` is ignored.

### Added
- `critical_delays` detects degenerate equations (`|D(iw)| = |N(iw)|` for
  every `w`) and roots that stay on the imaginary axis for every delay
  (common factors of `N` and `D`), with warnings and the fields
  `info.degenerate` and `info.persistentFrequencies`.
- `info.cancellationComplete` distinguishes total from partial pole-zero
  cancellation; `cpzmap` marks only total cancellations as "Cancelled".
- A mode-specific warning for exploratory pole searches of a single
  handle, recommending `ndpair(N,D)`.
- Regression tests for every issue above (`tests/regression/test_review_fixes.m`).
- `ROADMAP.md`.

### Changed
- Documentation and manual: when `unstable_root_count` is conclusive and
  when `Z = 0` implies stability; finiteness of zeros in a region;
  analyticity of the distributed-model factors in the search region; Padé
  crossing frequencies stated as candidates; softer wording on completeness.
- `buildtool manual` (and therefore `release`) regenerates every study first.
- The release packages also contain `CONTRIBUTING.md` and `ROADMAP.md`.
- The documentation test runner no longer closes figures it did not create.

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
