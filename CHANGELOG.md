# Changelog

All notable changes to ContourRoots are documented here. The project
follows [semantic versioning](https://semver.org/).

## [0.6.1] — 2026-09-27

### Fixed
- `'Method','dehoog'` could report convergence on a false plateau for
  responses with many lightly damped poles near the imaginary axis (the
  continuous-wing tip step at 120 m/s: 8.0987 mm instead of 8.0553 mm at
  t = 0.3 s, four times the tolerance). Successive degrees were only 1.4
  times apart and shared a truncated series that omitted higher modes. The
  degree ladder is now 10, 20, 40, 80, 160, 320: each check doubles the
  sampled bandwidth. Regression test `testDeHoogNoFalseConvergence`, which
  fails on 0.6.0.

### Documentation
- `SingularityBound` values supported by finite-rectangle root counts
  (Tutorials 10 and 11, quickstart, manual) are now stated as assumptions,
  not as guaranteed bounds.
- Tutorial 11 and manual chapter 12: correct explanation of condensation.
  A reduction without inversion, $\det E_{pp}(s)$, is equivalent to
  $\det K$; only condensations that invert $s$-dependent blocks create
  artificial poles. The full matrix is kept for conditioning.
- Links from packaged documents to `docs/development` (not shipped) now
  point to GitHub.
- ROADMAP lists the deferred items of stage M1 of the continuous-wing
  migration (hybrid-element demonstration, `ckernel` reuse, local MIMO
  superposition).

## [0.6.0] — 2026-09-27

### Added
- Tutorial 11 and `examples/continuous_wing`: poles, stability, flutter and
  time responses of the continuous Goland wing (bending–torsion cantilever
  with exact Theodorsen strip loads), with no structural discretization, no
  modal truncation and no rational aerodynamic fit. Migrated from a research
  prototype after review. Validated against closed-form beam frequencies, an
  independent finite-element model with Hankel-form loads (second-order
  convergence to the continuous flutter speed, 136.984 m/s) and an analytic
  modal time series. Manual chapter 12.
- Example models `wing_model`, `wing_propagator`, `wing_matrix`,
  `wing_delta`, `wing_transfer` (adaptive multiple shooting at high
  frequency), `wing_fem`, `wing_fem_matrix`, `theodorsen_loads_hankel`.
- `theodorsen_loads`: exact strip air loads, now shared by
  `aeroelastic_matrix` (no change in results).
- Regression tests `test_continuous_wing`.

### Documentation
- Tutorial 10 and ROADMAP: the de Hoog method can report convergence on a
  wrong plateau for responses with many lightly damped modes; prefer FFT and
  cross-check with quadrature.

## [0.5.0] — 2026-09-26

### Added
- `ckernel(G,t,...)` prepares the step (and, for FOH, ramp) kernels of a
  system once; `clsim(K,u,t)` then simulates any input on the same grid and
  hold without evaluating $G$ or inverting again. Direct terms and known
  delays are kept. `K` is a read-only snapshot (no live function handle, no
  hidden cache) and can be saved to a MAT file.
- Separate kernel tolerances (default `1e-8`/`1e-6`) and output tolerances;
  the propagated output error is checked for every input, and an input
  that needs more accuracy is reported unresolved, never silently refined.
- `info.kernelReused` and `info.kernelPreparationEvaluations` in `clsim`.
- `ckernel` reference page, Tutorial 10 Section 10.5.1 (a sine sweep that
  reproduces $|G(i\omega)|$), manual section, example `reuse_kernels.m`.
- Regression tests: agreement with fresh responses for all methods and
  holds, zero model evaluations on reuse, an exact delayed reference,
  superposition, option/grid rejection, read-only properties, MAT round trip.

### Documentation
- `clsim` explains why MATLAB's `lsim` with FOH and an `InputDelay` differs
  from `clsim` just after the delay when $u(0) \ne 0$.

## [0.4.0] — 2026-09-26

### Added
- Time responses computed directly from nonrational transfer functions,
  without Padé or modal approximation: `cstep`, `cimpulse`, `clsim` (the
  analogues of `step`, `impulse` and `lsim`) and `cinvlaplace`. Numerical
  Laplace inversion on a Bromwich line by shifted FFT, with independent
  de Hoog and adaptive-quadrature methods; separate period, bandwidth and
  line-shift convergence checks; exact ZOH/FOH input reconstruction through
  integrated step/ramp kernels; Dirac terms reported as metadata.
- Tutorial 10, reference pages for the four functions and their options,
  manual chapter 11, four example scripts and a validation study (delay
  feedback, diffusion, heated rod, coupled beam and aeroelastic pitch
  response against independent references).
- `THIRD_PARTY_NOTICES.md`: BSD notice for the adapted mpmath de Hoog
  recurrence, included in the release packages.

### Fixed
- A sample exactly at the transport delay of an LTI model uses the exact
  right limit instead of a rounding-error positive time.

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
