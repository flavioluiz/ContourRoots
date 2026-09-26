# Validation

ContourRoots is tested at four levels. Run `buildtool test` and
`buildtool docs` from the repository folder; see
[CONTRIBUTING.md](../CONTRIBUTING.md).

## Unit tests (`tests/unit`)

The public interface: `croots`, `cpoles`, `czeros`, `cpzmap` and `ndpair`
agree with the core solver; polynomial input agrees with `roots`;
cancellations remove the right locations; exploratory and incomplete
searches warn; invalid input gives clear errors.

## Regression tests (`tests/regression`)

- **Analytic spectra:** a triple root; three known roots under a scaling of
  the function by $10^{\pm100}$; $e^s$ with no roots; roots on the boundary
  and exhausted budgets must *not* be reported as complete.
- **Distributed models:** heat (three boundary conditions), string and
  duct poles against closed forms; total and partial cancellations.
- **Symbolic and LTI input** (skipped when the toolbox is absent).
- **Time delays:** the closed form $T_c = 2\pi/(3\sqrt3)$ of P1; the
  absence of crossings for P2; crossing residuals of P3; Padé predictions;
  bound-based unstable-root counts, including the stability window of P3.
- **Coupled beam:** static compliance, the decoupled limit (hidden modes),
  the bare-beam limit, invariance under dimensional scaling, and
  convergence of an independent finite-element model.
- **Baseline:** results recorded before the code was reorganized into this
  toolbox (research code at commit `57af5bf`) must be reproduced within
  $10^{-9}$ relative tolerance.
- **Aeroelasticity:** NASA and dimensional two-DOF flutter benchmarks;
  Bessel/Hankel identity, independent air-load assembly, still-air added
  mass, contour refinement, modal residuals, right-half-plane stability
  counts below and above flutter, a real unstable root above divergence,
  flutter by `fzero` on the spectral abscissa, an independent harmonic
  flutter solver, and the Jones state-space check.
- **Search domain:** a function that raises an error outside the search
  rectangle is never evaluated there (Newton steps and finite differences
  stay inside the current cell).
  See [Tutorial 9](tutorials/09_aeroelasticity.md) for the validation scope.

## Zero-state time responses

`test_response_api` and `test_time_response` test shifted FFT, de Hoog and
direct Bromwich quadrature against analytical first-order, delayed, diffusive,
fractional, unstable and oscillatory responses. Integrated ZOH/FOH kernels
are checked against exact recurrences and optional MATLAB `lsim`; tests
also cover direct Diracs, unknown initial limits, scalar-only evaluators,
conjugacy, numerical scale, refinement exhaustion and causal input handling.
The full [time-response study](tutorials/10_time_response.md) adds independent
method-of-steps, rod modal and beam FEM references plus an aeroelastic
input/output-channel cross-check. Convergence remains conditional on the
declared inversion half-plane and finite numerical resolution.

The study's maximum absolute errors are approximately 1.2e-8 (delay
feedback), 4.6e-8 (diffusion), 1.1e-7 (finite rod), and 2.5e-5 (beam relative
to the refined FEM reference). NASA pitch-step FFT/quadrature differences
are 3.0e-11 at 160 ft/s and 6.7e-11 at 180 ft/s. These are observed errors
for the specified test grids and parameters, not general accuracy promises.

## Documentation tests (`tests/docs`)

Every `matlab` block of the README and of `docs/` is executed, in order,
per document (except installation commands); the quick-start scripts are
executed; and code shown in a tutorial as the content of an example file
must be identical to that file.

## Reproducible studies (`examples/`)

The studies behind the tutorials and the manual compare results with
independent references (closed forms, modal formulas, finite elements,
the direct critical-delay method, time-domain simulation with `dde23`).
They are run with `buildtool examples` and write to `output/`.

## Environment

This version was tested with MATLAB R2023b on macOS, with and without the
Symbolic Math and Control System toolboxes (optional tests are skipped
when they are missing). Other releases have not been tested yet.
