# Validation

ContourRoots is tested at four levels. Run `buildtool test` and
`buildtool docs` from the repository folder; see
[CONTRIBUTING.md](../CONTRIBUTING.md).

## Unit tests (`tests/unit`)

The public interface: `croots`, `cpoles`, `czeros`, `cpzmap` and `ndpair`
agree with the core solver; polynomial input agrees with `roots`;
cancellations remove the right locations; exploratory and incomplete
searches warn; invalid input gives clear errors.

`test_matrix_modes` additionally checks LU winding against determinant
roots, permutation signs, scales from 1e-200 to 1e200, trace integration,
semisimple versus defective roots, close unresolved clusters, hidden modes,
domain guards and exhausted budgets. `test_matrix_modes_wing` compares
fixed-size wing assemblies against closed-form dry frequencies and the
scalar determinant route, and checks the exact flutter mode shape obtained
from the null vector (clamped root, free tip, independence of the number of
pieces). The reproducible full comparison is recorded in
[MATRIX_MODES_VALIDATION.md](../examples/continuous_wing/MATRIX_MODES_VALIDATION.md).

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

**Kernel reuse (`ckernel`).** `test_kernel_reuse` checks that:
- prepared and ordinary `clsim` responses agree (to $2\times10^{-6}$) for
  both holds and all three inversion methods, and for nonrational
  (diffusion, delay feedback) and unstable models;
- reuse evaluates the model zero times: the test evaluator counts its
  calls and is disabled after preparation; changing a parameter it
  captures does not change the prepared object;
- a delayed `tf` with direct transmission matches the delay-free `lsim`
  response shifted by the delay to $5\times10^{-8}$ (ramp, sine and step
  inputs, ZOH and FOH), and superposition holds to $10^{-12}$;
- a large, fast input is reported unresolved rather than silently
  refined; other grids, holds and model options are rejected; a failed
  preparation raises an error; properties are read-only; a MAT-file round
  trip gives identical results.

An additional check against `dde23` for $\dot x = -x - 0.5x(t-1) + \sin 2t$
agreed to $2\times10^{-5}$ (the ODE solver's accuracy), and the sine sweep
of `examples/time_response/reuse_kernels.m` reproduces $|G(i\omega)|$ to
$10^{-4}$.

The study's maximum absolute errors are approximately 1.2e-8 (delay
feedback), 4.6e-8 (diffusion), 1.1e-7 (finite rod), and 2.5e-5 (beam relative
to the refined FEM reference). NASA pitch-step FFT/quadrature differences
are 3.0e-11 at 160 ft/s and 6.7e-11 at 180 ft/s. These are observed errors
for the specified test grids and parameters, not general accuracy promises.

## Continuous wing (`test_continuous_wing`)

The Goland wing of [Tutorial 11](tutorials/11_continuous_wing.md) is a
continuous bending–torsion beam with exact Theodorsen strip loads. The
tests and the study `examples/continuous_wing/run_continuous_wing_study.m`
check that:
- the vacuum frequencies match the closed-form cantilever formulas
  (relative error at most $4\times10^{-14}$);
- the static tip flexibility matches the beam formulas ($10^{-12}$);
- $\Delta$ is the same function for 1 to 8 uniform strips
  ($1.4\times10^{-14}$), and for $n$ and $2n$ multiple-shooting pieces;
- the Laplace strip loads equal Theodorsen's lift and moment formulas
  written with Hankel functions ($10^{-12}$);
- $\Delta$ satisfies the Cauchy–Riemann equations and conjugate symmetry;
- the right-half-plane counts are 0 at 130 m/s and 2 at 145 m/s;
- the flutter speed from `fzero` on the spectral abscissa matches a local
  Newton solution ($10^{-9}$);
- an independent finite-element model with Hankel-form loads converges to
  the continuous flutter speed at second order (8 to 128 elements, from
  $2.4\times10^{-3}$ to $9.3\times10^{-6}$);
- the vacuum step and pulse responses match an analytic modal series;
- the FFT and adaptive-quadrature tip responses agree within tolerance.

The values are recorded in
[examples/continuous_wing/VALIDATION.md](../examples/continuous_wing/VALIDATION.md).

## Matrix models and hybrid assembly

`test_matrix_models`, `test_matrix_response` and `test_wing_hybrid` cover
the explicit matrix API, no constructor probing, batched shapes, shared
factorizations, domain/budget errors, optional symbolic/LTI adapters,
singular-descriptor rejection, both holds and all inversion methods,
delays/direct terms, singleton and rectangular models, inactive-input
reactivation, cancellation-aware errors, serialized immutable kernels,
and optional MIMO MATLAB `lsim` agreement, including a 2×2 model with a
different transport delay in each channel (ZOH, $1.4\times10^{-9}$). Scalar
kernel compatibility is also exercised by the existing suite.

The hybrid tests compare an element with two joined halves, check interface
compatibility/equilibrium and static compliance, compare every tip channel
against implicit and structured assembly, and verify explicit failure of
the hybrid chart at a cantilever eigenfrequency of a piece (an artificial
pole: the assembled transfer is finite there). The full
`run_hybrid_comparison` study adds both temporal assemblies and independent
scalar-channel superposition; results are saved with tolerances and runtime
metadata. It requires no pre-existing research artifacts.

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

## Shared MIMO kernel preparation (0.9.0)

`test_shared_grid` compares shared and independent FFT/de Hoog preparation,
including per-output tolerances, distinct channel delays, direct terms,
FOH/ZOH reuse, batched/structured callbacks, failure budgets, and loading
actual 0.8.0 snapshots. See the [benchmark report](benchmarks/shared_grid_preparation.md)
for measured node/factorization counts and the opt-in default decision.
