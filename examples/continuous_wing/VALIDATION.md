# Validation record: continuous Goland wing

Values from `run_continuous_wing_study` (ContourRoots 0.6.1,
MATLAB R2023b, macOS). They are observed agreements for the stated
parameters and tolerances, not general accuracy guarantees.

| Check | Result |
|---|---|
| 8 vacuum frequencies below 900 rad/s vs closed form | max relative error 3.2e-14 |
| $\Delta$ with 1, 2, 4, 8 uniform strips (4 test points, U = 130 m/s) | max relative difference 1.4e-14 |
| flutter: `fzero` on the spectral abscissa | U = 136.983977449 m/s, ω = 70.033012953 rad/s |
| flutter: local Newton on $\Delta(i\omega,U) = 0$ | same to 1e-9 |
| FE + Hankel loads, 8/16/32/64/128 elements | 137.308301 / 137.065051 / 137.004259 / 136.989048 / 136.985245 m/s (second-order convergence) |
| larger window [-80 50 1 700] at 0.95 and 1.05 U_f | 6 modes; 0 and 1 unstable upper-half-plane mode |
| right-half-plane counts, 10–180 m/s every 10 m/s | 0 up to 130 m/s, 2 from 140 m/s |
| tapered wing, 2/4/8/16/32 strips | 151.58 / 155.81 / 156.81 / 157.05 / 157.11 m/s |
| vacuum 1 kN half-sine pulse vs 150-mode series | 1.9e-6 mm (response 13.3 mm) |
| 1 kN tip step at 120 m/s: FFT vs adaptive quadrature | 7.0e-4 mm |
| same, FFT vs de Hoog (0.6.1) | 6.7e-3 mm, within the tolerance of about 9e-3 mm (0.6.0: 4.3e-2 mm, false plateau reported as converged) |

Regression tests: `tests/regression/test_continuous_wing.m` (10 tests, ~20 s).

The research prototype (in the author's `delay_systems` repository, with a
Portuguese report) additionally validated the section loads on the NASA
typical-section benchmark and the hybrid 6×6 port formulation of each strip.
The section benchmark is already part of Tutorial 9. The port formulation
is not needed to compute the poles and will be revisited with a generic
network assembly.
