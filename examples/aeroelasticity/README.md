# Aeroelasticity example: flutter with exact Theodorsen aerodynamics

A two-degree-of-freedom wing section (plunge and pitch) with Theodorsen's
unsteady aerodynamics, evaluated exactly (Bessel functions, no rational
approximation). Explained step by step in
[Tutorial 9](../../docs/tutorials/09_aeroelasticity.md) and in the manual.

| File | Content |
|---|---|
| `flutter_quickstart.m` | stability at a given speed, modes, root locus in $U$ and flutter speed, in a few lines (about 15 s) |
| `run_aeroelastic_study.m` | complete study: two published benchmarks, stability counts, divergence, Jones/quasi-steady/p-k comparisons, parameter study, figures and CSV files in `output/aeroelasticity/` (about 2 minutes) |
| `VALIDATION.md` | record of the numerical results and checks of the study |
| `../models/aeroelastic_*.m`, `../models/theodorsen_laplace.m` | the model: section data, Theodorsen's function, dynamic stiffness matrix, characteristic function, flutter solvers and the rational comparison models |

```matlab
repo = fileparts(fileparts(which('croots')));
run(fullfile(repo, 'examples', 'aeroelasticity', 'flutter_quickstart.m'))
```

Benchmarks: B. Perry III, NASA/TP-2015-218765, Appendix C (standard case of
NACA Report 496); C. Kaiser and D. Quero, Aerospace 9(3):127, 2022, Eq. (33)
and Table 1. The source documents are not redistributed; all figures are
computed by this example.
