# ContourRoots documentation

ContourRoots finds roots of scalar analytic functions and poles and zeros of
scalar nonrational transfer functions inside a rectangle of the complex
plane. This page is the map of the documentation.

## Start here

1. [Getting started](getting_started.md): install, find your first roots,
   plot them and read the diagnostics. (10 minutes)
2. [Tutorial 1 — First roots](tutorials/01_first_roots.md): polynomials,
   exponentials and a time-delay equation solved by hand and by ContourRoots.

## Tutorials

The tutorials build on each other. Each one has runnable code; paste the
blocks into MATLAB in order.

| # | Tutorial | You will learn |
|---|---|---|
| 1 | [First roots](tutorials/01_first_roots.md) | `croots`, the search rectangle, `AssumeAnalytic`, reading `info` |
| 2 | [Poles and zeros](tutorials/02_poles_and_zeros.md) | modes versus transfer poles, `ndpair`, cancellations, `cpzmap` |
| 3 | [Describing a model](tutorials/03_model_inputs.md) | handles, coefficient vectors, symbolic input, `tf`/`ss`, options |
| 4 | [How it works](tutorials/04_how_it_works.md) | argument principle, subdivision, Newton, cancellation checks |
| 5 | [Time-delay systems](tutorials/05_time_delay_systems.md) | critical delays, crossing directions, stability maps |
| 6 | [The pitfalls of Padé](tutorials/06_pade_pitfalls.md) | what a Padé approximation keeps and loses |
| 7 | [Distributed-parameter systems](tutorials/07_distributed_systems.md) | heat, string, duct and beam transfer functions |
| 8 | [Beam coupled to an oscillator](tutorials/08_coupled_beam.md) | deriving N/D for a PDE–ODE system, parameter studies, FEM check |
| 9 | [Flutter without rational approximations](tutorials/09_aeroelasticity.md) | stability of a wing section with exact Theodorsen aerodynamics: unstable-root counts, root locus in the airspeed, flutter and divergence; Jones, quasi-steady and p-k comparisons |

## Reference

- [Function reference](api/index.md): every public function, its inputs,
  options, outputs and every field of `info`.
- [Diagnostics and limits](diagnostics_and_limits.md): what "numerically
  complete" means, and when the method can fail.
- [Troubleshooting](troubleshooting.md): symptoms and fixes.
- [Validation](validation.md): how the code is tested.
- [Related software](related_software.md): other root finders and when to
  prefer them.
- [ContourRoots manual (PDF)](ContourRoots_manual.pdf): the mathematical and
  algorithmic background, with proofs and case studies, in article form.

## Conventions used in the documentation

- $s = x + iy$ is the complex variable; `region = [xmin xmax ymin ymax]`
  is the rectangle $x_{\min} < x < x_{\max}$, $y_{\min} < y < y_{\max}$.
- Code blocks assume that ContourRoots is on the MATLAB path (run
  `setup_contourroots` once per session, or install the Add-On).
- Polynomial coefficients are in descending powers, as in `roots` and `tf`.
- "Exact" in this documentation means that the original function is
  evaluated without Padé approximation or modal truncation. The computed
  roots are still floating-point numbers.
