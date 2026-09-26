# Validation record

Executed with MATLAB R2023b on 2026-09-26 with ContourRoots after version
0.2.0. The example needed one change in the core solver: Newton steps no
longer evaluate the function outside the search cell (before, a
right-half-plane search could step onto Theodorsen's branch cut and stop
with an error). This is a record of the checks performed, not a proof of
global spectral completeness.

## Reproduction

From the repository root:

```matlab
run('setup_contourroots.m');
r = runtests('tests','IncludeSubfolders',true);
assertSuccess(r);
addpath('examples/aeroelasticity');
results = run_aeroelastic_study;
```

All tests passed (57 unit and regression tests, including twelve
aeroelastic tests, plus the documentation tests). The complete study also passed
its contour, benchmark, independent harmonic-solver and parameter-sweep
assertions.

## Numerical results

| Case | Exact flutter speed | Angular frequency [rad/s] | Reduced frequency |
|---|---:|---:|---:|
| NASA | 173.262350488 ft/s | 75.4620174044 | 0.435536151923 |
| Kaiser--Quero | 212.172922162 m/s | 58.4376826161 | 0.275424790405 |

These round to the published speeds, 173.26 ft/s and 212.2 m/s. NASA's
tabulated reduced frequency is 0.4355. The independent harmonic Hankel
formulation and the Laplace-domain determinant agreed within 1e-6 in
speed and 1e-8 in reduced frequency.

The exact sweeps used 51 NASA speeds (0:5:250 ft/s) and 61 dimensional
speeds (0:5:300 m/s), with two upper-half-plane roots at every speed in
`[-180 100 0.1 230]`. At six speeds per case, an enlarged rectangle
`[-220 130 0.05 300]` and 64 contour points reproduced both roots within
1e-6. Maximum modal residuals, measured as the smallest/largest singular
value ratio of the unscaled dynamic matrix, were 6.51e-16 and 1.40e-13.

Stability counts in the right-half-plane window `[1e-3 60 -250 250]`
(`ContourRefinements` 10) gave 0 unstable roots below and 2 above flutter
at every speed where the count was conclusive. Inconclusive counts occurred
only with a mode within 0.1 of the imaginary axis: for the dimensional case
at 5-30 m/s (almost no aerodynamic damping); none for NASA. Above the static
divergence speed (353.553 ft/s and 394.686 m/s), a real unstable root was
found at 1.05 U_d (s = 1.8876 and 6.8724).

The stiffness study used eleven frequency ratios from 0.30 to 0.80 and
checked contour roots at 98% and 102% of every computed flutter candidate.
Counts are numerical and local to these windows; the branch cut is excluded.

Jones flutter-speed errors were -0.234% and -0.560%. Setting C=1 gave
100.01064226 ft/s for NASA. For the dimensional case it was already
unstable at 0.01 m/s, so no positive onset is reported. The p-k neutral
speed agreed with the exact speed, but off-neutral damping differed.

Generated CSV, MAT and PNG files are in `output/aeroelasticity/` (ignored
by Git). The tutorial figures are retained in `docs/assets/`; regenerate
them through `buildtool examples` and `tools/make_doc_figures.m`.
The tutorial lists the exact publication locations.
