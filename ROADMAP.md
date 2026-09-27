# Roadmap

Planned work and known limitations, roughly in order of priority. Items
come from the review of version 0.1.0 (the critical findings were fixed in
0.2.0, see [CHANGELOG.md](CHANGELOG.md)) and from the original
[repository plan](https://github.com/flavioluiz/ContourRoots/blob/main/docs/development/repository_plan.md). Contributions are
welcome; see [CONTRIBUTING.md](CONTRIBUTING.md).

## Next (0.3)

- **Continuous integration.** GitHub Actions with
  `matlab-actions/setup-matlab` running `buildtool test docs`, with and
  without the optional toolboxes, and a small matrix of MATLAB releases
  (at least one release older than R2023b) to establish a tested minimum
  version.
- **Stability map from crossings.** A function (e.g. `delay_stability_map`)
  that returns $Z(T)$ on an interval from `critical_delays` and the roots of
  $D+N$ ($Z(0)$), with explicit handling of tangential crossings and
  persistent imaginary roots. This is exact and cheap for large delays,
  where the contour count of `unstable_root_count` becomes inconclusive.
- **Quantitative documentation tests.** The documentation tests only check
  that code blocks run. Add assertions for the values stated in the text
  (counts, statuses, critical delays) to the main tutorials and examples.
- **Conditioning of Padé models.** `pade_characteristic` returns the expanded
  polynomial, which is badly conditioned for orders above about 10. Add a
  balanced state-space form (and use it in `delay_roots` seeding and in the
  studies), or evaluate in the scaled variable $z = sT$.
- **Frequency scaling in `critical_delays`.** Coefficient normalization and
  relative tolerances handle common scalings; an explicit frequency scaling
  (balancing the polynomial in $\omega^2$) would make the method robust to
  crossing frequencies spread over many orders of magnitude.

## Later

- **Multiple roots on the imaginary axis and tangential crossings** in
  `critical_delays`: higher-order direction analysis instead of the label
  `"tangent/degenerate"`.
- **Neutral systems** ($\deg N = \deg D$): diagnostics for essential
  spectra and for strong stability; the bound of `rhp_root_bound` does not
  apply.
- **Several delays and matrix delay equations**: at least an interface to
  pass a scalar determinant safely, and pointers to spectral methods.
- **Figures.** Larger labels in multi-panel figures for print; one style
  for all studies.
- **Provenance of generated results.** Record code version, MATLAB version
  and parameters with every study output, and check them when the manual
  is built.
- **Benchmarks** against GRPF, cxroots, QPmR and TDS-CONTROL on matched
  problems (missed and spurious roots, multiplicities, errors against
  references, evaluations and time), published with failures and tuning.
- **Performance.** Reuse of contour samples between parent and child cells;
  optional parallel evaluation for expensive models.
- **Distribution.** MATLAB File Exchange entry linked to GitHub releases;
  Octave compatibility only after the supported subset is tested.

## Not planned

- A package namespace (`+contourroots`) before 1.0, to keep unqualified
  calls stable.
- Interval-arithmetic certification, a GUI, and general MIMO nonlinear
  eigenvalue problems are out of scope for now. A bounded MIMO extension
  (matrix characteristic values, transfer poles, transmission zeros) is
  under study: see the proposal in
  [docs/development/mimo_implementation_plan.md](docs/development/mimo_implementation_plan.md)
  and its [review](docs/development/mimo/review.md). Nothing is scheduled.
