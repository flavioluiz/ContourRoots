# Related software

ContourRoots is not the only way to find roots of analytic functions. This
page lists alternatives and suggests when to prefer them. No systematic
benchmark has been published yet; the comparisons below are qualitative.

| Software | Language | Approach | Prefer it when |
|---|---|---|---|
| [cxroots](https://github.com/rparini/cxroots) | Python | contour integrals and moments, subdivision | you work in Python |
| [GRPF](https://github.com/PioKow/GRPF) | MATLAB | phase analysis on an adaptive triangular mesh; finds zeros **and poles** of a single meromorphic function | you only have $G$ as a single function and need both poles and zeros |
| [Root finding by Cauchy integration](https://www.mathworks.com/matlabcentral/fileexchange/68750-root-finding-cauchy-integration-method) | MATLAB | contour integrals | a small, simple code is enough |
| [RootsAndPoles.jl](https://github.com/fgasdia/RootsAndPoles.jl) | Julia | GRPF in Julia | you work in Julia |
| QPmR ([Vyhlídal and Zítek](https://doi.org/10.1007/978-3-319-01695-5_22)) | MATLAB | mapping of quasi-polynomial zeros | quasi-polynomials with many delays |
| [TDS-CONTROL](https://twr.cs.kuleuven.be/research/software/delay-control/) | MATLAB | spectral discretization of delay equations | **matrix** delay systems, stabilization and $H_\infty$ design |
| DDE-BIFTOOL | MATLAB | continuation and bifurcation of delay equations | nonlinear delay equations and bifurcations |

What ContourRoots adds, in the context of these tools:

- a MATLAB-style interface (`croots`, `cpoles`, `czeros`, `cpzmap`) that
  accepts handles, coefficient vectors, symbolic expressions and `tf`
  models;
- explicit handling of pole-zero cancellations for $N/D$ pairs;
- diagnostics that report unresolved regions instead of dropping them;
- tools for scalar time-delay stability (critical delays, crossing
  directions, bound-based counts) and a documented study of the limits of
  Padé approximations.
