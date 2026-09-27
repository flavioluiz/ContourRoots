# Continuous wing: flutter without discretizing the structure or the air

The Goland wing, a uniform cantilever that bends and twists, with exact
Theodorsen strip aerodynamics. Its poles, its stability and its flutter
speed are computed from one exact characteristic function
$\Delta(s,U) = \det K(s,U)$. There are no finite elements, no modal
truncation and no rational approximation of Theodorsen's function. Explained in
[Tutorial 11](../../docs/tutorials/11_continuous_wing.md).

| Script | What it does | Time |
|---|---|---|
| `wing_quickstart.m` | vacuum check, stability counts, flutter speed, tip step responses | ~30 s |
| `run_continuous_wing_study.m` | every figure and table of Tutorial 11, with independent checks | ~90 s |

Run `setup_contourroots` first. The model functions are in
`examples/models`:

| Function | Role |
|---|---|
| `wing_model` | Goland (`'goland'`), vacuum (`'dry'`) and synthetic `'tapered'` wings |
| `wing_propagator` | exact transfer matrix $e^{A\ell}$ of one strip |
| `wing_matrix`, `wing_delta` | boundary-value matrix $K(s,U)$ and its determinant |
| `wing_transfer` | tip transfer functions (force/moment/torque → deflection/slope/twist) |
| `theodorsen_loads` | exact strip air loads, shared with the section model of Tutorial 9 |
| `wing_fem`, `wing_fem_matrix`, `theodorsen_loads_hankel` | independent finite-element / Hankel-form reference (validation only) |

The study checks:
- vacuum frequencies against closed-form formulas;
- invariance to the number of strips;
- a finite-element model with independent air loads, which converges to the
  continuous result;
- right-half-plane counts in larger windows;
- time responses against an analytic modal series and against adaptive
  quadrature.

See [VALIDATION.md](VALIDATION.md) for the recorded values.
See [HYBRID_VALIDATION.md](HYBRID_VALIDATION.md) for the matrix/hybrid extension.

The development extension [Tutorial 12](../../docs/tutorials/12_hybrid_mimo.md)
adds `hybrid_quickstart.m` and `run_hybrid_comparison.m`: explicitly connect
6-by-6 hybrid ports, compare the implicit `cdyn` model, and simulate/reuse
all tip channels through the matrix response API. The full comparison also
checks independent scalar-channel superposition and writes
`output/hybrid_comparison/` (about 2 min). The optional all-port run
(`[1 2 3]`) takes much longer; the fast algebraic tests need no stored
results.

This example was first developed as a research prototype, then reviewed
and migrated into the toolbox. Tutorial 11 uses the public scalar functions
`croots`, `cstep` and `clsim`. Tutorial 12 adds native matrix responses;
matrix-valued spectra remain future work; see
the [MIMO proposal](https://github.com/flavioluiz/ContourRoots/blob/main/docs/development/mimo_implementation_plan.md).
