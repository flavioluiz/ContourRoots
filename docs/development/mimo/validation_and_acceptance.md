# MIMO validation and acceptance

Status: proposed. Parent document: [implementation plan](../mimo_implementation_plan.md).
Methods: [numerical methods](numerical_methods.md).

Fixtures below are specifications for future tests, not results of tests
already run. No MIMO implementation is delivered by these planning files.

## 1. Baseline and independence

Before implementation, record the then-current commit, MATLAB version,
available optional toolboxes, dirty files and full baseline results. Preserve
the current kernel-reuse work and documentation conventions. Rerun rather
than relying on a historical test count from a previous release.

Use three reference levels:

1. Hand-derived small rational/entire examples, including local orders.
2. Independent polynomial linearizations, minimal rational `pole`/`tzero`
   comparisons and stored high-precision transcendental reference values.
3. Existing physical models and independent scalar formulations, with
   grid/contour refinement and agreement across representations.

Calling the same scalar backend twice is a compatibility test, not an
independent mathematical oracle. Optional Symbolic/Control toolboxes may
generate or cross-check fixtures, but core runtime/tests must not require them.
Store reference provenance, parameter values, precision and generation method.

## 2. Mathematical fixture matrix

All reported counts below include multiplicities unless stated otherwise.
Choose rectangles with reference points strictly inside, not on boundaries.

| ID | Model | Required result / failure mode caught |
|---|---|---|
| F01 | `G=eye(2)` and `G=[1 0;0 1;1 1]` | Normal rank 2; no finite poles or transmission zeros |
| F02 | Constant zero/rank-deficient matrices, plus individual zero channels | Structural zeros accepted; zero channel is `identically_zero`; all-zero model is `degenerate`, not an arbitrary set of roots |
| F03 | `G=ones(2)/(s+1)` versus `G=eye(2)/(s+1)` | One pole location `-1` in each; multiplicity 1 versus 2; normal rank 1 versus 2; no finite zeros |
| F04 | `G=diag(1/(s+1)^2,1/(s+1))` | At `-1`: pole order 2, pole multiplicity 3; residue rank alone is insufficient |
| F05 | `G=[1 1;1 s+1]/(s+2)` | Pole `-2` multiplicity 2; transmission zero `0` multiplicity 1; channel `(2,2)` zero `-1` is not a transmission zero |
| F06 | `G=diag(1/(s+1),(s+1)/(s+2))` | At `-1` exponents `[-1,+1]`: pole and zero both retained; additional pole `-2`; determinant hides the coincidence |
| F07 | `G=[s;2*s]/(s+1)` and its transpose | Tall/wide: pole `-1` and zero `0`, both multiplicity 1; generic nullspace does not generate extra zeros |
| F08 | `G=[1;s]` and its transpose | No finite poles/zeros; selected minor `s` may produce a candidate at 0 that full-matrix validation must reject |
| F09 | `G=((s-2)/(s+1))*[1;2;3]*[1 2]` | Normal rank 1; pole `-1`, zero `2`, both multiplicity 1; persistent rank deficiency is not a continuum of transmission zeros |
| F10 | `H=diag(s+1,s+2), B=[1;0], C=[1 0], D=0` | `cmodes`: `-1,-2`; transfer pole: `-1`; no transmission zero; system-matrix zero `-2` is hidden-state structure |
| F11 | Quadratic `H=s^2*M+s*C+K`; also `(s+1)*eye(2)` and `[s+1 1;0 s+1]` | Polynomial linearization reference; last two have algebraic count 2 but geometric multiplicity 2 versus 1 |
| F12 | `H=diag(sin(s),1)` in `[-0.4,10] + i[-1,1]`; also `diag(exp(s),1)` | Four values `0,pi,2*pi,3*pi`, more than dimension 2; exponential fixture has none |
| F13 | `H=diag(s+1+0.5*exp(-s),s+2)` | Delay roots `-1+W_k(-0.5*exp(1))` plus `-2`; use independently generated high-precision branch references in a fixed rectangle, stating the branch set: in `[-6 1 -25 25]` MATLAB's `lambertw` branches `k=-4..3` give 8 delay roots (`k=0,-1` are the rightmost conjugate pair), so a symmetric `k=-3..3` reference misses one |
| F14 | Constant invertible left/right mixing of delay and distributed-parameter transfers | Same finite local pole/zero orders; changed channel zeros/directions; no rational approximation |
| F15 | Existing exact-Theodorsen analytic `H`; `cdyn(H,I,I,0)` | Within an analytic domain, modes equal transfer poles; `H^-1` has no finite transmission zeros there; never substitute the nonanalytic p-k model |
| F16 | Spectral points on rectangle edges; branch cut intersecting region; guarded evaluator | Unresolved/rejected as appropriate; no outer-boundary shift or out-of-domain derivative evaluation |
| F17 | Close distinct poles/zeros and nearly multiple characteristic values | Resolve with evidence or return an unresolved cluster; no proximity-only cancellation/multiplicity |
| F18 | Malformed dimensions, changing shape, nonfinite regular-node values and tiny budgets | Clear error or unresolved numerical result; never false complete/empty answers |
| F19 | Optional minimal/nonminimal `ss`, MIMO `tf/zpk` and delayed models | Faithful adapter semantics; hidden modes distinguished; no automatic Padé or discarded delay |
| F20 | Repeated runs, explicit axes, held plots, different probe settings | Reproducible results, unchanged global RNG/graphics state and correct MIMO labels |

Supplement F10 with `B=[1;1]`, `C=[1 delta]`: both poles exist for nonzero
`delta`, even when one is weakly visible. Sweep `delta` toward roundoff;
transition to uncertainty is acceptable, silent confident cancellation is not.

Supplement F11 with more than one eigenvalue per cell, insufficient probe
columns, repeated/defective values and exact count conservation after splits.
Do not let the moment extractor pass by returning only `n` values in F12.

Supplement F14 with an analytic unimodular mixing matrix such as
`[1 s;0 1]`, and with nonvanishing exponential factors. Finite local orders
remain invariant where transformations and their inverses are analytic.
Restrict distributed-parameter branches to declared domains. Overflow must
trigger diagnostics, not a switch to an approximate rational model.

## 3. Planned test files

Place files in the appropriate then-current runner-discovered directories;
verify test discovery rather than assuming a new folder is scanned.

| Proposed suite | Coverage |
|---|---|
| `tests/unit/test_matrix_models.m` | Constructors, sizes, numeric constants, structural zeros, symbolic/LTI optional adapters |
| `tests/unit/test_matrix_options.m` | Channel indices, incompatible options, budgets, rank assertions and domain contracts |
| `tests/regression/test_matrix_modes.m` | F11-F13; LU parity, trace count, moment capacity, left/right residuals |
| `tests/regression/test_mimo_poles.m` | F01-F04, F10, weak visibility, channel coverage and transfer multiplicities |
| `tests/regression/test_transmission_zeros.m` | F05-F09, full-rank checks, rectangular projections and normal rank |
| `tests/regression/test_mimo_local_structure.m` | F04/F06, defective roots, local partial multiplicities and high-order limits |
| `tests/regression/test_mimo_domains.m` | F16-F18, guards, branch cuts, derivative stencils, ambiguous neighborhoods |
| `tests/regression/test_mimo_compatibility.m` | Scalar calls, explicit 1-by-1 models, output shapes, warnings, RNG and optional adapters |
| `tests/graphics/test_mimo_cpzmap.m` | F20, axes/hold semantics, channel labels, simultaneous pole/zero markers |
| Existing documentation suite | New runnable examples and tutorial links after implementation; exclude proposal snippets |

Adapt the proposed graphics/unit paths to the repository's actual conventions
at implementation time. Add focused unit tests for LU permutation phase,
local coordinate scaling, Cauchy coefficient weights and determinantal orders.

## 4. Quantitative acceptance criteria

These are qualification targets for documented well-conditioned fixtures,
not guaranteed error bounds for arbitrary user functions.

| Quantity | Target |
|---|---|
| Location error | Matched-point `abs(lambda-lambda_ref)/(1+abs(lambda_ref)) <= 1e-7` |
| Coverage | Zero missed or spurious locations for the reference interior; exact expected algebraic totals |
| Matrix residual | Scaled right/left residual `<=1e-9`, using a published nonzero reference scale |
| Contour count | Integer discrepancy `<=1e-7`, two successive refinement agreements, and child/parent conservation |
| Independent paths | Agreement between direct/structured/channel forms where their mathematical problems coincide |
| Rank decisions | Stable classifications for threshold factors `0.1,1,10` on well-conditioned fixtures; otherwise explicit ambiguity |
| Local coefficients/orders | Stable under at least two contour resolutions and two admissible radii; all required order bounds accounted for |
| Degenerate fixtures | Correct status; no fabricated finite root list or known multiplicity where evidence is insufficient |
| Budgets | Stop within documented limits, preserving unresolved regions and candidates; no false completeness |
| Compatibility | Entire historical SISO and time-response/kernel suite passes unchanged |

Near ill-conditioned roots, a small residual need not imply a small location
error. Report isolation radii and conditioning. Compare invariant subspaces
rather than individual vectors when eigenvectors/directions are not unique.
Use permutation-aware one-to-one location matching; nearest-neighbor matches
alone can mask a missed root by matching two candidates to one reference.

Test global scalar factors `1e-100` and `1e100` where finite, plus fixed
diagonal channel scaling up to `1e-8`/`1e8`. Exact local orders are unchanged.
Well-conditioned equilibrated cases should pass; numerically indistinguishable
rank cases must become uncertain, not falsely validated. Also compare results
with scaling disabled and record the expected conditioning limitations.

Performance qualification should record dimension, contour nodes, evaluations,
factorizations, right-hand sides, peak planned storage and elapsed time on
one documented machine. Set regression ceilings only after baseline runs;
do not impose a speculative universal time limit. Hard evaluation/memory
budgets are testable independently of machine speed.

## 5. Adversarial and contract tests

- A handle that throws outside the region verifies all hidden probes and
  finite-difference stencils remain admissible. Constructor must not call it
  at `s=0` simply to infer dimensions.
- A known branch cut through the region is not fixed by listing only its
  endpoint. Reject the analytic assertion or require a valid subregion.
- `@(s) H(1i*imag(s))` is nonanalytic. Do not recommend it for contour searches;
  a user assertion cannot be automatically verified from finite samples.
- An assertion of deficient normal rank contradicted by a higher-rank witness
  is an error. Sampling alone without an upper bound cannot establish the
  deficient-rank completeness contract.
- A fixed projected determinant with a spurious zero must not produce a
  transmission zero after full-matrix validation (F08).
- Insufficient `ProbeSize`/`MomentOrder` must trigger adaptation or unresolved
  status, especially for F12. `MaxMinors`/`MaxLocalOrder` exhaustion cannot
  silently omit unexamined minors or higher local orders.
- Close pole/zero pairs of separation above/below requested tolerance must
  be resolved or clustered, never automatically cancelled as a pair.
- Repeated calls must not mutate user matrices, global RNG, graphics defaults
  or persistent model parameters. Preserve caller axes and hold state.
- Empty spectra must have evidence of resolved coverage; a caught evaluation
  failure cannot be converted to a successful empty answer.

## 6. Stage-specific delivery gates

| Gate | Required evidence before exposing the feature |
|---|---|
| P0 | Reviewed definitions/fixtures, current baseline, API collision check, preservation of unrelated changes |
| P1 | F01/F02/F18 adapters and channel tests; explicit 1-by-1 equivalence; legacy suite green |
| P2 | F11/F12/F13 analytic counts and extraction; boundary tests; count conservation; sufficient-capacity enforcement |
| P3 | F03 multiplicities, F10 hidden/weak modes, complete channel candidate coverage; F04 accurately flagged if unsupported |
| P4 | F05 square zeros and nonminimal F10 distinction; F06 solved or explicitly unresolved under a documented restricted scope |
| P5 | F07-F09 with rank provenance, artifact rejection and coverage; all rectangular limitations disclosed |
| P6 | F04/F06 fully resolved, partial-order oracle agreement, defective/local ambiguity tests, finite-work termination |
| P7 | Physical examples, optional adapters, docs/graphics, full regression report and reproducible benchmark metadata |

A restricted release after P4 may omit P5/P6 only with explicit unsupported
classes and incomplete statuses; it cannot advertise general MIMO spectra.
For the full planned feature, unresolved F04/F06 is a failed acceptance gate.
No gate passes because a plot looks plausible or Newton residuals are small.

## 7. Documentation and examples

Follow the current English reference style: Syntax, Description, Inputs,
Main options, Outputs, Examples, Errors and See also. Keep teaching narratives
in tutorials and derivations in the manual, linked to reference pages.

Planned deliverables once implementation is validated:

1. API pages for `cmimo`, `cdyn`, `cmodes`, `ctzeros`; targeted updates to
   `cpoles`, `czeros`, `cpzmap` and option/result references.
2. A tutorial using the next available number, without renumbering existing
   time-response/kernel documentation. Begin with the fully worked F05.
3. A second section/example showing hidden modes (F10), then the determinant
   counterexample (F06) once supported. Explain all three spectral objects.
4. `examples/mimo/` scripts for collective zeros, modal visibility, rectangular
   transfer matrices, delay systems and the validated aeroelastic model.
5. Manual derivation of matrix counting and local multiplicities; pseudocode,
   domain assumptions, runtime costs and known limitations.
6. README/API index/Contents/ROADMAP/validation updates consistent with actual
   supported stages. No feature claims before tests pass.

For the first tutorial, explicitly show

$$G(s)=\frac{1}{s+2}\begin{bmatrix}1&1\\1&s+1\end{bmatrix},\qquad
\det G(s)=\frac{s}{(s+2)^2}.$$

At `s=0`, direction `[1;-1]` is blocked even though all entries are nonzero.
At `s=-1`, the `(2,2)` channel is zero but the full matrix is invertible.
This demonstrates why channel zeros and transmission zeros need separate APIs.

Promote this proposed usage to runnable documentation only after the functions exist:

```matlab
G = cmimo({ndpair(1,[1 2]), ndpair(1,[1 2]); ...
           ndpair(1,[1 2]), ndpair([1 1],[1 2])});
region = [-3 1 -2 2];
[p,ip] = cpoles(G,region);                 % -2, multiplicity 2
[z,iz] = ctzeros(G,region);                %  0, multiplicity 1
[z22,i22] = czeros(G,region,'Output',2,'Input',2); % -1
assert(ip.complete && iz.complete && i22.complete);
```

Explain why these channel pairs supply analytic factor provenance, whereas
the mathematically equivalent opaque handle `@(s) [1 1;1 s+1]/(s+2)` alone
does not establish complete pole coverage. Explain the returned completeness
diagnostics; do not suppress warnings for presentation. MIMO time simulation
is explicitly a separate future topic.

## 8. Final implementation review checklist

- [ ] Public semantics and supported stages match code, help, Markdown and manual.
- [ ] Existing user/release edits were preserved; implementation diff is scoped.
- [ ] All tests are discovered and run; optional-dependency skips are explicit.
- [ ] Baseline SISO, simulation and kernel-reuse behavior is unchanged.
- [ ] Known references give correct locations, counts and partial multiplicities.
- [ ] Boundary, rank, domain and budget ambiguity cannot yield `complete=true`.
- [ ] No Padé, hidden modal truncation or unrequested numerical cancellation.
- [ ] Graphics preserve axes/hold state and distinguish matrix/channel meanings.
- [ ] Tutorial snippets execute; links and generated numerical values are checked.
- [ ] Regenerated PDF is visually reviewed if/when manual changes are authorized.
- [ ] Third-party licenses/references and optional dependencies are documented.
- [ ] Remaining restrictions are listed prominently before any release decision.

Creating these planning files does not authorize implementation, a commit,
version changes or publication; each belongs to its own requested workflow.
