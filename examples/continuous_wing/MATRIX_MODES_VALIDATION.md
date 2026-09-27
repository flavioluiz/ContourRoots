# M3 matrix-mode validation

Recorded on 2026-09-27, MATLAB R2023b, working tree based on 0.7.0.
Run `run_matrix_modes_comparison` from this directory with the toolbox on
the path. No Padé fit, modal truncation or finite-element approximation is
used by either compared route. Both share the exact strip propagator;
agreement between them is **not** independent validation of the air loads.
The dry closed-form frequencies provide an independent spectral reference.

## Fixed-size wing comparisons

One uniform physical strip, split into 1/2/4 exact numerical pieces:

| Model | Pieces | Matrix size | Count | Relative determinant discrepancy | Relative dry closed-form error |
|---|---:|---:|---:|---:|---:|
| dry | 1 | 12 | 4 | 2.25e-13 | 2.25e-13 |
| dry | 2 | 18 | 4 | 2.82e-15 | 5.63e-16 |
| dry | 4 | 30 | 4 | 2.60e-15 | 4.22e-16 |
| Goland, 150 m/s | 1 | 12 | 3 | 5.88e-16 | — |
| Goland, 150 m/s | 2 | 18 | 3 | 9.21e-16 | — |
| Goland, 150 m/s | 4 | 30 | 3 | 2.50e-16 | — |

Regions: dry `[-2 2 1 320]`; aerodynamic `[-40 15 1 320]`, above the branch
cut (left edge moved from -30 in review: a mode lies near -29.6+57i). Errors are `abs(lambda-reference)/(1+abs(reference))`. Acceptance:
both searches complete, equal counts, discrepancy below 1e-8; changing
pieces and the dry independent reference also below 1e-8.

Maximum fixed-scaled directional residual: 5.6e-15. Worst subdivision
discrepancy: 2.3e-13 (dry). The minimum LU pivot ratio across all count
contours improved from 1.7e-12 to 3.9e-11 for the dry example. This mixes
outer and near-mode local contours and is **not** a matrix condition number.
The independent frequency comparison is the stronger evidence that splitting
helps this case. No general monotonic conditioning or speed claim is made.

## Demonstrable advantage over explicit determinants

For `H=c*[0 s+2; s+1 0]`, c=1e-200, 1, 1e200, `cmodes` returns -1 and -2
with complete local and outer counts. At s=.5, the raw determinant is
respectively 0, -3.75 and -Inf. `croots(det(H))` reports unresolved for the
two extreme scales; it is complete only for c=1. The small regular case
also agrees with the scalar route to the stricter 1e-10 test threshold.

## Degeneracy and failure tests

Unit fixtures separately exercise a semisimple double root (count 2,
nullity 2), a defective double root (count 2, nullity 1, unresolved local
multiplicity), a close distinct cluster, scalar sinh with seven roots,
pivot parity, strongly unequal diagonal scales, hidden cdyn modes, guards,
boundary singularities and finite resource limits. Supplied derivatives
enable a trace-integral cross-check; numerical derivatives are used only
for refinement when no analytic derivative is supplied.

Defective clusters, Jordan-chain extraction and block moments are not
claimed solved by M3. See `docs/api/cmodes.md` for the qualified contract.
Generated MAT/CSV/PNG files live under ignored `output/matrix_modes/`; the
figure copied into `docs/assets/` is tracked with the tutorial.

## Flutter mode shape (added in review)

`wing_mode_shape` propagates the root state of the `cmodes` null vector
with the exact strip propagators. At U = 136.984 m/s the mode is
70.0330i; the propagated shape has clamped-root displacements below 1e-12
and free-tip loads 5e-14 relative to the root loads. Tip ratio
b|alpha|/|w| = 1.034 and phase -62.2 deg, identical with 1 and 2 pieces.
Test: `test_matrix_modes_wing/testFlutterModeShape`.

## Test and documentation run

The full unit/regression/documentation run passed all 137 tests then present,
with zero failures or skipped tests. Three additional independent-reference
and conservation tests were subsequently added and passed as part of the
10-test M3 unit suite; the two wing regressions and all new documentation
snippets were rerun successfully after the final local-cell containment
check. This covers 140 distinct tests across these runs, not a claim that
all were executed in a single 140-test invocation. No scalar solver or
time-response backend was changed. The updated PDF compiled successfully;
all three pages of the new chapter were visually inspected.
