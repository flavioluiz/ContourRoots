# cmodes — characteristic values of an analytic matrix

New in version 0.8.0. A tutorial introduction is [Tutorial 13](../tutorials/13_matrix_modes.md).

## Syntax

```matlab
H = @(s) [s+1 .1*exp(-s); 0 s+2];
[lambda,info] = cmodes(H,[-3 0 -1 1],'AssumeAnalytic',true);
assert(info.complete && max(abs(lambda-[-1;-2])) < 1e-10)
M = cdyn(H,[1;0],[1 0],0,'Dimensions',[2 1 1]);
[internalModes,im] = cmodes(M,[-3 0 -1 1],'AssumeAnalytic',true);
assert(im.complete && numel(internalModes)==2)
```

## Meaning and input contract

Find singularities of **H**, not poles of `G=C*(H\B)+D`. In the example,
`G=1/(s+1)` but `cmodes` also finds the hidden mode at -2. No observability,
controllability, cancellation or transmission-zero analysis is performed.

`H` is a scalar-node handle returning a finite square numeric matrix, a
constant numeric matrix, or a `cdyn` model. For `cdyn`, only `Factors{1}`
and its domain guard are used; B, C and D are not evaluated. Raw symbolic
matrices should be converted to a handle or supplied through `cdyn`.
`cmimo` is rejected. Dimensions must stay fixed; a direct handle's dimension
is established at the first permitted reference node in the search rectangle.
Sparse input is accepted but this initial backend uses **dense** LU/SVD.

`AssumeAnalytic=true` is required. It asserts regularity (det H is not
identically zero) and analyticity on a neighborhood of the closed rectangle.
The boundary must be nonsingular. This does **not** support an arbitrary
meromorphic transfer or hybrid port matrix. Avoid entire branch cuts, not
just their endpoints. Numerical completeness is conditional on adequate
contour sampling; it is not formal certification.

## Options

| Option | Default | Meaning |
|---|---|---|
| `AssumeAnalytic` | false | Required assertion for H; false raises an actionable error |
| `Derivative` | `[]` | Analytic H'(s); enables independent outer trace-integral check |
| `DomainCheck` | `[]` | Scalar logical guard, combined with a cdyn guard; rejected nodes raise `MatrixDomain` |
| `Scaling` | `'fixed'` | Fixed diagonal row/column equilibration; `'none'` disables it |
| `RootTolerance` | `1e-8` | Half-width of the local isolation square, times `1+abs(lambda)`; not a guaranteed point error |
| `RankTolerance` | `1e-8` | Rank threshold relative to the fixed scaled reference norm |
| `SeedPoints` | `[]` | Additional permitted refinement starts; never replace counts |
| `ContourPoints` | 12 | Initial samples **per edge**, doubled during refinement |
| `ContourRefinements` | 7 | Maximum doublings; two successive resolved agreements required |
| `MaxDepth` / `MaxCells` | 16 / 2048 | Subdivision depth and visited-cell bounds |
| `MaxIterations` | 60 | Newton iterations per seed |
| `MaxEvaluations` | 100000 | Total H and supplied-derivative evaluations, including repeats |
| `MaxMemoryMB` | 256 | Conservative dense workspace and contour-array gates; see limits below |
| `Display` / `Plot` / `Warn` | false / false / true | Table, new mode figure, and unresolved warning if info was not requested |

Unknown options are errors, including the not-yet-supported moment method,
`Singularities` exclusions and transfer-rank options. Derivative stencils and
Newton trials stay in the requested rectangle and respect both guards.

Without `Derivative`, a four-direction finite difference is used **only for
refinement**. The trace check is omitted explicitly, not silently claimed to
have passed. With `Derivative`, trapezoidal integration of `trace(H\Hprime)`
must agree with the LU winding to `1e-4*max(1,abs(count))`. An incorrect
derivative can therefore leave the outer contour unresolved.

## Outputs

`lambda` is a column of isolated numerical locations or **cluster
representatives**, sorted by decreasing real part then increasing imaginary
part. Consult `multiplicity` before interpreting a representative as one root.

| Field | Meaning |
|---|---|
| `kind` | `'matrix_characteristic'`, never transfer poles |
| `count`, `countComplete` | Outer algebraic count and whether it was resolved; NaN when unknown |
| `localCounts` | Algebraic count in each local isolation square |
| `multiplicity` | Numerical semisimple multiplicity; NaN for unresolved cluster/defective structure |
| `nullity` | Numerical nullity of scaled H at each representative; distinct from local count |
| `locationComplete`, `multiplicityComplete`, `complete` | All locations and local multiplicities accounted for; unresolved clusters prevent completeness |
| `status`, `certified` | `'numerically_complete'` or `'unresolved'`; certification always false |
| `locationRadius` | Circumradius of the local square, not a rigorous error bound |
| `unresolvedBoxes`, `clusters` | Unresolved search cells and counted local clusters, respectively |
| `rightVectors`, `leftVectors` | Cell arrays of numerical nullspace bases in original coordinates, columns normalized to unit 2-norm |
| `singularValues` | Singular values of the **fixed scaled** H at each representative |
| `residuals`, `leftResiduals` | Scaled directional residuals divided by the fixed regular-reference Frobenius norm |
| `unscaledResiduals`, `unscaledLeftResiduals` | Absolute residual norms with original H and unit original-coordinate bases |
| `simpleModeSlope` | Magnitude of scaled w'*H'*v/referenceNorm for simple modes; NaN otherwise, not an invariant physical condition number |
| `scaling`, `residualNormalization` | Fixed diagonals, reference node/norm, and normalization explanation |
| `traceCheck` | Availability, integral, agreement and omission/check reason |
| `history` | Count attempts, cells, sample counts and minimum LU pivot ratios (not singular-value condition numbers) |
| `evaluations`, `factorizations`, `linearSolves`, `svds`, `cells` | Callback calls, explicit factorization/solve/SVD calls and cells; excludes work hidden in callbacks and conditioning routines such as `rcond` |
| `stopReason`, `warnings`, `assumptions`, `options` | Termination and provenance |

Repeated semisimple roots are classified using local counts, numerical
nullity, machine-level coincidence and a nonsingular projected derivative.
A defective double root can have local count 2 but nullity 1. It is returned
as a counted unresolved cluster with `multiplicity=NaN`, not silently declared
a simple mode. Very close roots can likewise remain unresolved; try tighter
`RootTolerance` and a sensitivity study of `RankTolerance`. No partial
multiplicities or Jordan chains are computed in this first implementation.

## Algorithm and limits

At one regular reference node, freeze diagonal L and R and use A=L*H*R.
Pivoted LU supplies the sum of diagonal logarithmic magnitudes and phases,
including permutation parity. No product of pivots and no explicit
determinant of H is formed. Count acceptance checks phase increments below
pi/3 and log-magnitude increments below 2, and reconciles parent/child counts.
Internal split lines are shifted if needed; the user's outer boundary is
never moved. Simple refinement uses a bordered Newton solve; multiple-count
cells also use the logarithmic derivative `trace(A\A')` before local checks.

Resource exhaustion returns partial unresolved diagnostics, not successful
empty output. Shape errors, nonfinite matrix values and user callback errors
remain errors. A regular reference that cannot be found also leaves the
search unresolved. Scaling cannot recover information already lost inside a
callback or repair an intrinsically ill-conditioned eigenproblem.

Memory options gate estimated dense workspaces and contour arrays, not a
hard process-memory cap: callbacks, vendor LU/SVD workspaces, retained
history/vectors and MATLAB overhead are outside that estimate. The global
evaluation/cell limits bound total work. Nodes are not persistently cached;
no model state survives a search. Block contour moments and a sparse
large-scale backend are later enhancements.

See [Tutorial 13](../tutorials/13_matrix_modes.md) for determinant and wing
comparisons, and [matrix models](matrix_models.md) for structured inputs.
