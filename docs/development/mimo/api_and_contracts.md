# MIMO API and mathematical contracts

Status: proposed. Parent document: [implementation plan](../mimo_implementation_plan.md).

## 1. Proposed public functions

| Call | Meaning and output |
|---|---|
| `M = cmimo(channels)` | Transfer matrix from a rectangular cell array of scalar transfers/constants |
| `M = cmimo(G,'Size',[p m])` | Transfer evaluator: scalar complex `s` -> numeric `p`-by-`m` matrix |
| `M = cmimo(Gsym)` | Scalar-variable symbolic transfer matrix; optional Symbolic Toolbox |
| `M = cmimo(sys)` | Continuous MIMO `tf`/`ss`/`zpk`; optional Control Toolbox |
| `M = cdyn(H,B,C,D,...)` | Structured transfer `C*(H\B)+D` with original nonrational factors |
| `[lambda,info] = cmodes(H,region,...)` | Characteristic values of square analytic `H`; vectors in `info` |
| `[lambda,info] = cmodes(M,region,...)` | Same operation on `H` of a `cdyn` model, not its transfer poles |
| `[p,info] = cpoles(M,region,...)` | Poles of the input/output transfer matrix |
| `[p,info] = cpoles(M,region,'Output',i,'Input',j,...)` | Scalar poles of channel `(i,j)` |
| `[z,info] = czeros(M,region,'Output',i,'Input',j,...)` | Scalar zeros of channel `(i,j)` |
| `[z,info] = ctzeros(M,region,...)` | Finite transmission zeros of the matrix |
| `cpzmap(M,region,...)` | Transfer poles and transmission zeros, with explicit legend |
| `cpzmap(M,region,'Output',i,'Input',j,...)` | Existing SISO semantics for the selected channel |

`czeros(M,region)` without channel indices is an actionable error directing
the user to `ctzeros`. Require both channel indices, each scalar; do not
introduce cell-array output shapes through vector selectors in the first API.

Keep scalar `croots`, `ndpair`, `complex_spectrum`, `cpoles`, `czeros` and
`cpzmap` behavior unchanged. `cmodes` is distinct because singularity of an
internal dynamic matrix and poles of an observed transfer are different tasks.
Do not overload a raw numeric matrix in existing scalar functions.

An explicit 1-by-1 `cmimo` still has MIMO metadata, but its spectral locations
and orders must agree with the equivalent scalar operation. Ordinary SISO
inputs retain the original diagnostic schema and dispatch path.

## 2. Transfer model forms

### Channel representation

```matlab
% Proposed API: one scalar pair per nonconstant channel.
M = cmimo({ndpair(1,[1 1]), 0; ...
           ndpair(@(s) exp(-s),[1 2]), 2});
[p,info] = cpoles(M,[-4 1 -12 12],'AssumeAnalytic',true);
```

Cells accept existing scalar `ndpair`, handles, scalar symbolic expressions,
continuous SISO LTI objects and finite scalar constants. A numeric vector
inside a cell is ambiguous and must be wrapped in `ndpair`; zero is an
explicit structural zero. Do not relax `ndpair(0,D)` validation globally.

Store dimensions, representation type and per-factor analytic provenance.
Constructors validate declarations but must not evaluate an opaque model at
an arbitrary point such as `s=0`. First evaluation is a permitted point in
the search region. Constant matrices enter as `cmimo(numericMatrix)` and
are never interpreted as polynomial coefficient arrays.

### Matrix handle

`G(s)` receives a scalar complex number and returns a fixed-size matrix.
Use a dedicated matrix evaluator. Do not reuse the scalar adapter's
reciprocal, vectorization or shape-fallback semantics. Batched evaluation
is an optional later optimization with an explicit layout, not guessed
from a returned multidimensional array.

Reject changing dimensions, nonnumeric values or nonfinite values at
required regular evaluation nodes. A singular candidate is analyzed from
nearby admissible contours, not by replacing `Inf` or `0/0` with zero.

Opaque quotients have limited completeness provenance. A finite list of
sampled values cannot establish all their poles. Without analytic factors,
known analytic matrix structure or equivalent domain information, results
remain exploratory; failed function evaluation is an error, not a root.

### Structured representation

For `cdyn`, require `H`: `n`-by-`n`, `B`: `n`-by-`m`, `C`: `p`-by-`n`,
`D`: `p`-by-`m`. Each factor may be a numeric constant, symbolic matrix or
scalar-node matrix handle. Supply `'Dimensions',[n p m]` if handles prevent
dimension inference. All factors have one common spectral variable.

```matlab
% Proposed API: two coordinates, two forces and two measured displacements.
H = @(s) s^2*eye(2) + s*diag([0.2 0.3]) + [2 -1;-1 2] ...
         + 0.1*exp(-0.4*s)*eye(2);
M = cdyn(H,eye(2),eye(2),zeros(2),'Dimensions',[2 2 2]);
[modes,im] = cmodes(M,[-3 1 -15 15],'AssumeAnalytic',true);
[poles,ip] = cpoles(M,[-3 1 -15 15],'AssumeAnalytic',true);
```

Evaluate with linear solves (`H\B`), not `inv(H)`. Keep the original factors
available during the search for derivatives, visibility and system-matrix
checks. Do not infer minimality from matrix dimensions or a finite set of
samples. State-dependent/nonlinear-time models are not accepted.

### Optional LTI/symbolic adapters

- Symbolic matrices: substitute parameters, preserve quotients/factors and
  exact structural zeros; do not numerically rationalize transcendental terms.
- Delay-free `ss`: preserve realization for mode candidates, but transfer-pole
  output excludes hidden modes. Never claim `eig(A)` alone is the answer.
- `tf/zpk`: extract per-channel factors and exact external delay metadata.
- Internally delayed models: use a faithful supported representation or
  explicitly exploratory evaluation; do not set delays to zero or call Padé.
- Regular descriptor pencils with singular `E` need a separate support gate.
  Until then reject them clearly; infinite structure is not represented by
  arbitrary huge finite numbers.

## 3. Analytic and rank contracts

`AssumeAnalytic=true` asserts analyticity on a neighborhood of the closed
search region of the **named underlying factors**, not of a quotient known
to contain poles. Record exactly which object the assertion concerns:

- `cmodes`: `H`;
- channel pairs: all scalar numerators and denominators;
- `cdyn`: `H,B,C,D`; invertibility of `H` is not asserted everywhere;
- direct transfer zero search: `G`, only if the region is pole-free.

Do not add a universal "assume no difficulties" switch. `Singularities`
continues to denote excluded branch/accumulation points, not ordinary target
poles. A branch cut requires an analytic-domain declaration that avoids
the whole cut; a list of endpoints alone is insufficient. A proposed
`DomainCheck(s)` guard may prohibit evaluations, but cannot prove analyticity.
Derivative stencils, refinement and local validation must respect the region
and guard. Never move the user's outer boundary silently.

Normal rank is the generic maximum rank, not the rank at a zero. Full normal
rank has a numerical nonsingular-minor witness and a dimensional upper bound.
For deficient rank, accept a structural factorization/upper bound or explicit
`NormalRank=r` assertion, checked for contradictions. Sampling a maximum rank
alone does not prove a global upper bound; disclose `normalRankSource`.

## 4. Options

Reuse existing contour/root option names and meanings wherever applicable.
Separate location tolerance, rank tolerance and spectral-count tolerance.

| Option | Proposed role |
|---|---|
| `Output`, `Input` | Explicit scalar channel selection, as above |
| `AssumeAnalytic`, `Singularities` | Existing domain contract with matrix-specific provenance |
| `DomainCheck` | Optional scalar-node logical guard for a declared domain; not a proof |
| `Derivative` | `H'(s)` for `cmodes`, or `G'(s)` for a direct analytic transfer search |
| `NormalRank` | Optional asserted transfer normal rank, integer 0..min(p,m) |
| `Method` | Problem-specific explicit method; first `contour`, later `moments`; no opaque fallback |
| `RootTolerance`, `ContourPoints`, `ContourRefinements` | Inherit current isolation and contour controls |
| `MaxDepth`, `MaxCells`, `MaxIterations` | Existing finite-work controls |
| `RankTolerance` | Relative numerical-rank threshold, default calibrated by rank/scaling fixtures |
| `MaxEvaluations`, `MaxMemoryMB` | Bound matrix work, including probes/moments/local expansions |
| `ProbeSize`, `MomentOrder` | Advanced extraction controls; insufficient capacity remains unresolved |
| `MaxLocalOrder`, `MaxMinors` | Bounds on Laurent/Smith and minor-enumeration work |
| `Scaling` | `none` or fixed row/column equilibration; record applied transformations |
| `SeedPoints` | Optional accelerators; never replace contour coverage evidence |
| `Plot`, `Parent`, `Display`, `Warn` | Preserve existing function-specific graphics/warning conventions |

Numerical defaults other than existing scalar defaults are qualification
decisions, not arbitrary constants to freeze now. Each accepted default
must have a sensitivity/scaling fixture. Reject options that do not apply
to the selected operation, rather than ignoring them.

## 5. Results and diagnostics

Return sorted unique isolated locations as columns, not a copy for every
channel. Repeated orders are in diagnostics, following existing SISO style.
Unresolved clusters are tagged as clusters, not mislabeled repeated roots.

| Field | Contract |
|---|---|
| `kind` | `matrix_characteristic`, `transfer_poles`, `transmission_zeros`, or `channel_zeros/poles` |
| `status`, `complete`, `certified` | Existing vocabulary; certification always false |
| `locationComplete`, `multiplicityComplete`, `countComplete` | Distinct resolved aspects of the requested problem |
| `region`, `unresolvedBoxes`, `locationRadius` | Requested coverage and unresolved cells/clusters |
| `count`, `multiplicity` | Algebraic totals/orders, NaN if unknown; never a channel-count surrogate |
| `partialMultiplicities` | Local Smith-McMillan exponents when resolved; empty otherwise |
| `poleOrder`, `residueRank` | Largest pole order and rank of the simple-pole residue, distinct from multiplicity |
| `normalRank`, `normalRankSource`, `rankHistory` | Transfer rank and evidence/tolerance sensitivity |
| `residuals`, `leftResiduals`, `singularValues` | Scaled matrix/directional checks, with normalization described |
| `rightVectors`, `leftVectors`, `inputDirections`, `outputDirections` | Bases and normalization; empty/uncertain for unsupported degeneracies |
| `channelPresence`, `channelResults` | Which channels contain each pole, with scalar provenance |
| `hiddenModes`, `visibilityUncertain` | Structured modes absent/ambiguous in the transfer at numerical resolution |
| `coincidentLocations`, `clusters` | Do not cancel matrix poles against zeros merely by proximity |
| `analyticSource`, `assumptions`, `warnings` | Precise conditional contracts |
| `evaluations`, `factorizations`, `linearSolves`, `history`, `stopReason` | Reproducible computational work and termination reason |
| `scaling`, `unscaledResiduals` | Mapping back to the original input/output coordinates |

`complete` requires all applicable aspects to be resolved, including normal
rank and local classification. Missing multiplicities keep `count=NaN` and
`complete=false`, even if `locationComplete=true`. Stable contour totals
with unresolved individual points may set `countComplete=true` only for
the counted object; that does not establish transfer counts after reduction.

For zeros of deficient/rectangular matrices, return newly lost directions
relative to the generic nullspace, not arbitrary pre-existing null vectors.
For channel zero functions, report `identically_zero` and an empty set of
isolated zeros; never report an arbitrary mesh as roots. An all-zero matrix
has normal rank zero and no isolated rank-drop problem (`degenerate`).

## 6. Failure behavior and plots

Proposed error identifiers: `ContourRoots:MatrixShape`, `MatrixDomain`,
`MatrixRankContract`, `ChannelSelection`, `MatrixRepresentation`,
`MatrixBudget`. Numerical failure returns available candidates with
`unresolved`/`exploratory` status and an actionable warning; constructor,
dimension and evaluation errors remain errors.

A MIMO `cpzmap` labels "Transfer poles" and "Transmission zeros"; channel
plots include `(output,input)`. A coincident pole/zero displays both markers,
never a scalar "Cancelled" dot. Keep axes/hold/outputs behavior unchanged.
Do not offer global MIMO plots before the underlying zero semantics exist.
