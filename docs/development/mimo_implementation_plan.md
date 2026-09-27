# MIMO poles, transmission zeros and matrix characteristic problems

Status: implementation proposal; the APIs described here are not implemented.

Prepared on 2026-09-26. Initial inspection found HEAD `a30f607` (0.4.0)
and uncommitted kernel-reuse/release work declaring version 0.5.0. During
planning, that work was committed separately as `e29c5be` (0.5.0); the final
inspection found only these new planning files untracked. The current code
and documentation must be preserved. This proposal does not authorize a
version bump, commit, release build, or publication.

## Reading order

| Document | Purpose |
|---|---|
| This plan | Scope, decisions, implementation sequence and delivery gates |
| [API and contracts](mimo/api_and_contracts.md) | Functions, model forms, options, outputs and compatibility |
| [Numerical methods](mimo/numerical_methods.md) | Counting, extraction, poles, zeros, multiplicities and failure rules |
| [Validation and acceptance](mimo/validation_and_acceptance.md) | Fixtures, tests, quantitative criteria and release checklist |
| [Review](mimo/review.md) | Independent check of the fixtures, one correction, staging recommendations |

These documents form one plan. Examples are proposed interfaces, not runnable
documentation for an existing release. Keep them under `docs/development/`,
which the current executable-documentation runner excludes.

## 1. Objective

Find the finite poles and transmission zeros of a nonrational transfer
matrix inside a user-supplied rectangle, while retaining the existing SISO
interfaces and the toolbox's original-function evaluation philosophy.

Also support finite-dimensional, matrix-valued analytic characteristic
problems originating from infinite-dimensional physics:

$$H(s)x=B(s)u,\qquad y=C(s)x+D(s)u,$$
$$G(s)=C(s)H(s)^{-1}B(s)+D(s).$$

The entries may contain delays, hyperbolic functions or consistently chosen
branches of special functions. Do not replace these entries by rational
approximations or a finite modal model to implement the search.

The matrix dimension is finite; the spectrum need not be. Only isolated
finite spectral points in the requested analytic/meromorphic domain are
in scope. Numerical agreement is not interval-arithmetic certification.

## 2. Non-negotiable distinctions

| Object | What it means | Must not be confused with |
|---|---|---|
| Channel zero | Zero of one scalar `G(i,j)`, after scalar cancellation | A transmission zero of the matrix |
| Transfer pole | A pole in at least one reduced transfer entry | Every internal mode of a realization |
| Matrix characteristic value | A spectral parameter `lambda` at which `H(lambda)` loses invertibility | A singular value from `svd`, or necessarily an observable transfer pole |
| Transmission zero | A positive local Smith-McMillan exponent; at analytic points, loss of normal rank | Union of channel zeros or an always-present nullspace |
| Invariant zero of a realization | Rank loss of its system matrix | A transmission zero without a minimality/local-reduction justification |
| Multiplicity | Local matrix pole/zero multiplicity | Number of channels containing the location |

Use distinct public operations and diagnostic labels. In particular:

- `czeros` remains scalar/channel-wise; it must not silently become a
  transmission-zero search when passed a matrix.
- `ctzeros` explicitly means transmission zeros, not all invariant zeros
  of a supplied nonminimal realization.
- `cmodes` searches `H`; `cpoles` searches the input/output transfer.
- A matrix determinant's order is a net zero-minus-pole order. It can hide
  a pole and a zero at the same location in different directions.

See the mathematical definitions and counterexamples in the companion files.

## 3. Scope and staged support

### Foundation

- Base MATLAB runtime, initially validated on R2023b.
- Explicit `cmimo` and `cdyn` constructors; no ambiguous new meaning for a
  numeric coefficient vector or the scalar `ndpair` constructor.
- Scalar-channel extraction, channel-wise compatibility and transfer pole
  candidate locations from all channels.
- Existing zero channels, constant channels and rank-deficient matrices are
  legal model data, even though scalar root finding on zero is degenerate.

### First coherent MIMO capability

- Small/dense analytic square matrix eigenproblems with `cmodes`.
- Transfer pole locations plus validated simple-pole residue multiplicities.
- Square, full-normal-rank transmission zeros in regions where the
  representation supplies sufficient analytic/pole information.
- Explicit incomplete results for unresolved coincident pole/zero points,
  higher-order poles, deficient normal rank or unsupported rectangular work.
- MIMO pole/transmission-zero plots only once both searches are available.

Do not call the foundation alone "general MIMO poles and zeros". A restricted
release is acceptable only if its supported matrix class is stated in the
README, reference pages, diagnostics and error messages.

### Full planned extension

- Rectangular and rank-deficient transfer matrices, with normal-rank evidence.
- Coincident poles/zeros and higher-order local structure.
- Structured-model mode visibility and local input/output reduction.
- Conditioning-aware matrix scaling, resource budgets and physical examples.

### Not included

MIMO time simulation, nonlinear feedback, nonzero initial conditions,
infinite-dimensional operator discretization, arbitrary singular descriptor
systems, infinite eigenvalues/zeros, neutral-system essential spectra,
automatic branch-cut discovery, and a promise of all roots in the plane.
Large sparse/distributed-memory eigensolvers are a later optimization.

## 4. Repository integration strategy

The current scalar stack is
`cpoles/czeros -> cr_search -> complex_spectrum -> spectrum_model/spectrum_solve`.
It assumes scalar values and uses analytic numerator/denominator factors
for cancellation checks. Keep it as the scalar backend.

Add matrix dispatch before entering that backend. Do not make
`spectrum_values` accept arbitrary matrix outputs and hope that the scalar
algorithm generalizes. Reuse region validation, plotting conventions,
scalar factor searches and test infrastructure where their contracts fit.

Proposed new public files:

```text
matlab/cmimo.m          transfer-matrix constructor
matlab/cdyn.m           H,B,C,D structured-model constructor
matlab/cmodes.m         analytic matrix characteristic values
matlab/ctzeros.m        transmission zeros
```

Extend `cpoles`, explicit-channel `czeros`, and later `cpzmap`; retain every
existing SISO signature and diagnostic field. Private components are listed
in [Numerical methods](mimo/numerical_methods.md#9-private-file-layout).

## 5. Implementation work packages

### P0. Lock semantics and baseline

1. Re-read the then-current tree and any release work; capture its commit,
   version, dirty-file list and baseline test results.
2. Agree the API table, local multiplicity definition and diagnostic schema.
3. Create hand-derived fixtures before numerical engines.
4. Record this deliberate change in direction: `ROADMAP.md` currently lists
   general MIMO nonlinear eigenproblems as out of scope. Update that wording
   only when implementation begins, distinguishing this bounded matrix
   scope from arbitrary operator problems.

**Exit:** interface tests and mathematical fixtures are reviewed; no existing
behavior has changed. Resolve any future naming collision before coding.

### P1. Models and channel operations

Implement `cmimo`, `cdyn`, dimensions, scalar-node matrix evaluation,
analytic provenance, structural zeros and channel selection. Preserve LTI
delays exactly when representable; do not silently approximate internal delays.

Run the scalar solver on selected channels, with channel identifiers and
per-channel completeness. Deduplicate pole candidates only when local
evidence supports a common location; proximity alone produces a cluster.

**Exit:** adapter/channel fixtures pass with and without optional toolboxes;
all historical SISO tests pass. No transmission-zero claim is made yet.

### P2. Analytic square matrix engine

Implement LU-based log-determinant contour counting, derivative/trace
cross-checks, adaptive subdivision, bordered refinement and matrix residuals.
Add block contour moments (Beyn-type extraction) with probe-rank and moment
capacity diagnostics; subdivide or increase capacity when insufficient.

Expose `cmodes` only after count conservation, repeated-root, boundary and
domain tests pass. Count/extraction discrepancies must remain unresolved.

**Exit:** matrix-polynomial and nonrational characteristic fixtures agree
with independent references; no determinant overflow or false empty regions.

### P3. Transfer poles and simple local multiplicity

Combine channel candidates or structured `H` candidates with transfer
Laurent/residue checks. Distinguish internal modes from visible poles.
Handle residue rank greater than one. Preserve channel provenance and modal
directions where identifiable. Higher-order/ambiguous clusters are explicit
partial results, not forced into the simple-pole rule.

**Exit:** pole fixtures distinguish rank-one versus rank-two residues,
hidden modes and weakly visible modes; completeness flags are honest.

### P4. Square transmission zeros

For analytic full-rank transfer matrices, use the matrix engine directly.
For structured models, use the analytic system matrix as a candidate engine
and classify candidates as transmission/invariant/hidden/uncertain.
For meromorphic input, partition and analyze pole neighborhoods explicitly.

Do not subtract net determinant orders from channel counts as a universal
zero algorithm. Add MIMO `cpzmap` only after zero semantics are validated.

**Exit:** collective zero, channel-only zero and coincidence fixtures pass,
or unsupported local structure is explicitly unresolved. Optional `tzero`
comparisons use known-minimal realizations.

### P5. Rectangular and deficient-rank zeros

Determine/assert normal rank, select fixed analytic minors or compressions
for candidate generation, and validate candidates against the full matrix.
Always-present null directions must not generate zeros. Track candidate
coverage independently from acceptance and multiplicity.

**Exit:** tall, wide and deficient-rank fixtures pass; projection artifacts
are rejected; budget-limited minor searches cannot report completeness.

### P6. Higher-order and coincident local structure

Implement a bounded local Smith-McMillan analysis using contour/Taylor-Laurent
information and rank decisions, checked against small-matrix determinantal
divisors. Resolve aggregate and partial multiplicities separately.
Review this algorithm as an independent numerical feature, not a minor
patch to scalar cancellation logic.

**Exit:** higher-order, defective and simultaneous pole/zero fixtures pass;
uncertain rank staircases return unknown orders, never invented integers.

### P7. Examples, documentation and qualification

Add rational teaching examples, mixed-channel delay/PDE examples and the
existing exact-Theodorsen section. Show channel zeros versus transmission
zeros, modal visibility and the failure of a determinant-only shortcut.
Regenerate the manual only after results are validated and review its pages.

**Exit:** all applicable gates in the validation document pass; new tutorial
snippets run; package contents are inspected only in an explicitly authorized
release workflow. No automatic publication follows implementation.

Dependency order: P0 -> P1 -> P2 -> P3 -> P4 -> P5 -> P6 -> P7.
Fixtures/documentation drafts can accompany each package. No calendar or
version-number promise is attached to the numerical research gates.

## 6. Acceptance policy

- Zero missed/spurious points on well-conditioned reference fixtures.
- Separate location, count, multiplicity and normal-rank confidence.
- `complete=true` only for the advertised problem and full requested region,
  under explicit analytic assumptions, with no unresolved subregions.
- `certified=false` throughout; finite sampling is not a rigorous proof.
- Preserve root/pole coincidence in different directions and all existing
  scalar behavior, including warnings and plotting state.
- No silent rationalization, minimal realization, regularization or dropping
  of weak channels; any numerical reduction has recorded tolerances.
- Optional Control/Symbolic toolboxes are adapters/oracles, not dependencies.

The detailed numeric thresholds, negative tests and release checklist are
in [Validation and acceptance](mimo/validation_and_acceptance.md).

## 7. Documentation standard

Keep the current English documentation style. Each new API page contains
Syntax, Description, Inputs, Main options, Outputs, Examples, Errors and
See also. Add a numbered didactic tutorial after the existing sequence,
with equations, executable assertions and links to reproducible scripts.
Extend the current manual chapters/reference tables rather than replacing
the rewritten time-response documentation. Preserve all ongoing release
edits. Cite original numerical methods and review third-party licenses
before copying implementation code.
