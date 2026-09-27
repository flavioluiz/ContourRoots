# MIMO numerical methods

Status: proposed. Parent document: [implementation plan](../mimo_implementation_plan.md).
Contracts: [API and contracts](api_and_contracts.md).

This is an implementation design, not a claim that the algorithms already
exist. All completeness statements are numerical and conditional on the
declared analytic domain. Arbitrary black-box samples cannot certify an
arbitrary nonrational function.

## 1. Mathematical objects and local orders

Let a meromorphic transfer matrix `G` have normal rank `r`. Near an isolated
point `lambda`, analytic locally invertible row and column transformations
reduce it to a rectangular diagonal block with nonzero entries

$$ (s-\lambda)^{\alpha_1},\ldots,(s-\lambda)^{\alpha_r},
   \qquad \alpha_1\leq\cdots\leq\alpha_r. $$

Analytic nonvanishing factors are absorbed into the transformations. The
remaining rows/columns are zero. Define

$$ m_p=\sum_j\max(-\alpha_j,0),\qquad
   m_z=\sum_j\max(\alpha_j,0). $$

The largest pole order is `max(0,-alpha_1)`, not `m_p`. At a regular analytic
point, a transmission zero is a loss of normal rank. At a pole, do not try
to define it by the rank of a matrix containing `Inf`; use local orders.

For a square, full-normal-rank matrix,

$$ \operatorname{ord}_{\lambda}\det G=m_z-m_p. $$

Consequently, a determinant alone misses the pole and zero at `s=-1` of
`diag(1/(s+1),(s+1)/(s+2))`. Neither the union of channel zeros nor scalar
cancellation of determinant factors gives the general MIMO answer.

For the structured form `G=C*(H\B)+D`, use the system matrix

$$ P(s)=\begin{bmatrix}H(s)&-B(s)\\C(s)&D(s)\end{bmatrix}. $$

Where `H` is invertible, `rank(P)=n+rank(G)`. For square `G`,
`det(P)=det(H)*det(G)`. These identities do not justify classifying every
singularity of `P` as a transmission zero where `H` is singular.

Example: `H=diag(s+1,s+2)`, `B=[1;0]`, `C=[1 0]`, `D=0`. The internal modes
are `-1,-2`, the transfer is `1/(s+1)`, and the system matrix has an invariant
zero at `-2`. The transfer has no finite transmission zero. See the
[MathWorks distinction between invariant and transmission zeros](https://www.mathworks.com/help/control/ref/dynamicsystem.tzero.html).

## 2. Evaluation, domains and scaling

Implement a matrix-only evaluation layer with these responsibilities:

1. Evaluate one permitted scalar complex node; validate fixed shape and type.
2. Preserve the representation and analytic provenance of each factor.
3. Count all evaluations, factorizations and solves against shared budgets.
4. Cache evaluations within one search, never across unidentified changes to
   a function handle or its captured parameters.
5. Respect the outer rectangle and domain guard during contour evaluation,
   derivative estimation, refinement and local expansions.

Never use `abs(det(H))`, `svd(H)`, `conj(s)` or continuously changing QR
pivots as an analytic characteristic function. Singular values are useful
for validation and conditioning, not holomorphic root counting.

Allow fixed nonsingular left/right scaling `Htilde=L*H*R`, selected from
admissible reference samples and held fixed for the search or cell.
Transform derivatives consistently and map vectors back to original
coordinates. Use the same fixed transformations throughout a contour.
Pointwise norm normalization must not become an undocumented analytic
target or contaminate derivative formulas.

Residual normalization needs a nonzero reference scale, for example the
sum of norms of known terms `f_j(s)*A_j`, or a fixed local reference scale
from surrounding regular nodes. Dividing only by `norm(H(lambda))` fails
even for a scalar root; using `max(1,norm(H))` hides global tiny scaling.
Publish the normalization and retain unscaled residuals.

## 3. Counting characteristic values of an analytic matrix

Initial contract: square analytic `H`, regular (`det(H)` not identically
zero), and invertible on the closed integration contour. Then

$$ N=\frac{1}{2\pi i}\oint\operatorname{tr}(H(s)^{-1}H'(s))\,ds
      =\operatorname{wind}(\det H(\Gamma),0). $$

`N` counts characteristic values with algebraic multiplicity. A meromorphic
`H` instead gives a net count; do not accept it through the analytic-only
path without explicitly accounting for its pole divisor.

### Primary count: LU log-determinant winding

- Use pivoted LU, including permutation parity, to compute log magnitude
  and phase without multiplying diagonal entries into an overflowing determinant.
- Adapt contour nodes where phase increments, log-magnitude changes or
  conditioning indicate under-resolution. Reuse samples on nested contours.
- Require integer agreement and two successive converged refinements;
  the current scalar integer-error target `1e-7` is a starting qualification
  criterion, not a universal guarantee.
- A near-singular contour, inconsistent winding or nonfinite evaluation
  does not imply zero roots. Return an unresolved boundary/cell.
- Track linear algebra conditioning and uncertainty even when the rounded
  winding happens to be an integer. Do not round away disagreement.

### Cross-check: trace integral

Use a supplied analytic derivative when available; otherwise guarded
finite differences at interior/permitted nodes. Compute `H\Hprime`, not
`inv(H)*Hprime`. Independently refine the quadrature and compare totals.
Two formulas using the same samples are a diagnostic, not mathematical
certification. If unavailable, record why this cross-check was omitted.

Reuse scalar subdivision ideas, but keep matrix evaluation separate.
Parent and child counts must agree. Shift problematic internal split lines
only within the parent; never change the user's outer region. Outer-boundary
roots remain unresolved with guidance to enlarge or shift the rectangle.

## 4. Extraction and refinement

### Isolated simple characteristic value

Refine `(lambda,v)` with a fixed normalization vector `c` using bordered Newton:

$$ \begin{bmatrix}H(\lambda)&H'(\lambda)v\\c^*&0\end{bmatrix}
\begin{bmatrix}\delta v\\\delta\lambda\end{bmatrix}
=-\begin{bmatrix}H(\lambda)v\\c^*v-1\end{bmatrix}. $$

Use damping/guards to stay in the isolated admissible cell. Validate both
right and left residuals and the conditioning quantity `w'*Hprime*v`.
Do not apply the simple-root sensitivity/residue formula when that quantity
is unresolved or effectively zero. Retain clusters for defective/multiple
values until local structure is determined.

### Block contour moments

For a fixed probe matrix `V` with `ell` columns and a centered/scaled local
coordinate `z=(s-center)/radius`, compute moments

$$ A_k=\frac{1}{2\pi i}\oint z(s)^k H(s)^{-1}V\,ds. $$

Implement a Beyn-type extraction first for well-conditioned simple spectra
with sufficient block rank. Extend through block-Hankel/higher moments and
subdivision when needed. The reference methodology is
[Beyn's contour-integral nonlinear eigenvalue method](https://arxiv.org/abs/1003.1580);
the [systems-theoretic contour-method formulation](https://arxiv.org/abs/2012.14979)
is useful for the later block-Hankel design.

Important implementation requirements:

- A size-`n` nonrational matrix can have more than `n` characteristic values
  inside one rectangle. Never truncate the answer to matrix dimension.
- Probe rank, moment order and numerical SVD rank must reconcile with the
  independent contour count. A rank-deficient probe is not an empty spectrum.
- Increase quadrature/capacity or subdivide when extraction undercounts.
  Budget exhaustion returns the missing count and unresolved cells.
- Use deterministic local probe generation without modifying MATLAB's
  global random stream. Recheck against another probe where appropriate.
- Validate every extracted value against the original `H`, not just the
  reduced moment pencil. Compare extracted algebraic totals with the count.
- Avoid large powers of unscaled `s`; record mapping and quadrature weights.

```text
search_matrix(H, region, budgets):
    validate representation and domain contract
    queue <- region
    while queue is not empty and budgets permit:
        cell <- next cell
        count <- adaptive_matrix_count(cell)
        if count is unreliable:
            subdivide safely or mark unresolved
        else if count == 0:
            record a resolved empty cell
        else:
            candidates <- moments/refinement with adequate capacity
            validate candidates, isolation and local orders
            if candidates reconcile with count:
                accept cell and its diagnostics
            else:
                increase capacity or subdivide or mark unresolved
    preserve unresolved queue and reconcile parent/child totals
    return locations, counts, vectors, assumptions and work history
```

## 5. Transfer poles: candidates, visibility and multiplicity

### Candidate coverage

For channel-factor input, search all nonconstant scalar channels with the
existing reduction/cancellation backend. A transfer pole occurs in at least
one reduced entry. Completeness requires all channel searches to resolve.
Retain channel provenance; nearby points form a cluster until a shared
isolating neighborhood validates one location.

For `cdyn` with analytic `H,B,C,D`, characteristic values of `H` cover all
possible finite transfer poles. Hidden modes must be removed using transfer
evidence. If any other factor is meromorphic, either account for its poles
as additional candidates or reject the analytic-factor contract.

Opaque transfer handles without pole provenance provide exploratory
candidates only; do not infer all poles from a successful zero search of
`det(G)` or random entry sampling.

### Simple poles

At an isolated first-order pole, compute the residue matrix by a local contour:

$$ R=\frac{1}{2\pi i}\oint G(s)\,ds. $$

The pole multiplicity is `rank(R)`, not the number of affected channels.
For a simple characteristic value of analytic `H`, a second expression is

$$ R=\frac{C(\lambda)v\,w^*B(\lambda)}{w^*H'(\lambda)v}. $$

Compare the two when applicable. A regular `D` has zero residue. A mode with
provably vanishing residue is hidden in the simple-mode case; ambiguous
near-zero residues are `visibilityUncertain`, not silently cancelled modes.
Threshold sweeps, scaling and smaller valid radii must support rank decisions.

### Higher-order poles

For `k>=1`, local Laurent coefficients satisfy

$$ A_{-k}=\frac{1}{2\pi i}\oint
             G(s)(s-\lambda)^{k-1}\,ds. $$

Every contour must isolate that point and avoid all other poles and cuts.
Repeat with refined nodes and radii. A few small coefficients do not prove
absence of a higher-order pole: obtain an order bound from analytic channel
factors, structured local counts, or equivalent supplied information.
For higher orders, residue rank alone is insufficient; use Section 7.

## 6. Transmission zeros and candidate coverage

### Analytic square, full-normal-rank transfer

Apply the analytic matrix engine to `G`. Here determinant orders are zero
multiplicities because there are no poles in the region. Validate loss of
rank and original-coordinate input/output directions.

### Meromorphic and structured transfer

Resolve pole-free cells using the analytic path. Analyze every excluded pole
neighborhood separately, including coincident zeros; unexamined holes make
the requested regional answer incomplete.

For square structured input with analytic factors, use `P` from Section 1
as an analytic candidate engine. Classify candidates through the local
transfer, visibility and pole information. Neither subtracting lists of
roots nor calling every `P` root a transmission zero is sufficient for
nonminimal systems. A simultaneous pole and zero must survive as two labels.

### Rectangular and deficient-rank transfer

1. Establish the normal-rank contract `r`. A structural upper bound or
   explicit assertion is needed for deficient rank; finite samples alone
   only provide lower-bound evidence.
2. Select a fixed `r`-by-`r` submatrix, or constant compression `L*G*R`, with
   a nonsingular witness so its determinant is not identically zero.
3. On an analytic region, every rank drop of `G` is a zero of that chosen
   determinant. Search it with complete numerical coverage.
4. Reject projection/minor artifacts using the full matrix's rank relative
   to `r`, with local scaling and nearby normal-rank samples.
5. Validate newly lost input/output directions relative to the pre-existing
   nullspaces. `G*v=0` alone proves nothing for a wide matrix with a generic
   input nullspace.
6. Obtain multiplicities from local matrix structure, not the order of a
   single minor, which may have extra vanishing factors even at a true zero.

Freeze minor indices/compressions for each searched cell. Node-dependent
pivot selection is not an analytic determinant. Additional minors improve
conditioning and artifact rejection, but one nontrivial fixed minor with
resolved coverage can cover all rank-drop locations in an analytic cell.
Changing minors across cells requires complete coverage bookkeeping.

For meromorphic regions, first account for poles; a net winding of a minor
does not count all its zeros. Rank-deficient matrices use the `r`th singular
value, not their always-zero smallest singular value, as a rank-drop check.

## 7. Local Smith-McMillan analysis: separate research gate

Implement this only after simpler paths are validated. A bounded fallback
for small matrices provides both a method and an independent test oracle:

1. Obtain a justified local maximum entry pole-order bound `h`.
2. Form the analytic local matrix `M(s)=(s-lambda)^h*G(s)`. Evaluate it on
   admissible contours or from analytic factors; do not compute `0*Inf`.
3. Estimate Taylor coefficients by Cauchy integrals. Refine nodes/radii and
   record coefficient uncertainty; cancellation and roundoff limit order.
4. For each `k=1..r`, find the smallest vanishing order `d_k` among all
   nonidentically-zero `k`-by-`k` minors of `M`. Set `d_0=0`.
5. Recover `alpha_k=d_k-d_(k-1)-h` and separate positive/negative orders.

This determinantal-divisor enumeration is combinatorial and is not the
default for large dimensions. Bound work with `MaxMinors` and `MaxLocalOrder`.
An apparently zero minor cannot be declared identically zero just because
its first few coefficients vanish. Require structural evidence or return
unknown orders when the available bounds cannot distinguish the cases.

A production coefficient-based rank-staircase/block-Toeplitz algorithm
should replace enumeration after its derivation and termination criteria
are reviewed. Validate it against exhaustive small fixtures and exact
symbolic polynomial examples. This is a gated algorithm choice, not a
license to guess a staircase from one tolerance.

Require stability under quadrature, radius and rank-threshold changes.
When incompatible orders fit the numerical evidence, retain the location
but set multiplicities/counts to unknown. Preserve separate pole and zero
orders at coincident locations. No `minreal`-style proximity cancellation.

## 8. Completion, cost and reproducibility

- Maintain separate characteristic, transfer-pole and transfer-zero counts.
  A resolved `H` count does not establish the visible transfer count.
- Apply the same total evaluation budget across channels, candidate searches
  and local classification; do not reset the budget inside each subroutine.
- Report boundary ambiguity, insufficient moments, unresolved rank, unsupported
  local structure and exhausted budgets as distinct stop reasons.
- Dense factorizations cost approximately `O(n^3)` per matrix node; transfer
  evaluation includes `m` right-hand sides. Stream moment accumulation where
  possible; bound optional stored solves and block-Hankel allocations before
  allocating memory. Benchmark rather than promise large sparse scalability.
- Cache by model instance and exact node within the search. Do not assume two
  different handles with identical printed names represent the same function.
- A nonempty list with small residuals is not complete without regional
  coverage. Conversely, a converged total count does not resolve locations
  or individual orders automatically.

## 9. Private file layout

Names are proposed; consolidate tiny helpers when this improves readability.

| File under `matlab/private/` | Responsibility |
|---|---|
| `matrix_model.m` | Normalize transfer/structured forms and analytic provenance |
| `matrix_options.m` | Operation-specific options, dimensions and budgets |
| `matrix_values.m` | Scalar-node matrix evaluation, guards, shape and work accounting |
| `matrix_channel.m` | Extract one scalar channel without losing factor metadata |
| `matrix_scale.m` | Fixed equilibration and original-coordinate mappings |
| `matrix_count.m` | LU winding, trace cross-check and adaptive contour diagnostics |
| `matrix_moments.m` | Block contour moments and capacity diagnostics |
| `matrix_refine.m` | Bordered refinement and left/right residuals |
| `matrix_solve.m` | Cell queue, subdivision, count conservation and budgets |
| `matrix_poles.m` | Channel/structured candidates and pole classification |
| `matrix_normal_rank.m` | Rank witnesses, assertions and contradiction checks |
| `matrix_system.m` | Analytic system-matrix construction from `cdyn` factors |
| `matrix_zeros.m` | Transmission-zero candidates, full-matrix validation and coverage |
| `matrix_laurent.m` | Isolated local expansions and residue checks |
| `matrix_local_orders.m` | Bounded local partial multiplicities and ambiguity |
| `matrix_diagnostics.m` | Unified matrix result schema and completeness rules |

Keep public constructors thin. Do not put mathematical classification rules
inside plotting code or duplicate the scalar root engine inside each adapter.

## 10. References and implementation provenance

- [Guttel and Tisseur, The nonlinear eigenvalue problem (2017)](https://doi.org/10.1017/S0962492917000034):
  background and terminology for analytic matrix eigenproblems.
- [Beyn, An integral method for solving nonlinear eigenvalue problems](https://arxiv.org/abs/1003.1580):
  contour-based extraction design reference.
- [Brennan, Embree and Gugercin, Contour Integral Methods for Nonlinear Eigenvalue Problems: A Systems Theoretic Approach](https://arxiv.org/abs/2012.14979):
  block-moment/reduced-system perspective.
- [MathWorks `tzero`](https://www.mathworks.com/help/control/ref/dynamicsystem.tzero.html):
  comparison semantics for finite rational transmission/invariant zeros.
- [Authors' NEP MATLAB examples](https://github.com/nla-group/nep):
  potential independent comparison implementations, not a runtime dependency.

Review licenses and attribution before reusing any third-party code.
Prefer an implementation derived and tested in this repository; optional
reference scripts must not become undeclared mandatory dependencies.
