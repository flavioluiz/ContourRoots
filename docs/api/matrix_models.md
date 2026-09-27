# Matrix models: cmimo, cdyn, ceval, cchannel

New in version 0.7.0.

## Syntax

```text
M = cmimo(channelsOrConstantMatrix, Name, Value)
M = cmimo(G, 'Size', [ny nu], Name, Value)
M = cdyn(H, B, C, D, 'Dimensions', [n ny nu], Name, Value)
[G, info] = ceval(M, s, Name, Value)
g = cchannel(M, output, input)
```

## Description

These explicit value objects keep matrix inputs separate from the existing
scalar API. `cmimo` represents a transfer matrix. `cdyn` retains the factors
of $H(s)x=B(s)u$, $y=C(s)x+D(s)u$. Evaluation uses a single block solve
$C(H\backslash B)+D$ per node, sharing the factorization across inputs.
No finite-dimensional or rational approximation is introduced.

## Inputs

- A numeric matrix in `cmimo` is a constant gain, never polynomial coefficients.
- A cell array contains scalar handles, gains, `ndpair` models, scalar
  symbolic expressions or continuous SISO LTI models. Wrap coefficient
  vectors in `ndpair`; explicit zero gains are allowed.
- Matrix handles receive one complex scalar by default. `Size`/`Dimensions`
  are required when they cannot be inferred from constants. Construction
  never probes an arbitrary node such as zero.
- Symbolic matrices require Symbolic Math Toolbox, with parameters
  substituted and one common spectral variable. Continuous MIMO LTI input
  to `cmimo` is split into channels and requires Control System Toolbox.
  Time-response support of each LTI channel follows the scalar adapter;
  singular descriptor representations are not separately qualified.
- `cdyn` factors may independently be numeric, symbolic or handle-valued.
  `D(s)` is a factor, not a declaration of the transfer's high-frequency limit.

## Main options

| Constructor option | Meaning |
|---|---|
| `Size` / `Dimensions` | `[ny nu]` / `[n ny nu]`, positive integers |
| `Batched` | default false; if true, each matrix handle accepts `s(:)` and returns rows-by-columns-by-numel(s) |
| `DomainCheck` | scalar-node logical guard; a rejected node raises an error, not a zero transfer |
| `DomainDescription` | recorded provenance; not a proof of analyticity |
| `Feedthrough` | optional real `ny`-by-`nu` direct terms, before any declared delay |
| `InitialValue` | optional real `ny`-by-`nu` ordinary impulse right limits, before delays |
| `Delay` | additional nonnegative transport delay per channel; delays already in a scalar model are retained too |
| `InputNames`, `OutputNames`, `InputUnits`, `OutputUnits` | optional labels matching the channel counts |

Declared `Delay` multiplies the entire transfer channel, including direct
terms. Do not repeat a delay already present in a handle's expression.
For opaque models, omitted `Feedthrough` is an explicit zero-direct-term
assumption when simulating, just as in the scalar API.

`ceval` accepts `MaxEvaluations` (default 1e6 node evaluations) and
`MaxMemoryMB` (default 256). Batched evaluation is chunked to the available
budget. Sparse scalar-node factors are allowed; returned transfer pages
are dense. Symbolic conversion is scalar-node by default: only declare
batching if the resulting handles actually satisfy the page contract.

## Outputs

`M.Size`, `M.Dimensions`, `M.Representation`, `M.Factors` and `M.Options`
are inspectable, not assignable. Handles can still capture changing state;
this is a model, not a frozen snapshot. No cache survives a call to `ceval`.

`ceval` returns `ny`-by-`nu`-by-`numel(s)` in `s(:)` order. A scalar node
returns a matrix; MATLAB may omit trailing singleton dimensions. Never
infer channel meaning with `squeeze`.

`info` counts node `evaluations`, `factorEvaluations`, `factorizations`,
block `linearSolves` and `rhsColumns`. It does not certify spectral properties.
Solver counters cover toolbox-owned structured solves; operations inside an
opaque transfer handle are not observable. Memory budgets bound toolbox-owned
arrays/workspaces, not arbitrary allocations inside user callbacks.

`cchannel` returns a vectorized scalar evaluator. It does not retain the
model's response metadata in that handle; pass inversion-domain and
direct-term information explicitly when using it separately. Use matrix
response calls to preserve metadata and share evaluation among channels.

## Example

```matlab
H = @(s) [s+1 .1*exp(-s); 0 s+2];
M = cdyn(H,eye(2),eye(2),zeros(2),'Dimensions',[2 2 2]);
[G,info] = ceval(M,[1+2i 3+4i]);
assert(info.factorizations == 2 && info.rhsColumns == 4)
assert(norm(G(:,:,1)-H(1+2i)\eye(2),'fro') < 1e-12)
g12 = cchannel(M,1,2);
assert(abs(g12(1+2i)-G(1,2,1)) < 1e-12)
```

## Errors and limits

`ContourRoots:MatrixShape`, `MatrixRepresentation`, `MatrixOption`,
`MatrixDomain`, `MatrixEvaluation` and `MatrixBudget` distinguish invalid
declarations, forbidden nodes, failed evaluations and resource limits.
Handle errors are not turned into missing channels or zeros.

These constructors do not add matrix dispatch to `croots`, `cpoles`,
`czeros` or `cpzmap`. Matrix modes and transmission zeros are not implemented.

## See also

[MIMO responses](matrix_responses.md), [hybrid tutorial](../tutorials/12_hybrid_mimo.md).
