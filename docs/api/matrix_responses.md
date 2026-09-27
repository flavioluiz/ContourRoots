# MIMO time responses

New in version 0.7.0.

## Syntax

```text
[Y, t, info] = clsim(M, U, t, Name, Value)
[S, t, info] = cstep(M, t, Name, Value)
[h, t, info] = cimpulse(M, t, Name, Value)
[K, preparation] = ckernel(M, t, Name, Value)
[Y, t, info] = clsim(K, U, t, Name, Value)
```

## Description and inputs

`M` is an explicit [`cmimo` or `cdyn`](matrix_models.md) model. Raw matrix
inputs do not change the scalar API's meaning. `U` must be `Nt`-by-`nu`,
or a function handle returning exactly that shape on the column grid `t`.
No row-vector orientation is guessed. `clsim` starts from rest and uses
the selected FOH/ZOH interpolant, not a nonlinear feedback law.

The scalar inverse-Laplace engines are reused channel by channel. A
bounded, operation-local cache shares matrix evaluations at common nodes,
including across step/ramp preparation. Different channels may need
different adaptive grids; sharing is opportunistic, not a forced common
convergence decision. Structured nodes use one block solve for all inputs.
The cache is discarded after the operation and has a finite memory budget;
evicted nodes may be evaluated again. No cache key is based on printed
function-handle names.

## Options

Shared scalar options retain their meanings. Matrix additions are:

- `AbsTol`: scalar or one value per output (in that output's units).
- `Feedthrough`, `InitialValue`: full `ny`-by-`nu` matrices, consistent
  with any model metadata. `RegularImpulse=true` remains required for
  opaque impulse models even if the direct term is supplied.
- `MaxEvaluations`: optional shared matrix-node budget, default `Inf` (as in
  the scalar API, work is bounded per inversion by `MaxPoints`). `MaxPoints`
  bounds each scalar inversion/convolution; `MaxMemoryMB` reserves space
  for the complete bank/diagnostics and bounds the cache and scalar workspace.

Opaque channels require a common `SingularityBound` or `Abscissa` assertion;
rational channels can supply conservative denominator bounds themselves.
Regional pole counts do not establish a global inversion half-plane.

Kernel preparation defaults are `AbsTol=1e-8`, `RelTol=1e-6`; output
tolerances default to `1e-6`, `1e-4`. A bank prepares **every declared
channel**, regardless of the first input. Channel subsets are represented
by constructing a smaller model; implicit channel selection is not supported.

## Outputs and error accounting

| Call | Numeric output shape |
|---|---|
| `clsim` | `Nt`-by-`ny`, sum of all input contributions |
| `cstep`, `cimpulse` | `Nt`-by-`ny`-by-`nu`, one experiment per input |

Impulse arrays contain only the ordinary part. `info.singularTerms` stores
`output`, `input`, `time`, `order`, `weight` for known Dirac terms.
Unknown impulse right limits remain unresolved, not fabricated zeros.

For `clsim`, propagated channel error estimates are summed per output,
including a floating-point summation allowance. Acceptance compares that
sum with `AbsTol(i)+RelTol*abs(Y(:,i))`, **not** the largest channel signal.
Individual unresolved channels remain visible and conservatively prevent
an aggregate convergence claim. Strong cancellation can therefore require
tighter preparation even when every large channel response looks accurate.

`info.channelInfo` contains the scalar diagnostics. `errorEstimate` and
`resolvedMask` have the output shape. `converged`, `status`, `certified=false`,
`evaluations`, `factorizations`, `linearSolves`, `rhsColumns`, `cacheHits`
describe the operation. A scalar channel's evaluation count is a request
count; the top-level count reports actual shared matrix-node evaluations.
Factorization counters cover explicit structured models, not unknown linear
algebra performed inside opaque transfer callbacks.

## Frozen matrix kernels

`K` is a `ContourRootsMatrixKernel` value object with read-only `Time`,
`Interpolation`, `Size`, `Info`. Its private data contain numeric kernels
and metadata, never the original model handles or evaluation cache. It can
be saved/loaded. Changing a captured model parameter cannot change an
existing bank; prepare another bank explicitly.

Reuse requires the identical grid and hold. Only output tolerances,
convolution budgets and presentation may change. Every reuse recomputes
output errors, with `evaluations=0`, `factorizations=0`, `kernelReused=true`.
The separate scalar `ContourRootsKernel` class and its saved schema are
unchanged. Banks are snapshots for this implementation's schema, not a
promise of cross-version serialization compatibility.

## Example

```matlab
M = cmimo({ndpair(1,[1 1]), ndpair(2,[1 2]); 0, 1});
t = (0:.1:1).';
K = ckernel(M,t);
[Y,~,info] = clsim(K,[1+sin(t),cos(t)],t,'AbsTol',[1e-6 1e-6]);
assert(info.converged && info.evaluations == 0)
assert(max(abs(Y(:,2)-cos(t))) < 1e-12)
[S,~,si] = cstep(M,t);
assert(isequal(size(S),[numel(t) 2 2]) && si.converged)
```

## Errors and limits

The scalar `KernelGrid`, `KernelInterpolation`, `KernelOption`,
`KernelUnresolved`, `ResponseDomain` and `ResponseImpulse` contracts apply.
Matrix budget/domain/shape errors remain errors, even with `Warn=false`.
Unresolved inversion/output errors produce `ContourRoots:ResponseUnresolved`.
No bank is returned after unresolved preparation.

No matrix pole/zero solver, nonlinear feedback, nonzero initial conditions,
impulse derivatives or general singular-descriptor support is added.

## See also

[Models](matrix_models.md), [kernel reuse](ckernel.md),
[hybrid/implicit comparison](../tutorials/12_hybrid_mimo.md).
