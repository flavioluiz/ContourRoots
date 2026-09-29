# ckernel

Prepare the numerical kernels of a system once, then simulate many inputs
with `clsim` without inverting the Laplace transform again.

**Matrix models.** Explicit matrix models prepare a `ContourRootsMatrixKernel` bank for all channels. See [MIMO responses](matrix_responses.md); the scalar snapshot class and behavior described below remain unchanged.
For matrix models, opt into `SharedGrid=true` with FFT or de Hoog to prepare
all channels and both kernel orders on one adaptive spectral grid. Default
`false` keeps independent channel refinement. Shared setup uses the same
numeric snapshot schema and the same `clsim` error accounting. See the
[benchmark and limits](../benchmarks/shared_grid_preparation.md).

## Syntax

```text
K = ckernel(G, t, Name, Value)
[K, info] = ckernel(G, t, Name, Value)
[y, tOut, info] = clsim(K, u, t, Name, Value)
```

## Description

`clsim(G,u,t)` has two stages ([Tutorial 10, Section 10.5](../tutorials/10_time_response.md#105-arbitrary-inputs)):

1. **invert** $G(s)/s$ and $G(s)/s^2$ to obtain the step response $S$ and
   the ramp response $R$ on the grid `t` (expensive: many evaluations of $G$);
2. **combine** delayed copies of $S$ or $R$, weighted by the jumps or slope
   changes of the input (cheap: one FFT convolution).

Only stage 2 depends on the input. `K = ckernel(G, t)` performs stage 1 and
stores the result; `clsim(K, u, t)` then performs only stage 2. The answer
is the same as `clsim(G, u, t)` up to the tolerances, with **no new
evaluation of $G$**. Use it for sine sweeps, parameter studies of the
input, Monte Carlo runs, or comparing several reference signals.

What `K` stores depends on the hold:

| Hold | Kernels stored | Why |
|---|---|---|
| `'foh'` (default) | $S$ and $R$ | $R$ for the slope changes, $S$ for a nonzero initial value $u(0)$ |
| `'zoh'` | $S$ only | a held input is a sum of delayed steps |

The direct term $D$ and a known transport delay (from a `tf`/`ss` model)
are stored too. Each `clsim(K, ...)` call is an independent simulation
from rest; nothing carries over from the previous input.

**What `K` is, and is not.** `K` is a read-only object holding numbers for
**one model, one time grid and one hold**. It does not keep the function
handle `G`: if you change a parameter used by `G`, the old `K` still
describes the old system, so prepare a new one. There is no hidden cache;
reuse happens only when you pass `K` explicitly. The time samples must be
exactly the ones used for preparation (the same number of points is not
enough).

**Is it worth it?** The kernels are prepared at tighter tolerances than a
single `clsim` call needs (see below), so preparation can cost as much as
one or two ordinary calls. The saving grows with the number of inputs. In
the example below, seven inputs take 1.2 s with `ckernel` against 2.8 s
with seven independent calls.

## Two tolerances: kernel and output

The kernels carry a small numerical error. Each input combines those
errors with its own weights: for ZOH the estimated output error is

$$\epsilon_y(t_n) = \sum_{j \le n} |u_j - u_{j-1}|\;\epsilon_S(t_n - t_j),$$

and FOH does the same with $\epsilon_R$ and the slope changes (plus
$|u_0|\,\epsilon_S$). A large input, or one with many sharp changes,
amplifies the kernel error more than a small smooth one. Hence two
separate tolerances:

- **kernel tolerances**, given to `ckernel`: `AbsTol = 1e-8`,
  `RelTol = 1e-6` by default (tighter than usual, to leave room for
  amplification);
- **output tolerances**, given to `clsim(K,...)`: the usual `1e-6` and
  `1e-4`.

Every `clsim(K, u, t)` call recomputes $\epsilon_y$ for **its** input and
reports it in `info.errorEstimate`. If it exceeds the output tolerance, the
samples are marked unresolved and a `ContourRoots:ResponseUnresolved`
warning is issued, exactly as for an ordinary call. The kernels are **never
refined silently**. To fix it, relax the output tolerance if appropriate, or
prepare a new `K` with tighter kernel tolerances. The estimate is usually
conservative: in the validation, the true error was about $10^{-9}$ when the
estimate first exceeded $10^{-6}$.

## Inputs

**`G`** — the same model forms as [`cstep`](cstep.md#inputs): function
handle, `ndpair`, symbolic expression, polynomial vector or SISO `tf`/`ss`/`zpk`.

**`t`** — uniform time grid starting at 0, with at least two points. A row
or a column with the same values is accepted.

**`K`**, **`u`** (in `clsim(K,u,t)`) — the prepared object and the input
samples (or a handle, sampled on `t`), as in [`clsim`](clsim.md#inputs).

## Main options

For `ckernel` (fixed once prepared):

| Option | Default | Meaning |
|---|---|---|
| `Interpolation` | `'foh'` | hold that the kernels are prepared for |
| `AbsTol`, `RelTol` | `1e-8`, `1e-6` | tolerances of the **kernels** |
| `SingularityBound`, `Abscissa`, `AssumeStable` | as in [`cstep`](cstep.md#main-options) | where $G$ is analytic; required for a function handle |
| `Feedthrough` | exact or 0 | known direct term $D$ |
| `Method` | `'fft'` | inversion method: `'fft'`, `'dehoog'`, `'quadrature'` |
| `SharedGrid` | `false` | matrix models only: common FFT/de Hoog grid across channels and step/ramp orders |
| `MaxPoints`, `MaxRefinements`, `MaxMemoryMB` | `2^20`, 8, 256 | effort and storage budgets |
| `Display` | `false` | print the preparation cost |

For `clsim(K, u, t, ...)` only these may be given:

| Option | Meaning |
|---|---|
| `AbsTol`, `RelTol` | **output** tolerances (defaults `1e-6`, `1e-4`) |
| `Plot`, `Parent`, `Warn`, `Display` | presentation, as in `clsim` |
| `MaxPoints`, `MaxMemoryMB` | budgets of the convolution |
| `Interpolation` | allowed only if equal to `K.Interpolation` |

Anything that would change the model or the inversion (`Method`,
`SingularityBound`, `Feedthrough`, ...) is rejected: prepare a new `K`.
The full list of shared options is in [options and diagnostics](time_response_options.md).

## Outputs

**`K.Time`**, **`K.Interpolation`** — the grid and hold of the preparation.

**`K.Info`** (also returned as `info` by `ckernel`) — preparation
diagnostics: `evaluations` (the cost, paid once), `method`, `AbsTol`,
`RelTol`, `feedthrough`, `delay`, `singularityBound`, `domainSource`,
`assumptions`, and the full diagnostics of each kernel in `step` and
`ramp` (`ramp` is empty for ZOH). These properties can be read but not
assigned.

**`info`** from `clsim(K, ...)` — the usual [`clsim` diagnostics](clsim.md#outputs),
with:

| Field | Value on reuse |
|---|---|
| `kernelReused` | `true` (`false` for `clsim(G,...)`) |
| `evaluations` | `0`: no new evaluation of $G$ (the FFT convolution still runs) |
| `kernelPreparationEvaluations` | the one-time cost, `K.Info.evaluations` |
| `errorEstimate`, `resolvedMask`, `converged`, `status` | computed for **this** input |
| `kernelErrorEstimate`, refinement history | describe the preparation |

## Example

<!-- file: examples/time_response/reuse_kernels.m -->
```matlab
%% Several inputs, one preparation: a sine sweep through a delayed loop
% Run setup_contourroots once per session before this script.
% G is analytic in Re(s) > 0 because |s + 1| >= 1 > 0.5 there.

G = @(s) 1./(s + 1 + 0.5*exp(-s));
t = (0:0.02:20).';
K = ckernel(G, t, 'SingularityBound', 0);     % the inversions happen here

w = 0.5:0.25:2;                               % input frequencies (rad/s)
late = t > 10;                                % the transient has died out
measured = zeros(size(w));
for k = 1:numel(w)
    [y, ~, info] = clsim(K, sin(w(k)*t), t);  % no new inversion
    assert(info.converged && info.evaluations == 0)
    ab = [sin(w(k)*t(late)) cos(w(k)*t(late))] \ y(late);  % fit a sinusoid
    measured(k) = norm(ab);
end
predicted = abs(G(1i*w));
table(w.', measured.', predicted.', ...
      'VariableNames', {'omega', 'measured_gain', 'abs_G_jw'})
fprintf('%d transfer evaluations, all during preparation.\n', K.Info.evaluations);
```

The gains measured from the simulated steady state agree with $|G(i\omega)|$
to about $10^{-4}$, the size of the remaining transient
($e^{-1.1\cdot 10}$; the rightmost root is at $\mathrm{Re}\,s \approx -1.10$).

`K` can be saved to a MAT file and loaded later with the same toolbox
version; the function handle does not need to exist anymore. The snapshot
format may change between versions. Version 0.9.0 retains schema 1 and
loads scalar and matrix snapshots saved by 0.8.0.

## Errors

| Identifier | Cause |
|---|---|
| `ContourRoots:KernelUnresolved` | the step or ramp kernel did not converge during preparation (no object is returned, even with `Warn=false`) |
| `ContourRoots:KernelGrid` | `clsim(K,u,t)` with a grid different from `K.Time` |
| `ContourRoots:KernelInterpolation` | `clsim(K,...)` with a hold different from `K.Interpolation` |
| `ContourRoots:KernelOption` | an option that would change the model or inversion on reuse; `Plot`/`Parent` given to `ckernel` |
| `ContourRoots:KernelModel` | an array of kernel objects, or a snapshot from an unsupported version |
| `ContourRoots:ResponseBudget` | storage (`MaxMemoryMB`) or convolution (`MaxPoints`) budget exceeded |

If **preparation fails**, the usual cause is a long horizon combined with
the tight kernel tolerances: a corner of the response at $t = 0$ or at a
delay needs more frequencies than `MaxPoints` allows. Shorten the horizon,
increase `MaxPoints`, or relax the kernel tolerances (and accept looser
output tolerances later). Model and domain errors are as in
[`cstep`](cstep.md#errors).

## See also

[`clsim`](clsim.md), [`cstep`](cstep.md),
[options and diagnostics](time_response_options.md),
[Tutorial 10, Section 10.5.1](../tutorials/10_time_response.md#1051-several-inputs-one-preparation).
