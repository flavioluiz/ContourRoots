# cinvlaplace

Numerical inverse Laplace transform of a complete expression $F(s)$.

## Syntax

```text
f = cinvlaplace(F, t, Name, Value)
[f, tOut, info] = cinvlaplace(F, t, Name, Value)
```

## Description

`f = cinvlaplace(F, t)` returns $f(t) = \mathcal L^{-1}\{F\}(t)$ for a real,
causal $f$. Use it when the input has a known transform: $F(s) = G(s)U(s)$.
The inversion half-plane must then lie to the right of the singularities of
**both** $G$ and $U$ (e.g. of the pole at 0 of $U(s) = 1/s$).

The value at $t = 0$ is not computed by inversion (a transform only
determines $f(0^+)$ through a limit). For coefficient vectors and LTI
models it is known exactly; for a function handle it is `NaN` unless you
supply `'InitialValue'`. Use `'Method','dehoog'` or `'quadrature'` for
nonuniform or logarithmic time grids.

## Inputs

**`F`** — function handle, `ndpair`, symbolic expression or LTI model, as
in [`cstep`](cstep.md#inputs). **`t`** — nonnegative increasing times.

## Options and outputs

As in [`cstep`](cstep.md), plus `InitialValue`. All options:
[Options and diagnostics](time_response_options.md).

## Example

```matlab
F = @(s) exp(-0.7*sqrt(s))./s;                % diffusion kernel times a step
td = [0; logspace(-2, 1, 25).'];
[f, ~, info] = cinvlaplace(F, td, 'InitialValue', 0, ...
    'SingularityBound', 0, 'Method', 'dehoog');
info.status                                    % 'converged'
```

## See also

[`cstep`](cstep.md), [`clsim`](clsim.md),
[Tutorial 10](../tutorials/10_time_response.md#108-a-known-input-transform-cinvlaplace).
