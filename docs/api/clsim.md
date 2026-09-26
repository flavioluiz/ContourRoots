# clsim

Response of a nonrational transfer function to a prescribed input, without
a rational approximation (the analogue of `lsim`).

## Syntax

```text
clsim(G, u, t, Name, Value)
y = clsim(G, u, t, Name, Value)
[y, tOut, info] = clsim(G, u, t, Name, Value)
```

## Description

`y = clsim(G, u, t)` returns the zero-state response to the input `u`
sampled at the times `t`. Between samples the input is reconstructed by the
chosen hold:

- `'Interpolation','foh'` (default): linear between samples;
- `'Interpolation','zoh'`: constant on each interval, as a D/A converter.

The response is computed for exactly that reconstructed input. Internally
it is a sum of delayed **step** responses (ZOH, one per jump of the held
input) or **ramp** responses (FOH, one per change of slope), each obtained
by inverting $G(s)/s$ or $G(s)/s^2$. This avoids sampling an impulse
response that may be infinite at $t = 0$, and it matches MATLAB's `lsim`
with the same hold for rational models. Without output arguments it plots
input and output.

The third output is **diagnostics, not states** (unlike `lsim`).

## Inputs

**`G`** — as in [`cstep`](cstep.md#inputs).

**`u`** — vector of input samples, one per time, or a function handle
`u(t)`, which is sampled on `t` (features between samples are not seen).

**`t`** — uniform time grid starting at 0, with at least two points.

## Main options

`Interpolation` above, plus the options of [`cstep`](cstep.md#main-options).
`AbsTol`/`RelTol` apply to the output; the kernel tolerances are tightened
according to the size of the input.

## Outputs

**`y`** — response at `t`. **`tOut`** — times. **`info`** — diagnostics,
including `interpolation`, `kernelErrorEstimate` and the propagated
`errorEstimate`.

## Example

```matlab
G = @(s) exp(-0.7*sqrt(s));                  % diffusion
t = (0:0.02:6).';
u = sin(3*t).*exp(-0.4*t);
[y, ~, info] = clsim(G, u, t, 'SingularityBound', 0, 'Interpolation', 'foh');
info.status                                   % 'converged'
```

## Errors

As in [`cstep`](cstep.md#errors); `ContourRoots:ResponseTime` also when
`t` is nonuniform or does not start at 0.

## See also

[`cstep`](cstep.md), [`cinvlaplace`](cinvlaplace.md),
[Tutorial 10](../tutorials/10_time_response.md#105-arbitrary-inputs), MATLAB `lsim`.
