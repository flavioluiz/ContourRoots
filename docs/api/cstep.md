# cstep

Unit-step response of a nonrational transfer function, without a rational
approximation.

## Syntax

```text
cstep(G, t, Name, Value)
y = cstep(G, t, Name, Value)
[y, tOut, info] = cstep(G, t, Name, Value)
```

## Description

`y = cstep(G, t)` returns the zero-state response to a unit step applied at
$t = 0$, at the times in `t`. It inverts $G(s)/s$ numerically on a vertical
line of the complex plane (the Bromwich integral), so only values of $G$
are needed. Without output arguments it plots the response.

The step response is right-continuous: a known direct term $D$ contributes
$D$ from $t = 0$ (or from the transport delay of a `tf` model).

## Inputs

**`G`** — the transfer function:

| Form | Example | Inversion half-plane |
|---|---|---|
| function handle | `@(s) 1./(s + 1 + 0.5*exp(-s))` | you must give `SingularityBound` or `Abscissa` |
| `ndpair` of coefficient vectors | `ndpair([1 2], [1 0.4 4])` | from the poles, automatic |
| `ndpair` with handles | `ndpair(@(s) exp(-s), [1 1])` | you must give it |
| scalar gain | `2` | automatic |
| symbolic expression | `exp(-s)/(s+1)` | you must give it |
| `tf`, `zpk`, `ss` (SISO, continuous) | `exp(-0.3*s)/(s+1)` | automatic; delays are handled exactly |

**`t`** — column or row vector of nonnegative, strictly increasing times.
The default FFT method needs a uniform grid; `'dehoog'` and `'quadrature'`
accept any positive times.

## Main options

| Option | Default | Meaning |
|---|---|---|
| `SingularityBound` | from poles (rational) / required (handles) | $G$ is analytic for $\mathrm{Re}\,s >$ this value |
| `Abscissa` | chosen automatically | the inversion line itself (must be positive for a step) |
| `Method` | `'fft'` | `'fft'`, `'dehoog'` or `'quadrature'` |
| `AbsTol`, `RelTol` | `1e-6`, `1e-4` | convergence tolerance per sample |
| `Feedthrough` | exact, or 0 for handles | the direct term $D$ of $G$ |
| `Plot`, `Parent`, `Warn`, `Display` | | presentation |

All options: [Options and diagnostics](time_response_options.md).

## Outputs

**`y`** — step response at `t` (column). **`tOut`** — the times (column).
**`info`** — diagnostics: `status`, `resolvedMask`, `errorEstimate`,
`abscissa`, `domainSource`, `history`, ...
([all fields](time_response_options.md#diagnostics)).

## Examples

```matlab
G = @(s) 1./(s + 1 + 0.5*exp(-s));
[y, t, info] = cstep(G, (0:0.02:8).', 'SingularityBound', 0);
info.status                                   % 'converged'

y2 = cstep(ndpair(1, [1 -0.3]), (0:0.05:5).');  % unstable: bound from the pole
```

## Errors

| Identifier | Cause |
|---|---|
| `ContourRoots:ResponseDomain` | no bound for a handle, or a bound/abscissa that contradicts known poles |
| `ContourRoots:ResponseTime` | times negative, not increasing, or nonuniform with `'fft'` |
| `ContourRoots:ResponseEvaluation` | $G$ returned a nonfinite value on the line |
| `ContourRoots:ResponseSymmetry` | $G(\bar s) \ne \overline{G(s)}$: not the transfer of a real system |
| `ContourRoots:ResponseBudget` | `MaxPoints` or `MaxMemoryMB` exceeded |

A response that does not converge is returned with a warning
(`ContourRoots:ResponseUnresolved`), not an error.

## See also

[`cimpulse`](cimpulse.md), [`clsim`](clsim.md), [`cinvlaplace`](cinvlaplace.md),
[Tutorial 10](../tutorials/10_time_response.md), MATLAB `step`.
