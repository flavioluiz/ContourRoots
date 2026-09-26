# Time-response options and diagnostics

Shared by [`cimpulse`](cimpulse.md), [`cstep`](cstep.md), [`clsim`](clsim.md)
and [`cinvlaplace`](cinvlaplace.md). All responses are real, causal,
single-input single-output and zero-state (the system starts at rest).
[Tutorial 10](../tutorials/10_time_response.md) explains the ideas.

## Options

### Where the inversion line may be placed

| Option | Default | Meaning |
|---|---|---|
| `SingularityBound` | from the poles for coefficient vectors and `tf`/`zpk`/`ss`; **required** for handles and symbolic input | assertion that the transform is analytic for $\mathrm{Re}\,s >$ this value (right of every pole, branch point and cut) |
| `Abscissa` | chosen from the bound | the vertical line $\mathrm{Re}\,s = \sigma$ itself; must exceed the bound, and be positive for steps, ramps and `clsim` |
| `AssumeStable` | `false` | use the imaginary axis for impulse inversion; valid only if the impulse response is integrable (never for unstable or marginal systems) |

These are statements about your model. A bounded `croots`/`cpoles` search
helps choose them but cannot prove them; `info.domainSource` records where
the line came from.

### Accuracy and effort

| Option | Default | Meaning |
|---|---|---|
| `Method` | `'fft'` | `'fft'` (uniform times), `'dehoog'` or `'quadrature'` (any positive times) |
| `AbsTol`, `RelTol` | `1e-6`, `1e-4` | a sample is accepted when successive refinements differ by less than `AbsTol + RelTol*abs(y)` |
| `MaxRefinements` | 8 | maximum refinement rounds (de Hoog: at most 8 degrees) |
| `MaxPoints` | `2^20` | largest FFT or convolution size; evaluations per quadrature integral |
| `MaxMemoryMB` | 256 | memory budget for work arrays |

### Singular parts of the response

| Option | Default | Meaning |
|---|---|---|
| `Feedthrough` | exact for rational/LTI models, otherwise 0 | the direct term $D$ in $G = D + G_r$ (a split of the given $G$, not an extra gain) |
| `RegularImpulse` | `false` | required by `cimpulse` for handles: states that $G - D$ has no Dirac terms |
| `InitialValue` | exact for rational/LTI models, otherwise unknown (`NaN` at $t = 0$) | right limit $f(0^+)$ of the impulse response (`cimpulse`) or of the inverse (`cinvlaplace`) |

### Input and presentation

| Option | Default | Meaning |
|---|---|---|
| `Interpolation` | `'foh'` | `clsim` only: linear (`'foh'`) or held (`'zoh'`) input between samples |
| `Plot` | `true` without outputs | draw the response |
| `Parent` | new figure | axes to draw into |
| `Warn` | `true` | warn when a response is unresolved (also when `info` is requested) |
| `Display` | `false` | print a one-line summary |

## Diagnostics

The third output `info` is a structure (it is not a state vector).

| Field | Meaning |
|---|---|
| `status`, `converged` | `'converged'` when every sample passed; otherwise `'unresolved'` |
| `resolvedMask` | per-sample: passed the checks |
| `errorEstimate` | per-sample size of the last refinement differences (an estimate, not a bound) |
| `periodError`, `bandwidthError`, `shiftError` | the three separate differences (period doubled, bandwidth doubled, line shifted); for de Hoog, `bandwidthError` is the change with the degree |
| `certified` | always `false` |
| `method`, `abscissa`, `singularityBound`, `domainSource` | method, line used (maximum over times for de Hoog/quadrature), bound, and its origin |
| `singularTerms` | struct array of known Dirac terms: `time`, `order` (0), `weight` |
| `amplification` | factor $e^{\sigma t}$ applied when undoing the shift (large values reduce accuracy) |
| `symmetryError` | check that $G(\bar s) = \overline{G(s)}$ (real system) |
| `time`, `internalStep`, `fftPeriod`, `bandwidth`, `points`, `evaluations` | grids and work done |
| `history` | one row per refinement: FFT (period, points, period/bandwidth/shift differences); de Hoog (time, degree, period/degree/shift differences); quadrature (time, cutoff, tail/level/shift differences) |
| `stopReason` | why the refinement stopped |
| `interpolation`, `kernelConvention`, `kernelErrorEstimate`, `stepKernelInfo` | `clsim` only: hold, integrated-kernel method, and the kernel checks |
| `inputErrorEstimate` | `NaN`: what the input does between samples is not estimated |
| `assumptions`, `warnings`, `endpointPolicy` | the assumptions made, and how $t = 0$ and jumps were treated |

## When a response is unresolved

Invalid options, domains and model values are **errors**. A response whose
refinement does not converge is **returned with a warning**
(`ContourRoots:ResponseUnresolved`); the samples with
`resolvedMask = false` should not be used. Typical causes:

- a jump or corner at a time the solver does not know (a delay inside a
  function handle): only the sample at that time is affected;
- $t = 0$ for a handle without `InitialValue`;
- a very long horizon for a growing response (large `amplification`);
- a very narrow resonance or a very fast feature: try a finer time grid,
  a larger `MaxPoints`, or `Method='quadrature'` at selected times.

See also [Troubleshooting](../troubleshooting.md).
