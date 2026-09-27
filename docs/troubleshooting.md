# Troubleshooting

**A time response warns `ContourRoots:ResponseUnresolved`**
Some samples did not pass the convergence checks; `info.resolvedMask`
shows which. The usual causes and fixes:
- a single sample at a jump or corner the solver does not know (for
  example a delay inside a function handle): expected; that sample is the
  midpoint of the jump, the others are fine. With a `tf` model the jump
  time is known and handled exactly;
- the sample at $t = 0$ of a handle's impulse response: supply
  `'InitialValue'` if you know $g(0^+)$, or start the grid at a positive time;
- a long horizon for a growing response (see `info.amplification`):
  shorten the time interval;
- fast or very lightly damped dynamics: use a finer time grid, increase
  `MaxPoints`, or check selected times with `'Method','quadrature'`.
`'Warn',false` only hides the message; unresolved samples stay unresolved.

**`clsim(K,u,t)` is unresolved although `ckernel` succeeded**
The kernel errors are amplified by the jumps and slope changes of each
input, so a large or fast input can miss the output tolerance. The kernels
are not refined automatically. Prepare a new `K` with tighter `AbsTol`/
`RelTol`, or relax the output tolerance if the default was stricter than
you need. See [`ckernel`](api/ckernel.md#two-tolerances-kernel-and-output).

**`ckernel` fails with `ContourRoots:KernelUnresolved`**
The kernels are prepared at tight tolerances; over a long horizon a corner
at $t = 0$ or at a delay may need more frequencies than the budget. Shorten
the horizon, increase `MaxPoints`, or relax the kernel tolerances.

**`clsim(K,...)` rejects the grid, the hold or an option**
`K` is tied to the exact time samples and hold used in `ckernel`, and to
its model and inversion options. Prepare a new `K` for any of these changes.

**A time response asks for `SingularityBound`**
For a function handle, ContourRoots cannot know where the transfer function
is analytic, and the inversion line must lie to the right of all its
singularities. Supply `'SingularityBound',a` with a justified value
([Tutorial 10, Section 10.3](tutorials/10_time_response.md#103-the-one-assumption-where-is-g-analytic)).
For `cimpulse` with a handle, also state `'RegularImpulse',true` (no
hidden Dirac terms) and give `'Feedthrough'` if $G$ has a direct term.

**`Undefined function 'croots'`**
ContourRoots is not on the path. Run `setup_contourroots` from the
installation folder (e.g. `run("ContourRoots/setup_contourroots.m")`), or
install the `.mltbx` Add-On. With the ZIP installation, the path is reset
when MATLAB restarts unless you call `savepath`.

**Warning `ContourRoots:Exploratory`**
You passed a function handle without `'AssumeAnalytic',true`, or a single
quotient handle. If the function is analytic in the rectangle, add the
option; if it is a quotient, use `ndpair(N, D)`. To keep the exploratory
search but silence the warning, request `info` or pass `'Warn',false`.

**Warning `ContourRoots:Incomplete`, status `unresolved`**
Look at `info.unresolvedBoxes`. The most common cause is a root on an edge
of the rectangle (for instance, a root on the imaginary axis while
`xmax = 0`). Move the edges a little. See
[Diagnostics and limits](diagnostics_and_limits.md) for other causes.

**Error `complex_spectrum:AnalyticContract` with symbolic input**
The expression contains something the automatic check cannot prove
analytic (`sqrt`, `log`, a division inside a function, ...). If the
expression is analytic anyway (e.g. $\cosh\sqrt s$), rewrite it as
numerical handles with a series near removable points, and use
`ndpair(..., ...)` with `'AssumeAnalytic',true`
([Tutorial 7](tutorials/07_distributed_systems.md)).

**Error `complex_spectrum:Parameters`**
A symbolic expression still contains parameters besides $s$. Substitute
numbers with `subs` first.

**Error `complex_spectrum:SingularityInRegion`**
The rectangle contains a point declared in `Singularities`. Choose a
rectangle that excludes it.

**The function returns `NaN` at some point**
Typically $0/0$ at a removable singularity (e.g. $\sinh(s)/s$ at $0$). Use a
series near that point, as in `examples/models/distributed_model.m`.

**The search is slow**
- Vectorize the function handle (`.*`, `./`, `.^`).
- Use a smaller rectangle, or split a very tall one into several.
- Provide `SeedPoints` when sweeping a parameter
  ([Tutorial 3](tutorials/03_model_inputs.md#37-following-roots-as-a-parameter-changes)).
- Provide `'Derivative'` if you have an analytic derivative.

**Roots have tiny imaginary parts, e.g. `1.0000 - 0.0000i`**
The search works in complex arithmetic. For real roots, use `real(r)` after
checking `abs(imag(r))`.

**`cpoles` misses a pole I expected**
Check `info.cancelledLocations`: the pole may be cancelled by a zero of the
numerator at the same place. This is correct input/output behavior; use
`croots` on the denominator (or on the characteristic function) to see all
modes.
