# Troubleshooting

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
