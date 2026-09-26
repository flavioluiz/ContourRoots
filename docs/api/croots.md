# croots

Roots of a scalar analytic function inside a rectangle.

## Syntax

```text
r = croots(F, region)
[r, info] = croots(F, region)
[...] = croots(F, region, Name, Value)
```

## Description

`r = croots(F, region)` returns the roots of `F` with
`region(1) < Re(s) < region(2)` and `region(3) < Im(s) < region(4)`, as a
column vector ordered by decreasing real part (then increasing imaginary
part). A root of multiplicity $m$ appears once; its multiplicity is in
`info.multiplicity`.

`[r, info] = croots(...)` also returns the diagnostics described in
[Diagnostics and limits](../diagnostics_and_limits.md). When `info` is not
requested and the search is not numerically complete, `croots` issues the
warning `ContourRoots:Exploratory` or `ContourRoots:Incomplete`.

## Inputs

**`F`** — the function, in one of these forms:

| Form | Example | Notes |
|---|---|---|
| function handle | `@(s) s + exp(-s)` | vectorized is faster; declare `AssumeAnalytic` |
| coefficient vector | `[1 0 -1]` | descending powers, as in `roots` |
| `ndpair` | `ndpair(N, D)` | roots of `N` not cancelled by `D` (same as `czeros`) |
| symbolic | `s + exp(-s)` | Symbolic Math Toolbox; checked automatically |
| `tf`/`zpk`/`ss` | `(s+2)/(s+1)` | Control System Toolbox; zeros of the model |

**`region`** — `[xmin xmax ymin ymax]`, with `xmin < xmax` and
`ymin < ymax`. Roots on the boundary are not counted (the status becomes
`unresolved`); choose edges away from roots.

## Name-value options

| Option | Default | Description |
|---|---|---|
| `AssumeAnalytic` | `false` | Declare that a function handle (or both handles of an `ndpair`) is analytic on a neighborhood of the closed rectangle. Enables the completeness check. Not needed for coefficient vectors and checked symbolic input. |
| `Plot` | `false` | Plot the result in a new figure. |
| `Display` | `false` | Print a table of locations, multiplicities, residuals and radii. |
| `Warn` | `true` | Warn when the result is not numerically complete and `info` is not requested. |
| `SeedPoints` | `[]` | Extra starting points for Newton's method. They never replace the contour counts. |
| `Derivative` | `[]` | Handle to $F'(s)$ for a single function handle. Otherwise a complex finite difference is used. |
| `Singularities` | `[]` | Known branch or accumulation points. A region containing one is rejected with an error. |
| `RootTolerance` | `1e-7` | Relative half-width of the box that isolates each location; roots closer than this form one cluster. |
| Advanced | | `ContourPoints` (32), `ContourRefinements` (7), `MaxDepth` (24), `MaxCells` (5000), `MaxIterations` (80), `GridSize` ([25 41]): see [complex_spectrum](complex_spectrum.md). |

## Outputs

**`r`** — column vector of complex locations.

**`info`** — structure; the most used fields are `status`, `complete`,
`multiplicity`, `residuals` and `unresolvedBoxes`.

## Examples

```matlab
F = @(s) s.^2 + s + 1 + (2*s + 3).*exp(-s);
[r, info] = croots(F, [-8 2 -20 20], 'AssumeAnalytic', true);
info.status                           % 'numerically_complete'

croots([1 0 -1], [-2 2 -1 1])         % 1 and -1, like roots([1 0 -1])

[r, info] = croots(@(s) (s - 0.3).^3, [-1 1 -1 1], 'AssumeAnalytic', true);
info.multiplicity                     % 3
```

## Errors

| Identifier | Cause |
|---|---|
| `ContourRoots:Region` | missing or malformed `region` |
| `ContourRoots:Mode` | `'Mode'` passed; use `cpoles`/`czeros` instead |
| `complex_spectrum:AnalyticContract` | symbolic input that cannot be shown analytic |
| `complex_spectrum:SingularityInRegion` | a declared singularity lies in the region |

## See also

[`cpoles`](cpoles.md), [`czeros`](czeros.md), [`cpzmap`](cpzmap.md),
[`complex_spectrum`](complex_spectrum.md), MATLAB `roots`.
