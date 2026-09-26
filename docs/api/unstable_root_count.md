# unstable_root_count

Number of roots of $D(s) + N(s)e^{-sT} = 0$ with $\mathrm{Re}\,s \ge 0$,
for one delay.

## Syntax

```text
Z = unstable_root_count(N, D, T)
Z = unstable_root_count(N, D, T, nPerEdge)
```

## Description

Requires $\deg N < \deg D$ (retarded case). The radius $R$ of
[`rhp_root_bound`](rhp_root_bound.md) contains every root with
$\mathrm{Re}\,s \ge 0$ for every delay. The function counts, with the
argument principle, the roots in the rectangle $[0, R+1] \times [-(R+1),
R+1]$, sampling `nPerEdge` points per edge (default 4000).

The bound is a mathematical fact; the count is a floating-point contour
computation. It is reliable away from the critical delays; at or very near
a critical delay a root lies on the imaginary axis and the count is
ill-conditioned.

## Example

```matlab
Z = arrayfun(@(T) unstable_root_count(3, [1 0.8 4], T), [0.2 1 2.75 3])
% 0 2 0 2
```

## See also

[`rhp_root_bound`](rhp_root_bound.md), [`critical_delays`](critical_delays.md).
