# unstable_root_count

Number of roots of $D(s) + N(s)e^{-sT} = 0$ with $\mathrm{Re}\,s > 0$,
for one delay, or `NaN` when the count cannot be resolved.

## Syntax

```text
Z = unstable_root_count(N, D, T)
[Z, info] = unstable_root_count(N, D, T, nPerEdge, maxPerEdge)
```

## Description

Requires $\deg N < \deg D$ (retarded case). The radius $R$ of
[`rhp_root_bound`](rhp_root_bound.md) contains every root with
$\mathrm{Re}\,s \ge 0$ for every delay. The function counts, with the
argument principle ([`delay_root_count`](delay_tools.md#delay_root_count)),
the roots in the rectangle $[0, R+1] \times [-(R+1), R+1]$. The sampling is
adaptive: it starts at `nPerEdge` points per edge (default 256) and doubles,
up to `maxPerEdge` (default $2^{22}$), until the count is stable over three
levels and every phase increment is below $\pi/3$.

The left edge of the rectangle **is** the imaginary axis. A root on it
cannot be counted, and roots extremely close to it need more samples than
any budget allows. This happens at or near a critical delay, for roots that
stay on the axis for every delay (common factors of $N$ and $D$), and for
large delays, where many roots approach the axis. Then `Z` is `NaN`, a
warning `unstable_root_count:Inconclusive` is issued (unless `info` is
requested), and `info.reason` says why. An unresolved count is **never**
reported as zero.

`Z == 0` means asymptotic stability only when no root lies on the imaginary
axis. The bound is a mathematical fact; the count is a floating-point
computation, not an interval-arithmetic proof. For many delays, or large
ones, the crossing formula based on [`critical_delays`](critical_delays.md)
(Tutorial 5, Section 5.5) is exact and much cheaper.

**Outputs:** `Z` (a nonnegative integer or `NaN`) and `info` with fields
`resolved`, `samplesPerEdge`, `reason` and `radius`.

## Example

```matlab
Z = arrayfun(@(T) unstable_root_count(3, [1 0.8 4], T), [0.2 1 2.75 3])
% 0 2 0 2
[Z, info] = unstable_root_count(2, [1 1], 1000, [], 2^16);   % many roots near the axis
Z, info.reason                                             % NaN: inconclusive
```

## See also

[`rhp_root_bound`](rhp_root_bound.md), [`critical_delays`](critical_delays.md).
