# Other delay and Padé tools

These functions work with $D(s) + N(s)e^{-sT} = 0$, with coefficient
vectors `N`, `D` in descending powers. Type `help <name>` for details.

## delay_roots

```text
[r, info] = delay_roots(N, D, T, xlim, ylim, nx, ny, Name, Value)
```

Roots in the rectangle `xlim` x `ylim` by damped Newton iterations from an
`nx`-by-`ny` grid, plus Padé poles as starting points; the argument
principle (`delay_root_count`) checks the number found
(`info.countMatches`, false when the count is inconclusive). This is the original, delay-specific solver, kept
for compatibility. For new work, `croots` gives the same roots with the
subdivision-based completeness check.

```matlab
[r, info] = delay_roots(2, [1 1], 1, [-8 2], [-20 20], 36, 60);
info.countMatches
```

## delay_root_count

```text
[count, winding, info] = delay_root_count(N, D, T, xlim, ylim, nPerEdge, maxPerEdge)
```

Argument-principle count of the roots in a rectangle, with adaptive
sampling (from `nPerEdge`, default 256, doubling up to `maxPerEdge`,
default $2^{22}$ points per edge). The count is accepted when it is the same
at three consecutive levels and the phase increments are below $\pi/3$;
otherwise `count` is `NaN` (inconclusive) and `info.reason` explains why.
A negative or unresolved winding number is never rounded to zero.

## pade_delay and pade_characteristic

```text
[num, den] = pade_delay(T, n)
c = pade_characteristic(N, D, T, n)
```

Coefficients of the diagonal $[n/n]$ Padé approximation of $e^{-sT}$
(computed with gamma functions; no toolbox needed), and of the polynomial
$D(s)Q_n(sT) + N(s)P_n(sT)$, of degree $\deg D + n$ (only exactly zero
leading coefficients are removed). For $n \gtrsim 10$ the expanded
polynomial is badly conditioned; see [Tutorial 6](../tutorials/06_pade_pitfalls.md).

```matlab
[num, den] = pade_delay(1, 2)          % (1 - s/2 + s^2/12)/(1 + s/2 + s^2/12)
roots(pade_characteristic(2, [1 1], 1.5, 1))
```

## pade_critical_delays

```text
predicted = pade_critical_delays(N, D, Tmax, orders, Name, Value)
```

Critical delays predicted by Padé models of the given orders, at the
crossing frequencies of the exact equation (the only ones where a diagonal
Padé model can cross). The Padé phase is continuous and monotone, so each
admissible branch is solved directly (no sampling grid): enlarging the delay
interval never removes a crossing. Returns a table with `Order`,
`Frequency`, `Delay`, `PhaseDefect` and `Residual`. Options: `Tmin`,
`Tolerance` (`GridSize` is accepted and ignored).

## match_pade_roots

```text
comparison = match_pade_roots(exactRoots, padeRoots, nDominant)
```

Matches each of the `nDominant` first exact roots to the nearest unused
Padé root; returns a table with `ExactRoot`, `PadeRoot` and
`AbsoluteError`.

## See also

[`critical_delays`](critical_delays.md),
[Tutorial 6](../tutorials/06_pade_pitfalls.md).
