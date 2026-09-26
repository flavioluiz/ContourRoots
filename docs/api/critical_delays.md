# critical_delays

Imaginary-axis crossings of $D(s) + N(s)e^{-sT} = 0$ for all delays in an
interval.

## Syntax

```text
crit = critical_delays(N, D, Tmax)
[crit, info] = critical_delays(N, D, Tmax, Name, Value)
```

## Description

For real coefficient vectors `N` and `D` (descending powers),
`critical_delays` returns every pair $(\omega, T)$ with $\omega > 0$ and
`Tmin <= T <= Tmax` such that $s = i\omega$ is a root. It

1. solves $|D(i\omega)|^2 = |N(i\omega)|^2$, a polynomial in $\omega^2$,
   for the crossing frequencies;
2. generates every delay branch $T = (\theta + 2\pi k)/\omega$ from the
   phase condition;
3. refines each $(\omega, T)$ by a two-variable Newton correction and keeps
   it only if $|F(i\omega, T)|$ is below `ResidualTolerance`;
4. computes the crossing speed $\mathrm{Re}(ds/dT)$.

No Padé approximation and no sweep in $T$ are used. The method is exact up
to floating-point arithmetic for simple crossings. The case $\omega = 0$
(a root at $s = 0$) does not depend on $T$ and is reported in
`info.zeroRootForAllDelays` instead.

## Options

| Option | Default | Description |
|---|---|---|
| `Tmin` | `0` | lower bound of the delay interval |
| `RootTolerance` | `1e-9` | tolerance for accepting real positive roots of the frequency polynomial |
| `ResidualTolerance` | `1e-9` | accepted residual, relative to $1+|D(i\omega)|$ |

## Outputs

**`crit`** — table sorted by delay, with columns `Frequency`, `Delay`,
`Branch` ($k$), `CrossingSpeed` ($\mathrm{Re}\,ds/dT$), `Direction`
(`"destabilizing"`, `"stabilizing"` or `"tangent/degenerate"`) and
`Residual`.

**`info`** — structure with `magnitudePolynomial` (coefficients in
$x = \omega^2$), `candidateFrequencies`, `zeroRootForAllDelays` and
`delayRange`.

## Example

```matlab
crit = critical_delays(3, [1 0.8 4], 10);        % s^2+0.8s+4+3e^{-sT}
crit(:, {'Frequency', 'Delay', 'Direction'})
```

## Limitations

Single delay, real coefficients, simple crossings. Tangential or multiple
crossings are labelled `"tangent/degenerate"` and need further analysis.
For $\deg N = \deg D$ (neutral systems) the crossings are still computed,
but stability conclusions require extra care.

## See also

[`unstable_root_count`](unstable_root_count.md),
[Tutorial 5](../tutorials/05_time_delay_systems.md).
