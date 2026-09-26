# czeros

Zeros of a scalar transfer function inside a rectangle, after
cancellations.

## Syntax

```text
z = czeros(G, region)
[z, info] = czeros(G, region)
[...] = czeros(G, region, Name, Value)
```

## Description

`z = czeros(G, region)` returns the zeros of $N$ that are not cancelled by
zeros of $D$ of equal or higher order. For a single function handle (not an
`ndpair`), it is the same search as `croots`.

Inputs, options and outputs are those of [`cpoles`](cpoles.md), with the
roles of numerator and denominator exchanged.

## Example

```matlab
G = ndpair(@(s) (s + 1).*exp(-s), [1 3 2]);      % (s+1)e^{-s} / ((s+1)(s+2))
z = czeros(G, [-3 1 -3 3], 'AssumeAnalytic', true)   % empty: s = -1 cancels
```

## See also

[`cpoles`](cpoles.md), [`croots`](croots.md), [`cpzmap`](cpzmap.md),
MATLAB `zero`.
