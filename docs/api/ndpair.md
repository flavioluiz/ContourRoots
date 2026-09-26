# ndpair

Transfer function given by a numerator and a denominator.

## Syntax

```text
G = ndpair(N, D)
```

## Description

`G = ndpair(N, D)` represents $G(s) = N(s)/D(s)$. Each factor can be

- a function handle of one complex variable, e.g. `@(s) sinh(s)`;
- a coefficient vector in descending powers, e.g. `[1 3 2]`;
- a symbolic expression or `symfun` (Symbolic Math Toolbox).

Give $N$ and $D$ as separate analytic functions and **do not divide them
yourself**: the separation is what makes counting and cancellation
detection possible (see [Tutorial 2](../tutorials/02_poles_and_zeros.md)).
If both factors are coefficient vectors or symbolic expressions,
analyticity is established automatically; if either is a handle, pass
`'AssumeAnalytic', true` to the search functions after checking it.

`G` is a structure with fields `Numerator` and `Denominator`. Structures
with these fields created by hand are accepted as well.

## Examples

```matlab
G1 = ndpair(@(s) exp(-s), [1 1]);             % e^{-s}/(s+1)
G2 = ndpair([1 2], [1 4 3]);                  % (s+2)/((s+1)(s+3))
s0 = 0.3 + 0.2i;
G1.Numerator(s0) ./ polyval(G1.Denominator, s0)    % evaluate G1 at s0
```

## Errors

`ContourRoots:ndpair` — a missing factor, an unsupported type, or an
identically zero coefficient vector.

## See also

[`cpoles`](cpoles.md), [`czeros`](czeros.md), [`cpzmap`](cpzmap.md).
