# cpoles

Poles of a scalar transfer function inside a rectangle, after
cancellations.

## Syntax

```text
p = cpoles(G, region)
[p, info] = cpoles(G, region)
[...] = cpoles(G, region, Name, Value)
```

## Description

`p = cpoles(G, region)` returns the poles of $G = N/D$ in the rectangle.
It finds the zeros of $D$, then counts the zeros of $N$ in a small box
around each one. If the order of $N$ is at least the order of $D$ there,
the candidate is **cancelled** and not returned; otherwise the pole is
returned with the remaining multiplicity.

## Inputs

**`G`** — preferably `ndpair(N, D)` with $N$ and $D$ analytic in the
region. Also accepted: a symbolic quotient (split with `numden`), a
`tf`/`zpk`/`ss` model, or a single function handle. With a single handle,
poles are searched as zeros of `1./G`: the search is exploratory unless
`AssumeAnalytic` is set, which then asserts that `1./G` is analytic, i.e.
that `G` has **no zeros** in the rectangle.

**`region`** — `[xmin xmax ymin ymax]`.

## Name-value options

Same as [`croots`](croots.md#name-value-options).

## Outputs

**`p`** — column vector of poles.

**`info`** — diagnostics; in addition to the fields listed for `croots`:
`cancelledLocations` (zeros of $D$ where some order was removed by
cancellation), `cancellationOrders` (orders removed at each),
`cancellationComplete` (true where the candidate disappeared entirely;
false for a partial cancellation, still listed in `p`),
`cancellationUncertain`, and
`targetCountBeforeCancellations`.

## Examples

```matlab
G = ndpair(@(s) sinh(s/2), @(s) sinh(s));
[p, info] = cpoles(G, [-1 1 -10 10], 'AssumeAnalytic', true);
p.'                                   % +-pi i, +-3 pi i
info.cancelledLocations.'             % 0, +-2 pi i

G = ndpair([1 2], [1 4 3]);           % (s+2)/((s+1)(s+3))
cpoles(G, [-4 1 -1 1])                % -1 and -3
```

## See also

[`czeros`](czeros.md), [`cpzmap`](cpzmap.md), [`ndpair`](ndpair.md),
MATLAB `pole`.
