# cpzmap

Pole-zero map of a nonrational transfer function.

## Syntax

```text
cpzmap(G, region)
cpzmap(ax, G, region)
[p, z] = cpzmap(G, region)
[p, z, infoP, infoZ] = cpzmap(G, region)
[...] = cpzmap(..., Name, Value)
```

## Description

`cpzmap(G, region)` computes the poles (`cpoles`) and zeros (`czeros`) of
`G` in the rectangle and plots them: poles as `x`, zeros as `o`, cancelled
locations as small grey dots, the rectangle as a dashed box, and the axes
as dotted lines. The title shows the status of the searches.

`[p, z] = cpzmap(...)` returns poles and zeros without plotting, as `pzmap`
does. Add `'Plot', true` to plot as well. `[p, z, infoP, infoZ]` also
returns both diagnostics structures (and then no incompleteness warning is
issued).

`cpzmap(ax, ...)` plots into the axes `ax`. The current hold state is
respected.

All other name-value options are passed to both searches (see
[`croots`](croots.md#name-value-options)).

## Example

```matlab
G = ndpair(@(s) s + 0.5, @(s) (s + 0.5).*(s + 2).*(s.^2 + 1));
figure
cpzmap(G, [-3 1 -2 2], 'AssumeAnalytic', true)
[p, z] = cpzmap(G, [-3 1 -2 2], 'AssumeAnalytic', true);
```

## See also

[`cpoles`](cpoles.md), [`czeros`](czeros.md), MATLAB `pzmap`.
