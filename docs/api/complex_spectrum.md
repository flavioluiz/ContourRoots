# complex_spectrum

The core solver used by `croots`, `cpoles`, `czeros` and `cpzmap`.

## Syntax

```text
[locations, info] = complex_spectrum(model, region, Name, Value)
```

## Description

`complex_spectrum` finds zeros (`'Mode','zeros'`, default) or poles
(`'Mode','poles'`) of a scalar model in the rectangle `region`. The public
functions call it with the appropriate mode; use it directly only when you
need an option combination they do not expose. It never issues the
incompleteness warning.

The algorithm is described in [Tutorial 4](../tutorials/04_how_it_works.md)
and, with proofs, in the [manual](../ContourRoots_manual.pdf).

## Options

| Option | Default | Description |
|---|---|---|
| `Mode` | `'zeros'` | `'zeros'` or `'poles'` |
| `AssumeAnalytic` | `false` | see [`croots`](croots.md) |
| `Derivative` | `[]` | derivative of a single function handle |
| `SeedPoints` | `[]` | extra Newton starting points |
| `RootTolerance` | `1e-7` | isolation box half-width relative to `1+abs(s)` |
| `ContourPoints` | `32` | initial samples per edge of each contour |
| `ContourRefinements` | `7` | maximum number of doublings of the sampling (up to `32*2^7 = 4096` per edge) |
| `MaxDepth` | `24` | maximum subdivision depth |
| `MaxCells` | `5000` | maximum number of cells visited |
| `MaxIterations` | `80` | Newton iterations per start |
| `GridSize` | `[25 41]` | grid of starting points in exploratory mode |
| `Singularities` | `[]` | known branch or accumulation points; region containing one is rejected |
| `Plot`, `Display` | `false` | plot or print the result |

A contour count is accepted when two successive doublings give the same
integer, the computed winding number is within `1e-7` of it, and all
sampled increments of $\arg F$ are below $\pi/3$ and of $\log|F|$ below 2.

## Example

```matlab
[p, info] = complex_spectrum(ndpair(@(s) sinh(s/2), @(s) sinh(s)), ...
    [-1 1 -10 10], 'Mode', 'poles', 'AssumeAnalytic', true, ...
    'ContourPoints', 64);
info.status
```

## See also

[`croots`](croots.md), [`cpoles`](cpoles.md),
[Diagnostics and limits](../diagnostics_and_limits.md).
