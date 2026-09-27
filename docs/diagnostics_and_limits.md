# Diagnostics and limits

This page explains exactly what ContourRoots guarantees, what it assumes,
and how to probe a result you are not sure about.

## What the result means

When `info.status` is `'numerically_complete'`:

- the function was treated as analytic in a neighborhood of the closed
  rectangle (automatically for polynomials, `ndpair` of polynomials and
  checked symbolic input; by your `AssumeAnalytic` declaration otherwise);
- the argument-principle count on the outer boundary converged: three
  consecutive sampling levels (each doubling the number of points) gave the
  same integer, and at the last two levels the winding number was within
  `1e-7` of that integer and every sampled phase increment was below
  $\pi/3$;
- every cell of the subdivision was resolved, and the multiplicities of the
  returned roots add up to the outer count.

It does **not** mean:

- a mathematical proof (no interval arithmetic is used; `info.certified` is
  always `false`);
- anything about roots outside the rectangle;
- that the declared analyticity is true.

## Every field of `info`

| Field | Meaning |
|---|---|
| `status` | `'numerically_complete'`, `'unresolved'` or `'exploratory'` |
| `complete` | logical form of the above: counts resolved and reconciled |
| `certified` | always `false` in this version |
| `count` | total number of roots (poles/zeros) returned, with multiplicity, after cancellations |
| `targetCountBeforeCancellations` | outer count of the searched factor before cancellations |
| `multiplicity` | multiplicity of each returned location (see *clusters* below) |
| `locationRadius` | half-diagonal of the small box that isolated each location: a resolution estimate, not an error bound |
| `residuals` | absolute value of the searched factor at each location ($F$, or $D$ in pole mode); depends on the scaling of that factor, and is not $\lvert G(p)\rvert$ |
| `cancelledLocations`, `cancellationOrders` | candidates where some order was removed by pole-zero cancellation, and how many orders were removed |
| `cancellationComplete` | `true` where the candidate was removed entirely; `false` for a partial cancellation, which stays in the result with a reduced multiplicity |
| `cancellationUncertain` | `true` where the cancellation test itself could not be resolved |
| `unresolvedBoxes` | rectangles `[xmin xmax ymin ymax]` that could not be resolved |
| `contourResolved` | whether the outer contour count converged |
| `cellsVisited` | number of cells of the subdivision |
| `exploratory` | `true` in exploratory mode |
| `region`, `mode` | the rectangle and `'zeros'`/`'poles'` |
| `analyticSource` | why the function was treated as analytic (grammar check, polynomial, user assertion, ...) |
| `note` | a reminder of these caveats |

## Known failure modes and what happens

| Situation | What ContourRoots does | What you should do |
|---|---|---|
| A root on (or extremely close to) the outer boundary | count does not converge; status `unresolved` | shift or enlarge the rectangle slightly |
| Two roots closer than `RootTolerance*(1+|s|)` | returned as one location with multiplicity 2 | reduce `RootTolerance` if they must be separated |
| A pole and a zero closer than the resolution | treated as cancelled | reduce `RootTolerance`; confirm analytically |
| Overflow of the function on the contour (e.g. $e^{-sT}$ far to the left) | count fails for that contour; status `unresolved` | choose a less extreme rectangle, or rescale the function |
| A very oscillatory function | more contour points are used, up to `ContourPoints*2^ContourRefinements` per edge | increase `ContourRefinements` or split the rectangle |
| Budget exhausted (`MaxDepth`, `MaxCells`) | remaining cells listed in `unresolvedBoxes` | increase the budget or use a smaller rectangle |
| A handle that is not analytic but is declared so | counts are meaningless; results may be wrong **without warning** | never declare analyticity you have not checked |
| A branch point or accumulation point in the region | the count cannot converge near it | declare it with `Singularities` (region rejected), or exclude it |
| An opaque handle without declaration | exploratory search; `complete` is `false` | use `ndpair` or `AssumeAnalytic` |

## Clusters and multiplicity

Near a multiple root, floating-point evaluation cannot distinguish a root
of multiplicity $m$ from $m$ simple roots very close together. ContourRoots
reports a single location with `multiplicity` $m$ whenever $m$ roots fall
in a box of half-width `RootTolerance*(1+|s|)`. Treat the location as the
center of a cluster whose size is about `locationRadius`.

## How to probe a result

1. **Change the rectangle** slightly (shift each edge by a small, irrational
   amount). The roots inside both rectangles must not change.
2. **Change the sampling**: `'ContourPoints',64` or
   `'ContourRefinements',9`. The status and roots must not change.
3. **Enlarge the window** in the direction of interest (e.g. to the right,
   for stability). The rightmost root must stay the same.
4. **Evaluate the original function** at the roots (`info.residuals`) and,
   when possible, compare with an independent method (a modal formula, a
   finite-element model, `critical_delays` for delay equations).

## Scope

- Spectral searches accept scalar functions of one complex variable. For
  a matrix characteristic problem, pass an analytic scalar determinant.
  Explicit matrix models and MIMO time responses are supported separately;
  they do not add matrix spectral searches or transmission zeros.
- A rectangle, not an arbitrary curve.
- One evaluation of the function is done at a time on vectors of points;
  the cost is dominated by function evaluations on contours.
