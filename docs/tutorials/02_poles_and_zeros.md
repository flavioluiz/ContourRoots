# Tutorial 2 — Poles and zeros of transfer functions

**Goal:** understand the difference between the roots of a characteristic
function and the poles of a transfer function, and compute poles and zeros
with `ndpair`, `cpoles`, `czeros` and `cpzmap`.

## 2.1 Modes and poles are not the same thing

A linear system has *modes*: the values of $s$ for which $e^{st}$ solves the
homogeneous equations. They are the roots of a characteristic function
$\Delta(s)$. A transfer function $G(s) = N(s)/D(s)$, from one input to one
output, has *poles*: the values of $s$ where $|G(s)| \to \infty$.

Usually $D = \Delta$, and every mode is a pole. But a mode that the input
cannot excite, or that the output cannot see, is **cancelled** by a zero of
the numerator at the same place: it is a root of $D$ and of $N$, and it is
not a pole of $G$. Both notions are useful:

- the **modes** (`croots` on $\Delta$) decide internal stability;
- the **poles** (`cpoles` on $G$) describe the input/output behavior.

## 2.2 Give numerator and denominator separately

Consider $G(s) = \dfrac{\sinh(s/2)}{\sinh(s)}$. Both functions are entire.
$\sinh(s) = 0$ at $s = i\pi k$, and $\sinh(s/2) = 0$ at $s = 2i\pi k$. So
the candidates with even $k$ cancel, and the poles are $\pm i\pi$,
$\pm 3i\pi$, ... Build the transfer function with `ndpair`:

```matlab
G = ndpair(@(s) sinh(s/2), @(s) sinh(s));
[p, info] = cpoles(G, [-1 1 -10 10], 'AssumeAnalytic', true);
p                           % +-pi i and +-3 pi i
info.cancelledLocations     % 0 and +-2 pi i: zeros of D that are not poles
info.cancellationOrders     % one order cancelled at each
```

`'AssumeAnalytic',true` declares that **both** factors are analytic in the
rectangle. For coefficient vectors and symbolic expressions it is not
needed: ContourRoots checks them itself.

Why not simply pass `@(s) sinh(s/2)./sinh(s)`? A quotient has poles, and the
contour integral of a function with poles counts *zeros minus poles*. A
count of zero could mean "nothing here" or "one zero and one pole". With
the two factors separate, ContourRoots counts the zeros of $D$ and of $N$
independently, and subtracts common orders in small boxes around each
candidate. A single quotient handle is accepted, but the search is then
exploratory.

## 2.3 Zeros

`czeros` does the symmetric job: zeros of $N$ that are not cancelled by
$D$. For $G(s) = \dfrac{(s+1)e^{-s}}{(s+1)(s+2)}$ the zero at $s = -1$
cancels:

```matlab
G = ndpair(@(s) (s + 1).*exp(-s), [1 3 2]);    % denominator s^2 + 3s + 2
z = czeros(G, [-3 1 -3 3], 'AssumeAnalytic', true)   % empty
p = cpoles(G, [-3 1 -3 3], 'AssumeAnalytic', true)   % -2 only
```

Note that the numerator is a handle and the denominator a coefficient
vector: both forms can be mixed.

## 2.4 Partial cancellation

Cancellation is counted with multiplicity. If $D$ has a triple zero where
$N$ has a double zero, one pole remains:

```matlab
G = ndpair(@(s) (s + 0.2).^2, @(s) (s + 0.2).^3);
[p, info] = cpoles(G, [-1 1 -1 1], 'AssumeAnalytic', true);
p, info.multiplicity, info.cancellationOrders     % -0.2, 1 and 2
```

## 2.5 The pole-zero map

`cpzmap` draws poles (x) and zeros (o), like `pzmap`, together with the
search rectangle and the cancelled locations:

```matlab
G = ndpair(@(s) s + 0.5, @(s) (s + 0.5).*(s + 2).*(s.^2 + 1));
figure
cpzmap(G, [-3 1 -2 2], 'AssumeAnalytic', true)
```

With output arguments, `cpzmap` returns the poles and zeros instead of
plotting (as `pzmap` does):

```matlab
[p, z] = cpzmap(G, [-3 1 -2 2], 'AssumeAnalytic', true)
```

## 2.6 Modes that the transfer function hides

A sensor placed at a node of a vibration mode cannot see that mode. In the
heat equation example of [Tutorial 7](07_distributed_systems.md), a sensor
at the middle of a rod with Neumann boundary conditions hides the modes
$-\pi^2$ and $-9\pi^2$:

```matlab
root = fileparts(fileparts(which('croots')));
addpath(fullfile(root, 'examples', 'models'))      % example models
G = distributed_model('heat_neumann', 'Sensor', 0.5);
[p, info] = cpoles(G, [-100 2 -3 3], 'AssumeAnalytic', true);
p.'                                   % 0 and -4 pi^2
info.cancelledLocations.'             % -pi^2 and -9 pi^2
```

The characteristic roots of the rod are all four values; the transfer
function from the boundary input to the sensor has only two poles. If a
hidden mode were unstable, the transfer function would not show it:
**input/output poles do not replace an internal stability analysis.**

## Summary

- Use `ndpair(N, D)` with separate analytic factors.
- `cpoles` returns zeros of $D$ not cancelled by $N$; `czeros` the reverse.
- `info.cancelledLocations` lists the candidates where cancellation
  occurred; `info.cancellationComplete` tells full from partial cancellation.
- `cpzmap` plots everything in one call.

Next: [Tutorial 3 — Describing a model](03_model_inputs.md).
