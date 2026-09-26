# Tutorial 3 — Describing a model

**Goal:** know every way to give a function to ContourRoots, what each one
guarantees, and which options matter in practice.

## 3.1 The five input forms

| Input | Example | Analyticity | Toolbox |
|---|---|---|---|
| Coefficient vector | `[1 3 2]` | automatic | none |
| Function handle | `@(s) s + exp(-s)` | you declare it (`AssumeAnalytic`) | none |
| Numerator/denominator | `ndpair(@(s) sinh(s/2), @(s) sinh(s))` | automatic for vectors and symbolic factors; declared for handles | none |
| Symbolic expression | `s + exp(-s)` with `syms s` | checked automatically | Symbolic Math |
| LTI model | `exp(-2*s)/(s+1)` with `s = tf('s')` | automatic for `tf` with delays | Control System |

The same forms work in `croots`, `cpoles`, `czeros` and `cpzmap`.

## 3.2 Function handles

A handle is the most flexible form: anything you can compute in MATLAB.

```matlab
F = @(s) s.^2 + 0.8*s + 4 + 3*exp(-0.8*s);
r = croots(F, [-5 1 -20 20], 'AssumeAnalytic', true);
```

Good practice:

- **Vectorize** with `.*`, `./`, `.^`. Non-vectorized handles also work
  (ContourRoots falls back to `arrayfun`), only more slowly.
- **Remove removable singularities.** An expression such as
  $\sinh(\sqrt s)/\sqrt s$ is entire, but evaluating it at $s = 0$ gives
  `NaN`. Replace it by its Taylor series near the problem point, as in
  `examples/models/distributed_model.m` (see [Tutorial 7](07_distributed_systems.md)).
- **Branches of `sqrt` and `.^(1/4)`**: an expression can contain a root
  and still be entire when it depends only on even functions of that root,
  as $\cosh(\sqrt s)$ or $\cos(q^{1/4})\cosh(q^{1/4})$. Check this before
  declaring `AssumeAnalytic`.

## 3.3 Numerator/denominator pairs

For transfer functions, always prefer `ndpair(N, D)` with $N$ and $D$
analytic (see [Tutorial 2](02_poles_and_zeros.md)). Each factor can be a
handle, a coefficient vector or a symbolic expression:

```matlab
G1 = ndpair([1 2], [1 4 3]);                       % (s+2)/((s+1)(s+3))
G2 = ndpair(@(s) exp(-s), [1 1]);                  % e^{-s}/(s+1)
p1 = cpoles(G1, [-4 1 -1 1])                        % no assumption needed
p2 = cpoles(G2, [-4 1 -1 1], 'AssumeAnalytic', true)
```

## 3.4 Symbolic expressions (Symbolic Math Toolbox)

Symbolic input is the most convenient when available. ContourRoots splits
a quotient into numerator and denominator (`numden`), checks that both are
entire, and differentiates them exactly:

```matlab
if license('test', 'Symbolic_Toolbox')
    syms s
    Delta = 1 + s + s^2 + (2*s + 3)*exp(-s);
    [r, info] = croots(Delta, [-8 2 -20 20]);          % no AssumeAnalytic
    p = cpoles(cosh(s)/sinh(s), [-1 1 -10 10]);
end
```

The automatic check is conservative: it accepts polynomials, `exp`, `sin`,
`cos`, `sinh`, `cosh`, sums, products and nonnegative integer powers. For
anything else (`sqrt`, `log`, ...), ContourRoots raises the error
`complex_spectrum:AnalyticContract` and asks you to provide regularized
factors. Substitute numerical values for every parameter first (`subs`):
only the complex variable may remain.

## 3.5 Control System Toolbox models

Continuous-time SISO `tf`, `zpk` and `ss` models are accepted. Transport
delays of `tf` models (`InputDelay`, `OutputDelay`, `IODelay`, or
`exp(-T*s)` factors) are handled exactly, as an analytic numerator factor:

```matlab
if license('test', 'Control_Toolbox')
    s = tf('s');
    p = cpoles(exp(-2*s)/(s + 1), [-3 1 -3 3])        % -1
    z = czeros((s + 2)/(s + 1), [-3 1 -3 3])          % -2
end
```

State-space models with *internal* delays are evaluated with `evalfr` in
exploratory mode: ContourRoots cannot extract an analytic characteristic
function from a generic realization. For a complete search, write the
characteristic function or the pair $N/D$ yourself. MIMO models are not
supported.

## 3.6 Options you will use

| Option | Default | Use it when |
|---|---|---|
| `AssumeAnalytic` | `false` | your handle is analytic in the rectangle |
| `Plot`, `Display` | `false` | you want a quick figure or table |
| `Warn` | `true` | you want to silence the incompleteness warning |
| `SeedPoints` | `[]` | you know approximate roots (from a previous parameter value, or from an approximate model): Newton starts there first |
| `Derivative` | `[]` | you have $F'(s)$ as a handle; otherwise a complex finite difference is used |
| `Singularities` | `[]` | the function has known branch or accumulation points: regions containing them are rejected with an error |
| `RootTolerance` | `1e-7` | two roots closer than this (relative) are reported as one multiple root |

Seeds only help Newton converge faster; they never replace the contour
counts. The remaining options (`ContourPoints`, `ContourRefinements`,
`MaxDepth`, `MaxCells`, `MaxIterations`, `GridSize`) control the numerical
effort and are described in [`complex_spectrum`](../api/complex_spectrum.md).

## 3.7 Following roots as a parameter changes

A common pattern is to sweep a parameter and use the previous roots as
seeds. The contour check is repeated at every step:

```matlab
seeds = [];
region = [-6 1 -30 30];
Ts = linspace(0.5, 2, 16);
alpha = zeros(size(Ts));
for k = 1:numel(Ts)
    F = @(s) s + 1 + 2*exp(-s*Ts(k));
    [r, info] = croots(F, region, 'AssumeAnalytic', true, 'SeedPoints', seeds);
    assert(info.complete)
    alpha(k) = max(real(r));           % rightmost real part in the window
    seeds = r;
end
plot(Ts, alpha, 'o-'), yline(0, ':'), xlabel('T'), ylabel('max Re(s)')
```

Next: [Tutorial 4 — How it works](04_how_it_works.md).
