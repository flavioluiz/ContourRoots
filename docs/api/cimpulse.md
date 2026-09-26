# cimpulse

Impulse response of a nonrational transfer function, without a rational
approximation.

## Syntax

```text
cimpulse(G, t, Name, Value)
g = cimpulse(G, t, Name, Value)
[g, tOut, info] = cimpulse(G, t, Name, Value)
```

## Description

`g = cimpulse(G, t)` returns the zero-state response to a unit impulse at
$t = 0$. Without output arguments it plots the response.

An impulse response may contain Dirac terms: $G(s) = D + G_r(s)$ gives
$D\,\delta(t)$, and a transport delay $\tau$ moves it to $D\,\delta(t-\tau)$.
A vector of samples cannot represent a Dirac, so **`g` is the ordinary part**
$g_r(t)$, and the Dirac terms are listed in `info.singularTerms` (time, order,
weight).

For coefficient vectors and `tf`/`zpk` models, $D$, delays and the right
limit $g(0^+)$ are known exactly. For function handles and symbolic input
they are not:

- `'RegularImpulse',true` is **required**: it states that $G$ has no hidden
  Dirac terms (after subtracting a known `Feedthrough`);
- `'InitialValue',g0` gives $g(0^+)$ if you know it; otherwise `g(1)` is
  `NaN` when `t(1) = 0` (the other samples are unaffected).

Impulse responses that are infinite at zero (e.g. $1/\sqrt{\pi t}$) are
computed at positive times only.

## Inputs

Same as [`cstep`](cstep.md#inputs).

## Main options

`SingularityBound`, `Abscissa`, `Method`, `AbsTol`, `RelTol`, as in
[`cstep`](cstep.md#main-options), plus `RegularImpulse`, `InitialValue` and
`Feedthrough` above, and `AssumeStable` (use the imaginary axis as line,
valid only when the impulse response is integrable). All options:
[Options and diagnostics](time_response_options.md).

## Outputs

**`g`** — ordinary impulse response at `t`. **`tOut`** — times.
**`info`** — diagnostics; `info.singularTerms` holds the Dirac terms.

## Examples

```matlab
g = cimpulse(ndpair(1, [1 1]), (0:0.05:3).');       % exp(-t)

[g, ~, info] = cimpulse(ndpair([1 2], [1 1]), (0:0.1:1).');
info.singularTerms                                  % 1*delta(t): (s+2)/(s+1) = 1 + 1/(s+1)

tp = logspace(-3, 1, 20).';
g = cimpulse(@(s) 1./sqrt(s), tp, 'RegularImpulse', true, ...
    'SingularityBound', 0, 'Method', 'dehoog');     % 1/sqrt(pi*t)
```

## Errors

As in [`cstep`](cstep.md#errors), plus `ContourRoots:ResponseImpulse` when
`RegularImpulse` is missing for a handle, and `ContourRoots:ResponseModel`
for improper rational models (impulse derivatives are not supported).

## See also

[`cstep`](cstep.md), [`clsim`](clsim.md), [`cinvlaplace`](cinvlaplace.md),
[Tutorial 10](../tutorials/10_time_response.md), MATLAB `impulse`.
