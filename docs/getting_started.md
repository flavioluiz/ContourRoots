# Getting started

This page takes you from installation to a first complete root search, in
about ten minutes. No prior knowledge of complex analysis is assumed.

## 1. Install

Use any of the three options of the [README](../README.md#installation-no-git-needed):
the `.mltbx` Add-On (double-click), three lines in the Command Window, or
a ZIP file. If you used the ZIP file or the Command Window, run
`setup_contourroots` once in every new MATLAB session (or `savepath`).

Check that MATLAB finds the toolbox:

```matlab
contourroots
which croots
```

## 2. The problem ContourRoots solves

You have a function $F(s)$ of one complex variable and want the values of
$s$ where $F(s) = 0$. For a polynomial, MATLAB's `roots` gives all of them.
For functions such as

$$F(s) = s^2 + s + 1 + (2s+3)\,e^{-s},$$

which appear as characteristic equations of systems with a time delay,
there are **infinitely many** roots, and `roots` does not apply. You must
say *where* to look. ContourRoots looks inside a rectangle:

```matlab
F = @(s) s.^2 + s + 1 + (2*s + 3).*exp(-s);
region = [-8 2 -20 20];         % -8 < Re(s) < 2 and -20 < Im(s) < 20
r = croots(F, region, 'AssumeAnalytic', true)
```

Three things to notice:

- **The function handle** uses element-wise operators (`.*`, `.^`), so it
  can be evaluated at many points at once. It is faster, but not required.
- **The region** is a rectangle `[xmin xmax ymin ymax]`. To decide
  stability of a linear system you usually want a rectangle that extends a
  little to the right of the imaginary axis (here up to `Re(s) = 2`).
- **`'AssumeAnalytic',true`** tells ContourRoots that `F` has no poles,
  branch cuts or other singularities in (a neighborhood of) the rectangle.
  Functions built from polynomials, `exp`, `sin`, `cos`, `sinh` and `cosh`
  have this property. It is needed because MATLAB cannot inspect what is
  inside a function handle.

## 3. Plot the roots

```matlab
croots(F, region, 'AssumeAnalytic', true, 'Plot', true);
```

or, for more control, plot them yourself:

```matlab
figure
plot(real(r), imag(r), 'x', 'MarkerSize', 10, 'LineWidth', 2)
xline(0, ':'); grid on
xlabel('Re(s)'); ylabel('Im(s)')
```

## 4. Check that the search is complete

Ask for the second output, `info`:

```matlab
[r, info] = croots(F, region, 'AssumeAnalytic', true);
info.status          % 'numerically_complete'
info.multiplicity    % each root is simple here
info.residuals       % |F(r)|, tiny
```

`info.status` is the first thing to look at:

| Status | Meaning | What to do |
|---|---|---|
| `numerically_complete` | The contour counts agree with the roots found. | Nothing. |
| `unresolved` | Some part of the rectangle could not be resolved (e.g. a root on the boundary). | See `info.unresolvedBoxes`; shift or enlarge the rectangle slightly. |
| `exploratory` | The function was not declared analytic, so completeness cannot be checked. | If it is analytic, add `'AssumeAnalytic',true`; otherwise use `ndpair` (next tutorial). |

If you call `croots` with a single output and the result is not complete,
it prints a warning, so a missing root never goes unnoticed silently.

## 5. Try it on a polynomial

Because a polynomial is analytic, no assumption is needed, and the
results agree with `roots`:

```matlab
c = [1 2 5 4];                          % s^3 + 2 s^2 + 5 s + 4
rC = croots(c, [-3 1 -3 3]);
rM = roots(c);
max(arrayfun(@(z) min(abs(rC - z)), rM))   % agreement to about 1e-12
```

## 6. Where next

- [Tutorial 1](tutorials/01_first_roots.md) explains the example above in
  depth, and why this system is unstable.
- [Tutorial 2](tutorials/02_poles_and_zeros.md) introduces transfer
  functions, `cpoles`, `czeros` and `cpzmap`.
- Type `help croots` for the complete list of options.
