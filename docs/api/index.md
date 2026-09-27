# Function reference

Type `help <name>` in MATLAB for the same information in short form.

## Main functions

| Function | Purpose |
|---|---|
| [`croots`](croots.md) | Roots of a scalar analytic function in a rectangle (like `roots`). |
| [`cpoles`](cpoles.md) | Poles of a transfer function, after cancellations (like `pole`). |
| [`czeros`](czeros.md) | Zeros of a transfer function, after cancellations (like `zero`). |
| [`cpzmap`](cpzmap.md) | Pole-zero map (like `pzmap`). |
| [`ndpair`](ndpair.md) | Transfer function as a numerator/denominator pair. |
| [`complex_spectrum`](complex_spectrum.md) | Core solver behind all of the above; full option list. |

## Time responses

| Function | Purpose |
|---|---|
| [`cimpulse`](cimpulse.md) | Ordinary impulse response and known Dirac metadata. |
| [`cstep`](cstep.md) | Unit-step response from G(s)/s. |
| [`clsim`](clsim.md) | Zero-state response to sampled/handle input, ZOH or FOH. |
| [`ckernel`](ckernel.md) | Prepare reusable step/ramp kernels for several `clsim` inputs. |
| [`cinvlaplace`](cinvlaplace.md) | Inverse Laplace transform of a complete expression. |
| [Options and diagnostics](time_response_options.md) | Inversion-domain contracts, convergence and limits. |

## Time-delay systems

| Function | Purpose |
|---|---|
| [`critical_delays`](critical_delays.md) | All imaginary-axis crossings of $D(s)+N(s)e^{-sT}=0$ up to a maximum delay, with directions. |
| [`unstable_root_count`](unstable_root_count.md) | Number of roots with $\mathrm{Re}\,s > 0$ for one delay (`NaN` if inconclusive). |
| [`rhp_root_bound`](rhp_root_bound.md) | Radius containing every right-half-plane root, for every delay. |
| [Other delay and Padé tools](delay_tools.md) | `delay_roots`, `delay_root_count`, `pade_delay`, `pade_characteristic`, `pade_critical_delays`, `match_pade_roots`. |

## Utilities

| Function | Purpose |
|---|---|
| `setup_contourroots` | Add ContourRoots to the path for the current session (`'-quiet'` to suppress output). |
| `contourroots` | Print the version and an overview. |
| `contourroots_version` | Version string, e.g. `'0.6.0'`. |
| `characteristic_roots`, `transfer_poles` | Earlier names of `croots` and `cpoles`, kept for compatibility. |

## The `info` structure

All searches return the same diagnostics structure. Its fields are
described in [Diagnostics and limits](../diagnostics_and_limits.md#every-field-of-info).
