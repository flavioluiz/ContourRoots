# Time-domain simulation of nonrational transfer functions

Status: implemented and released in ContourRoots 0.4.0.

Prepared on 2026-09-26 against ContourRoots 0.3.0, commit `fc8da3b`.
This plan adds a simulation module without replacing or changing the root
solver. No version bump, release build or publication is part of this plan.
The original design below is retained as a development record. The current
API reference and Tutorial 10 describe the implemented contracts.

### Implementation record (2026-09-26)

- Added `cimpulse`, `cstep`, `clsim`, `cinvlaplace`, dedicated model adapters,
  shifted FFT, de Hoog and direct adaptive Bromwich quadrature. The latter
  was added after late oscillatory tests exposed Q-D conditioning limits
  in double precision; no silent method switching is performed.
- Implemented exact hold identities using integrated kernels, direct terms,
  explicit external delays, unknown-origin diagnostics, resource budgets,
  conjugacy checks and two successive refinement checks.
- Added executable tutorial/API documentation, four quickstarts, a full
  delay/PDE/beam/aeroelastic study and a LaTeX/PDF manual chapter.
- Deliberate differences: input handles use the user's fixed sampled
  interpolant (`inputErrorEstimate=NaN`); automatic input resampling and a
  reusable cross-call kernel cache are deferred. Direct impulse convolution
  remains an independent regression check rather than a second public API.
- Extra explicit options: `RegularImpulse`, `InitialValue`, `MaxMemoryMB`.
  The numerical methods expose unresolved cases and do not certify global
  causality, tail bounds or unseen resonances. The root solver is unchanged.
- Verification: all 86 unit/regression/documentation tests passed in MATLAB
  R2023b. With Control/Symbolic paths removed, the new module's 23 base-only
  tests passed and three optional adapter tests were skipped. The full study
  passed, and the rebuilt 56-page manual's affected pages were rendered and
  visually checked. No version change, commit or release build was performed.

## 1. Objective and scope

Compute the zero-state response of a causal, continuous-time LTI system
directly from its nonrational transfer function:

$$Y(s)=G(s)U(s),\qquad y(t)=\mathcal L^{-1}\{G(s)U(s)\}(t).$$

The primary workflow samples the transfer function on a vertical line in
the Laplace plane, uses an FFT-based inversion, and applies Duhamel's
convolution to a user-specified input. It introduces frequency truncation,
sampling and quadrature errors, but no obligatory Pade approximation,
modal truncation or finite-state realization.

The output is evaluated on a time grid; this is not an ODE time-marching
method. Do not describe the numerical result as exact or discretization-free.

### First public simulation release

- SISO, real-valued time responses, with zero structural/internal states and
  zero prehistory, including delay and aerodynamic memory.
- Explicit time vectors and prescribed external inputs.
- Stable, marginal and unstable models when a valid inversion half-plane
  is supplied or established for a supported model class.
- Handles and `ndpair` as primary inputs; optional symbolic and continuous
  SISO `tf`/`ss`/`zpk` adapters, without rationalizing nonrational terms.
- Base MATLAB runtime. Optional toolboxes are adapters and independent
  test oracles, not core dependencies. Preserve R2023b compatibility.
- Numerical convergence diagnostics, not certified error bounds.

### Out of scope initially

Nonzero initial conditions, state trajectories, MIMO, time-varying/nonlinear
models, discrete-time systems, Simulink blocks, feedback against a live
external simulator, arbitrary distributional outputs, automatic global
stability proofs, and automatic time-window selection for unknown models.
An arbitrary complex expression need not define a causal transfer function.

## 2. Public API

| Function | Purpose | Proposed data signature |
|---|---|---|
| `cimpulse` | Unit-impulse response | `[g,tOut,info] = cimpulse(G,t,options...)` |
| `cstep` | Unit-step response | `[y,tOut,info] = cstep(G,t,options...)` |
| `clsim` | Prescribed input response | `[y,tOut,info] = clsim(G,u,t,options...)` |
| `cinvlaplace` | Advanced inverse-transform interface, after core validation | `[f,tOut,info] = cinvlaplace(F,t,options...)` |

Use `clsim`, not `sim` or a second alias `csim`: this matches MATLAB's LTI
terminology and avoids confusion with Simulink. Do not shadow MATLAB names.
The third output is **diagnostics, not a state trajectory**; explicitly
document this difference from MATLAB's `lsim`.

Conventions:

- Without outputs, plot the response and emit any convergence warning.
- With outputs, return column vectors and do not create a figure by default.
- `Plot=true` allows plots together with returned data; accept a parent axes
  through `Parent` without changing unrelated figures or graphics defaults.
- Always warn on an unresolved response by default, even when `info` is
  requested; `Warn=false` suppresses presentation, not diagnostic flags.
- `t` is required in version one. No horizon inferred from a finite pole map.
- `cstep` and `cimpulse` accept strictly increasing, nonnegative times; the
  first FFT implementation requires a uniform grid. Independent inversion
  will later support arbitrary positive evaluation times explicitly.
- `clsim` requires uniform `t` beginning at zero. Reject missing prehistory,
  nonuniform sampling, mismatched lengths, nonfinite values and complex
  time-domain inputs in the first release.
- `u` may be sampled values or a handle `u(t)`. Both are represented using
  an explicit interpolation contract. A handle is sampled, not claimed to
  have been integrated exactly between nodes.

Proposed usage:

```matlab
G = ndpair(@(s) exp(-0.7*s), [1 1]);
t = (0:0.01:8).';
[ys,t,stepInfo] = cstep(G,t,'SingularityBound',0);
[g,~,impulseInfo] = cimpulse(G,t,'SingularityBound',0, ...
    'RegularImpulse',true,'InitialValue',0);

u = sin(2*t).*(t < 3);
[y,~,info] = clsim(G,u,t,'Interpolation','zoh', ...
    'SingularityBound',0);

clsim(G,@(t) sin(2*t).*exp(-t),t,'Interpolation','foh', ...
    'SingularityBound',0);  % plot only

% Advanced path when the input transform is known:
Gfun = @(s) exp(-0.7*s)./(s+1);
Ufun = @(s) 1./s;
[y,~,info] = cinvlaplace(@(s) Gfun(s).*Ufun(s),t, ...
    'SingularityBound',0);
```

Do not advertise these signatures in the main README/API index until their
implementations and executable documentation tests pass.

## 3. Model adapter and domain contract

Add a dedicated response adapter. `spectrum_model` currently converts a
model into a zero/pole search target, may invert a handle, and rejects
degenerate root problems. Those semantics must not be reused as transfer
evaluation semantics.

The response adapter must:

1. Evaluate the actual transfer function, including its gain, external and
   internal delays, and specified branches. Never simulate its determinant
   or denominator when the user requested an input/output response.
2. Accept scalar handles, `ndpair`, scalar symbolic expressions and supported
   continuous SISO LTI objects. Check identical symbolic variable names in
   both factors; reject unresolved parameters.
3. Support scalar-only handles through a bounded scalar-evaluation fallback.
   Validate dimensions, finite results and conjugate symmetry; do not use
   `real(...)` or `ifft(...,'symmetric')` to hide an inconsistent evaluator.
4. Allow a zero-gain model, e.g. `@(s) 0*s`, and constant gains. The existing
   `ndpair(0,D)` constructor rejects zero factors for spectral reasons;
   support zero via the response adapter without changing that constructor
   in this feature. Reject ambiguous raw numeric coefficient vectors; a
   scalar is a constant gain, coefficients belong in `ndpair`.
5. Preserve cancellations in the quotient. Do not silently replace a `0/0`
   sample or an `Inf/Inf` overflow by zero. Use an explicit removable-limit
   evaluator where known, select a valid sampling line, or report failure.
6. Extract exact metadata when available: rational poles, external delay,
   direct feedthrough and properness. Do not guess global properties from a
   few high-frequency samples of an opaque handle.

### Inversion assumptions

A valid vertical line lies in a right half-plane of analyticity and
convergence of the causal transform. Branches must be consistent there;
some analytic functions still do not have an admissible causal inverse.

Proposed options:

- `SingularityBound=a`: a user assertion that the relevant transform is
  analytic for Re(s)>a, with growth compatible with the requested inverse.
  The solver chooses a line to the right with a numerical safety margin.
- `Abscissa=sigma`: an explicit inversion line; for opaque models the user
  is responsible for its validity. If both options are supplied, require
  `sigma>a`. Record that this is an assertion, not a discovered fact.
- `AssumeStable=true`: explicitly permits the unshifted frequency-response
  route only under the documented impulse-integrability conditions. Do not
  enable this by default or equate pole locations with BIBO stability for
  every nonrational system.
- For step and ramp kernels, incorporate the additional zero-frequency
  singularity introduced by `1/s` or `1/s^2`; do not evaluate it at s=0.
- For generic `cinvlaplace`, any supplied bound applies to the full F(s),
  including the input transform, not just to the plant.

If no supported metadata or explicit assertion establishes a usable line,
stop with an actionable error. A search by `cpoles` may detect a contradiction
to a bound, but absence of poles in a finite rectangle cannot validate it.
Do not repurpose the root API's `AssumeAnalytic` flag for this contract.

## 4. Numerical architecture

### 4.1 FFT inversion on a vertical line: primary engine

For F(s), invert F(sigma+i*omega) to obtain exp(-sigma*t)*f(t), then undo
the exponential weighting. With N frequency bins:

$$\Delta\omega=\frac{2\pi}{T_{\rm FFT}},\qquad
\Delta t=\frac{T_{\rm FFT}}{N},\qquad
\omega_{\rm Nyquist}=\frac{\pi}{\Delta t}.$$

For frequency samples in MATLAB FFT order:

$$f_\sigma(t_n)\approx\frac{1}{\Delta t}\operatorname{IFFT}(F_k)_n.$$

Implementation checklist:

- Specify zero, positive and negative frequency bins and even-N Nyquist
  handling explicitly; test normalization and phase with analytical pairs.
- Keep the requested output interval separate from the longer internal FFT
  period. Discard the guard interval, not evidence of failed convergence.
- Refine period and bandwidth separately. Doubling N alone is ambiguous:
  it can increase either the period or the bandwidth depending on dt.
- Shifted inversion reduces periodic aliasing of growing/slowly decaying
  responses, but undoing the shift amplifies numerical errors. Bound the
  amplification, expose it in `info`, and return unresolved on budget or
  dynamic-range exhaustion. Consider blocks for wide time ranges later.
- Estimate errors in the final unweighted output, not only in f_sigma.
- Never claim that zero-padding an existing spectrum adds missing physical
  bandwidth or validates the response. Re-evaluate G at new frequency nodes.
- No default smoothing/window that silently changes the transfer function.
  If filtering is later offered, expose it as a model approximation.

For unstable simulation, prefer convolution of the weighted kernel with
the weighted input, followed by unweighting the output, instead of forming
an exponentially large unweighted kernel first. The identity must be
implemented consistently with the chosen ZOH/FOH input reconstruction;
weighting sample values and silently re-interpolating them is not equivalent.

### 4.2 Impulse response

Invert G(s) after handling any known singular component. Ordinary impulse
responses may be unbounded but locally integrable at t=0, so do not require
a finite initial sample or substitute an arbitrary spike.

If G(s)=D+G_r(s), return the regular response of G_r in `g` and describe
D*delta(t) in `info.singularTerms` with time, order and weight. Plot the
singular component with a labeled arrow, not a mesh-dependent amplitude.
The same metadata representation can cover known delayed impulses from
explicit transport-delay models. This is a documented extension of the
usual sampled-data interface, not a claim that a numeric vector contains
a Dirac distribution.

Opaque handles require an explicit regular-response/feedthrough contract;
sampling alone cannot rule out delayed deltas. `Feedthrough=D` asserts a
split of the supplied complete G, not an extra gain added to it. If the
requested singular structure is not supported or cannot be evaluated,
refuse or mark unresolved. Do not promise impulse simulation for arbitrary
improper transfers or arbitrary impulse trains in the first release.

At discontinuities, numerical inversion can produce midpoint values and
ringing. Use analytical right limits only when established by metadata or
the model contract. Report unresolved t=0 samples explicitly; do not label
the entire vector converged while silently inserting guessed values.

### 4.3 Step response

Compute the step response directly from G(s)/s, not solely by integrating
an already truncated impulse response. This is both more robust and a
cross-check of the impulse route. Treat direct feedthrough and known
transport-delay jumps explicitly. State the right-continuous convention
at known jump times and report numerical uncertainty around unknown jumps.

### 4.4 Arbitrary input: Duhamel with explicit reconstruction

Implement ZOH and FOH input reconstruction. Default to FOH for a supplied
handle or samples unless the user selects ZOH; document that discontinuous
signals should specify their breakpoints and/or select ZOH. No automatic
interpretation of a sampled step or unresolved high-frequency input.

For sampled data, numerical refinement must preserve the original chosen
interpolant. For a handle, finer sampling can also refine its approximation;
distinguish this additional input error in the diagnostics.

An important improvement over a naive dt*conv(g,u) is to use integrated
response kernels. Write G=D+G_r, with S_r=L^-1{G_r/s} and
R_r=L^-1{G_r/s^2}, both zero before time zero.

For ZOH samples at t_j, define u_-1=0. Then

$$y(t)=D u(t)+\sum_{j:t_j\le t}(u_j-u_{j-1})S_r(t-t_j).$$

For FOH samples, let m_j=(u_{j+1}-u_j)/(t_{j+1}-t_j). On the supplied
input interval,

$$y(t)=D u(t)+u_0S_r(t)+m_0R_r(t)
 +\sum_{j\ge1:\ t_j\le t}(m_j-m_{j-1})R_r(t-t_j),$$

where only defined segment slopes enter the sum. These identities are
exact for the reconstructed input and exact kernels; the implementation
still has kernel-inversion error. They avoid sampling singular g(0) and
provide a consistent causal treatment of endpoints and feedthrough.

Retain direct impulse convolution as an independently tested route for
regular kernels and for educational comparison. Use linear convolution,
with FFT length at least L_kernel+L_input-1 when accelerating finite
sequences. Such padding fixes circular convolution, not inversion aliasing.
Account explicitly for whether stored weights are samples of g or already
integrated interval weights; never apply dt twice.

No values of u after the requested final time should influence an earlier
response. Reuse kernels across multiple inputs through an explicit,
validated internal cache; do not rely on persistent hashes of opaque handles.

### 4.5 Independent inversion method

Add an accelerated Bromwich method, initially de Hoog, as an independent
reference and explicit alternative. Review its mathematical formulation,
double-precision limits, source license and attribution before implementing.
Do not blindly copy a multiprecision implementation into MATLAB doubles.
Its internal series acceleration can use a Pade construction; distinguish
that numerical device from replacing the physical transfer function G(s)
by a finite-order rational aerodynamic or delay model.

This path permits sparse/nonuniform positive times and selected-point
cross-checks. Do not install a runtime Python/mpmath dependency. External
high-precision implementations may produce offline reference data.
Do not choose fixed Talbot or Stehfest as universal fallbacks for delay and
highly oscillatory models; a second method can also fail and must report it.

## 5. Options, diagnostics and failure behavior

Keep the initial public option set small:

| Option | Intent |
|---|---|
| `Method` | `fft` initially; `dehoog` once validated; no opaque auto switching |
| `SingularityBound`, `Abscissa`, `AssumeStable` | Explicit inversion-domain contract |
| `AbsTol`, `RelTol` | Requested response convergence, not a certificate |
| `MaxRefinements`, `MaxPoints` | Bound computation and memory |
| `Interpolation` | `zoh` or `foh`, for `clsim` |
| `Feedthrough` | Exact supplied direct-term split where not known |
| `Plot`, `Parent`, `Display`, `Warn` | Presentation, consistent across wrappers |

Allow explicit period/bandwidth overrides only after the adaptive policy is
validated. Do not expose a large set of unstable tuning knobs in version one.

Return a dedicated response `info`, not the root solver's count/complete
schema. Minimum contents:

- `status`: `converged` or `unresolved`, conditional on model assertions.
- `converged`, `certified=false`, `assumptions` and `warnings`.
- `method`, `abscissa`, `singularityBound`, `domainSource`.
- Requested/internal grids, FFT period, bandwidth, points and evaluations.
- Per-sample error estimates and masks, plus absolute/scaled refinement
  differences separated into period, bandwidth, shift and input effects.
- Exponential amplification, symmetry error and nonfinite evaluations.
- Interpolation, kernel convention, singular terms and endpoint policy.
- Refinement history and reason for stopping.

Require agreement over successive refinements of period and bandwidth;
add a line-shift check where possible. Use absTol+relTol*abs(y) scaling
without division by near-zero responses. Stability of two coarse answers
alone is not proof that a narrow resonance has been resolved. Add narrow
resonance tests, use supplied spectral hints, and expose the finite-band
limitation. No silently clipped pre-response, regularized pole, replacement
of NaNs, automatic Pade fallback or spurious success on budget exhaustion.

Invalid domain/inputs are errors. Numerically unresolved results return
available data and diagnostics with a warning; document which entries are
unusable. No display of an unresolved result as a validated simulation.

## 6. Verification and acceptance tests

| Case | Independent reference / property | Main risk |
|---|---|---|
| 1/(s+1) | exp(-t); 1-exp(-t); analytic forced responses | IFFT scale, signs, endpoints |
| Lightly damped second-order model | Closed form; optional MATLAB `impulse/step/lsim` | Narrow resonance, long tails |
| exp(-T*s)/(s+1) | Shifted first-order impulse/step | Phase, non-grid-aligned delay, ringing |
| exp(-T*s) | Delayed input/step, explicit delayed Dirac | Singular impulse policy |
| D+1/(s+1) and constant/zero gain | D*u plus known smooth response | Feedthrough, zero-model handling |
| 1/(s-0.3) | exp(0.3*t), growing step | Invalid unshifted inversion, amplification |
| 1/s and 1/(s^2+w0^2) | Constant/ramp and undamped sine | Imaginary-axis singularities |
| exp(-a*sqrt(s)), a>0 | g=a exp(-a^2/(4t))/(2 sqrt(pi) t^(3/2)); step=erfc(a/(2 sqrt(t))) | Branch, diffusion, algebraic tail |
| 1/sqrt(s) | g=1/sqrt(pi*t); step=2 sqrt(t/pi) | Integrable singularity at zero |
| Delayed-feedback scalar model | Independently derived method of steps | Infinite spectrum, delay feedback |
| Existing heat/beam transfer model | Analytic modal solution where available; independently converged FEM otherwise | PDE input/output mapping |
| Aeroelastic transfer channel | Independent inversion; Jones only as approximation comparison | Branch cut and near-flutter behavior |

Acceptance gates:

1. Analytical smooth cases meet a predeclared mixed error target (initial
   target absTol=1e-6, relTol=1e-4, adjusted only with a documented reason).
   Check entire requested intervals, including amplitude and phase; for
   jumps report one-sided and away-from-jump errors separately.
2. ZOH/FOH responses match the corresponding continuously reconstructed
   input, including u(0) nonzero, pulses, ramps, chirps and smooth handles.
   Optional MATLAB comparisons must explicitly select the same hold.
3. All output samples outside documented singular/jump policies satisfy
   the stated error criteria. Test limitations at very short/long times.
4. Increase period with bandwidth fixed, then bandwidth with period fixed.
   Compare at the same output times and verify diagnostics detect failures.
5. Tests catch circular convolution, incorrect dt factors, conjugacy errors,
   overflow, excessive shift, hidden cancellations, missing-domain contracts,
   negative/duplicate times, unknown options and exhausted budgets.
6. Verify linearity and time-shift behavior, and invariance to input changes
   strictly after a given observation time.
7. Root locations are a qualitative growth/frequency cross-check only.
   Pole sums and Jones simulations are not exact transient oracles for a
   nonrational system with a branch contribution.
8. Base-MATLAB tests pass without optional toolboxes. Optional comparisons
   report skipped explicitly. All pre-existing tests remain passing.

Measure runtime, peak allocation estimates and transfer evaluations against
N and horizon. Enforce a memory budget before allocating FFT/convolution
arrays. Compare finite-sequence `conv` with padded FFT convolution.

## 7. Repository changes when implementation is authorized

```text
matlab/
  cimpulse.m, cstep.m, clsim.m
  cinvlaplace.m                 # expose after independent-engine validation
  private/
    response_model.m            # actual transfer evaluation + metadata
    response_options.m          # options and time/input validation
    response_fft.m              # shifted FFT inversion
    response_refine.m           # period/bandwidth/shift checks
    response_kernels.m          # regular impulse, step and ramp kernels
    response_convolve.m         # hold-consistent linear convolution
    response_dehoog.m           # independent inversion
    response_plot.m             # plots and singular-term annotations
examples/time_response/
  first_step.m
  arbitrary_input.m
  delays_and_diffusion.m
  unstable_response.m
  run_time_response_study.m
tests/unit/test_response_api.m
tests/regression/test_time_response.m
docs/api/                      # one page per new public function
docs/tutorials/10_time_response.md
docs/assets/                   # curated, regenerated comparison figures
```

Use existing private-function scoping and test fixtures. Keep installation
unchanged; do not recursively add private/examples/test directories.
Update `Contents.m`, `contourroots.m`, API index, main documentation index,
validation page, troubleshooting, README and changelog only as functions
become available. Integrate study generation and executable snippets with
the current build tools. Do not break pole/root options or result fields.

The tutorial should distinguish four approximations: rational-model error,
finite frequency band, periodic inversion aliasing, and input interpolation.
Show impulse, step and arbitrary-input plots together with convergence
panels, not only visually plausible response curves. Include at least one
intentional failure of an unshifted inversion for an unstable plant.
Carry the validated chapter into the LaTeX/PDF manual before its release.

## 8. Implementation milestones and release gates

1. **Contract and fixtures:** lock signatures, domain declarations, singular
   response policy and analytical test data. Review against current adapters.
2. **Inversion kernel:** shifted FFT, scaling, guard window, refinement and
   resource limits. Pass analytical smooth/delay/diffusion cases before plots.
3. **Independent oracle:** implement and validate de Hoog; compare selected
   positive times across stable, marginal, unstable and branch-cut cases.
4. **Impulse/step wrappers:** metadata handling, plotting, endpoint semantics
   and optional LTI/symbolic adapters. Do not release with unresolved policies.
5. **Arbitrary input:** ZOH/FOH integrated-kernel formulas, weighted arithmetic,
   causal endpoints, linear convolution and reusable internal kernels.
6. **Existing models and documentation:** heat/delay/beam/aeroelastic examples,
   complete regression suite, reproducible convergence figures and manual.
7. **Release review:** verify no claims of universal causality/stability or
   certified errors, optional-dependency behavior, resource limits, and all
   failure messages. Choose the version and publish only when requested.

Prioritize the three MATLAB-like wrappers as one coherent feature. Defer
`cfreqresp`, Bode/Nyquist plotting, response metrics, a public cache object,
MIMO and nonzero initial conditions rather than expanding the first release.
No `initial` analogue is planned: a transfer function alone does not specify
an arbitrary internal initial state, especially for distributed systems.

## 9. Sources to use during implementation

- MATLAB [`impulse`](https://www.mathworks.com/help/control/ref/dynamicsystem.impulse.html),
  [`step`](https://www.mathworks.com/help/control/ref/dynamicsystem.step.html)
  and [`lsim`](https://www.mathworks.com/help/control/ref/dynamicsystem.lsim.html):
  interface conventions and explicit input interpolation; do not claim full
  compatibility with all MATLAB syntax or its state-output contract.
- MATLAB [`ifft`](https://www.mathworks.com/help/matlab/ref/ifft.html) and
  [linear versus circular convolution](https://www.mathworks.com/help/signal/ug/linear-and-circular-convolution.html):
  normalization, ordering, symmetry and zero-padding.
- [Numerical inverse Laplace transforms, mpmath documentation](https://mpmath.readthedocs.io/en/latest/calculus/inverselaplace.html):
  algorithm descriptions and failure modes, not a runtime dependency.
- de Hoog, Knight and Stokes, *An Improved Method for Numerical Inversion of
  Laplace Transforms*, SIAM Journal on Scientific and Statistical Computing
  3(3), 357-366 (1982), [DOI 10.1137/0903022](https://doi.org/10.1137/0903022). Read the original method before
  implementation and record provenance for any reused code.
