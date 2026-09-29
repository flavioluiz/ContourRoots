# Shared-grid MIMO kernel preparation (0.9.0)

`ckernel(M,t,'SharedGrid',true,...)` prepares an explicit `cmimo` or `cdyn`
bank on common adaptive spectral grids. The option is **off by default**;
FFT and de Hoog are supported. Scalar `ckernel`, ordinary `cstep`,
`cimpulse`, fresh `clsim`, and quadrature keep their existing behavior.

```matlab
M = cdyn(@(s) [s+1 .1; 0 s+2],eye(2),eye(2),zeros(2), ...
    'Dimensions',[2 2 2]);
t = (0:.1:.5).';
[K,setup] = ckernel(M,t,'SingularityBound',0,'SharedGrid',true, ...
    'AbsTol',[1e-6 2e-6],'RelTol',1e-4);
[y,~,info] = clsim(K,[1+t sin(t)],t);
assert(setup.sharedGrid && setup.converged && info.converged)
assert(info.evaluations == 0 && setup.rhsColumns == 2*setup.factorizations)
```

## Why a common grid helps

The earlier matrix preparation calls the scalar inverse transform once per
channel and kernel order. Its bounded cache shares whole matrix evaluations
only when complex nodes match exactly. Output-dependent tolerances affect
the integration line, delays affect the inversion horizon, and adaptive
checks can choose different refinement depths. Consequently a model with
many outputs can perform many different factorizations despite the cache.

Shared preparation evaluates the undelayed full transfer matrix at each
requested node, subtracts each channel's direct term, and applies `1/s` and
`1/s^2` for step and ramp. A structured node uses `C*(H\B)+D`: one
factorization with all input columns as right-hand sides. All output
channels and both orders consume those samples immediately; cache eviction
cannot cause a separate evaluation for each consumer. The operation-local
bounded cache also reuses exact nodes between probes/refinement levels;
evicted nodes can still be reevaluated. No cache or callback is saved.

The common line uses the tightest absolute tolerance and the largest
singularity bound (including the integrator pole at zero). An explicit
`Abscissa` must still satisfy every channel's domain and be positive for
integrated kernels. Each channel retains its own absolute tolerance and
relative error test. Exact direct terms and channel delays are unchanged:
the ordinary inverse is evaluated at `t-delay`, with exact zero and
pre-delay values; delay factors are not approximated spectrally.

FFT shares period, bandwidth and shifted-line probes, and increases the
period whenever any active channel/order requires it. Two successive
checks must pass for every positive-time sample. De Hoog shares degree,
period and shifted-line probes at each output-time index; the period uses
the largest active shifted time. It retains the degree-40 guard and two
successive doubling checks, using the identical Q-D recurrence as the
scalar engine. Channels can require more work when their shifted times
are very different. The error floor uses the largest amplification on the
common grid, conservatively. These remain conditional error estimates,
not analyticity or accuracy certificates.

`info.channelInfo` preserves the error components, refinement history,
resolved masks, domain assumptions, delay and feedthrough. `info.sharedGrid`
identifies the setup path. Top-level `evaluations`, `factorizations`,
`linearSolves` and `rhsColumns` count actual work; per-channel evaluation
counts describe logical inversion requests and must not be summed to infer
factorizations. Shared cell-channel evaluation counts matrix nodes, with
`factorEvaluations` counting all scalar channel evaluations at those nodes.

## Measurements and default decision

Measurements use MATLAB R2026a on the 10-core Apple Silicon Mac, serial
`-singleCompThread` processes. The pre-change code is an immutable export
of commit `9b70d02` inside this worktree. The small and hybrid cases run a full warmup before
their recorded run; the dense case warms 16 block solves. Other authorized MATLAB computations were active on the
machine, so wall times are indicative, not an isolated hardware speedup
claim. Evaluation and factorization counts are deterministic.

The benchmark cases are:

- **Small structured:** 2 states, 2 inputs, 2 outputs, default kernel
  tolerances, `t=0:.1:1`.
- **Hybrid wing:** Goland wing at 120 m/s, force/torque inputs and
  deflection/twist outputs, scaled by `1e6`, `t=0:.01:.2`, absolute and
  relative tolerances `1e-3`, singularity bound 5. This is an opaque
  `cmimo`; internal hybrid factorizations are not counted by `ceval`.
- **Dense structured:** a real dense 300-state stable matrix with poles
  between -20 and -0.5, 2 inputs and 30 differently scaled outputs,
  de Hoog on `t=0:.5:1`, per-output absolute tolerances from `1e-6` to
  `1e-4`, relative tolerance `1e-5`. Every structured node performs a
  dense block solve. An initial dense FFT pilot was stopped during warmup
  for runtime; no completed FFT timing is claimed for this 300-state case.

| Case / path | Distinct nodes | Evaluations | Factorizations | Cache hit rate | Setup [s] | Max output difference |
|---|---:|---:|---:|---:|---:|---:|
| Small / independent | 98,308 | 98,308 | 98,308 | 79.22% | 10.859 | — |
| Small / shared | 73,732 | 73,732 | 73,732 | 49.65% | 2.103 | 4.48e-9 |
| Hybrid / independent | 16,386 | 16,386 | opaque | 93.94% | 12.214 | — |
| Hybrid / shared | 16,386 | 16,386 | opaque | 51.52% | 5.589 | 0 |
| Dense 300 / independent | 16,896 | 16,896 | 16,896 | 92.97% | 302.201 | — |
| Dense 300 / shared | 1,056 | 1,056 | 1,056 | 47.31% | 29.653 | 4.26e-10 |

Differences compare `clsim` outputs on the same sampled inputs and FOH.
For the dense case, an exact modal FOH state update gives maximum output
errors of **4.52e-10 independent** and **2.55e-11 shared**. The dense shared
bank reduces actual factorizations by **16×**. The lower cache hit rate is
expected: shared preparation eliminates redundant channel/order requests
before they reach the cache. Rate means `cacheHits/(cacheHits+evaluations)`.
Hybrid internal linear algebra is opaque: the library reports zero explicit
factorizations for this representation, which does not mean the callback
performs none. Its gain here comes from fewer cache lookups and requests,
not fewer distinct transfer evaluations.

The compact [machine-readable results](shared_grid_measurements.json) are
checked in. Full diagnostics and outputs are written locally to
`output/shared_grid/legacy.json`, `shared.json` and their MAT counterparts.

A separate cheap 2×2 channel-cell bank (two rational channels, a constant
and zero) took 0.141 s with independent preparation and 0.970 s shared
(medians of three warm runs), despite reducing node evaluations from
495,672 to 65,538. Scalar rational callbacks are vectorized and inexpensive;
whole-page caching and matrix coordination dominate this case. Thus shared
preparation does **not** meet the requirement of being no slower on all
existing kinds of model. `SharedGrid=false` remains the default.

Reproduce with the following commands. For a historical comparison, put
this benchmark script under `tools/` in an isolated export of 9b70d02 and
run its `legacy` mode there first, using the same absolute output folder.
The script does not modify global RNG state, and records actual distinct
callback nodes independently of the library's work counters. `repetitions`
defaults to three; the release table uses one measured run per case.
A fourth argument selects cases for resuming a run, for example
`benchmark_shared_grid('shared','output/shared_grid',1,{'dense300'})`.
The optional `cheapChannels` case reproduces the default-decision fixture;
run it with three repetitions for the quoted median comparison. A shared run started before the baseline finishes may record a
missing difference; subtract the saved output arrays once both finish.

<!-- no-test -->
```matlab
addpath('tools')
benchmark_shared_grid('legacy','output/shared_grid',1);
benchmark_shared_grid('shared','output/shared_grid',1);
```

## Limits and compatibility

- `MaxEvaluations` is a hard shared budget on actual matrix-node calls,
  including reevaluations after eviction. `MaxPoints` still bounds each
  probe. The common matrix grid needs storage proportional to the number
  of channels: `MaxMemoryMB` reserves storage for snapshots/diagnostics,
  the cache and matrix inversion workspace. Unresolved preparation returns
  no kernel, including when a common grid exceeds a point/memory budget.
- `SharedGrid=true` with quadrature is rejected explicitly; it is not
  silently ignored. `UseParallel` is not introduced in 0.9.0. Parallel
  callback evaluation and shared adaptive quadrature remain future work.
- Saved scalar and matrix data remain schema 1. The regression fixture
  `tests/fixtures/kernel_v080.mat` was created with unmodified 0.8.0 code
  and is loaded and simulated by the 0.9.0 tests. New shared snapshots are
  also saved/loaded and remain independent of mutable model callbacks.
- Reuse uses the same channel convolution, per-output sum of channel error
  estimates, and roundoff allowance. Cancellation does not weaken output
  tolerance checks. `SharedGrid` is a preparation option, not a reuse option.

Tests compare independent and shared FFT/de Hoog banks for explicit scalar,
rectangular structured, opaque, batched, channel-cell, unstable, delayed,
direct-term and optional LTI inputs, FOH/ZOH, serialization, exact budgets,
domain checks and symmetry failures. The old scalar engines and their
existing regression tests remain in the validation suite.

## Validation run

On MATLAB R2026a, the release run passed all **153 unit/regression tests**
and all **3 documentation tests**, with no skips. After a final loop-variable
scope cleanup, all seven shared-grid tests passed again and MATLAB's static
check was clean. The example on this page was also executed separately.
The full study/manual/package build is recorded when it finishes.
