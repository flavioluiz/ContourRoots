# Performance audit and release-cost study

Baseline: clean `4434195` (ContourRoots 0.8.0), MATLAB R2023b, macOS arm64.
This study distinguishes archive creation from the complete release pipeline.
No publication, installation, version bump or validation-gate removal is
part of the optimization. The original sources were frozen in a separate
scratch copy; the optimized core was measured in a fresh MATLAB process.

## Why a release is not just packaging

`buildtool release` has these dependencies, unchanged by this work:

```text
release
  test       unit + regression tests
  docs       executable Markdown + standalone quickstarts
  manual
    examples all full scientific studies
    figures  regenerate documentation assets
    LaTeX    rebuild the technical manual
  archives   stage -> ZIP -> MATLAB installer
```

All of these tasks execute on each full release. The executable documentation
and standalone quickstarts deliberately cover overlapping examples in
different contexts. The studies include parameter sweeps, independent
oracles and multiple inverse transforms. This is useful release validation,
not compression overhead. A release taking more than ten minutes cannot be
attributed to the archive step without separating these stages.

The isolated baseline packaging measurement (one warm MATLAB process,
fresh staging directory) was:

| Stage | Seconds |
|---|---:|
| Copy/filter released files | 0.158 |
| ZIP | 0.841 |
| ToolboxOptions construction | 2.957 |
| Set metadata/path | 0.027 |
| packageToolbox | 1.456 |
| Total measured stages | 5.439 |

The installer stage above was profiled, so its timing includes profiler
overhead. MATLAB startup is excluded. The installed R2023b implementation
of ToolboxOptions/packageFromOptions was inspected: path derivation is
directory enumeration, and the options-based path packages the declared
files. There was no observed ten-minute dependency scan or network delay
in this isolated run. That does not rule out filesystem, antivirus, startup
or environment delays on another machine.

## Changes that preserve numerical work and qualification

1. **Nested contour reuse.** `spectrum_count` and the matrix counter retain
   old values/logarithms on exact nested grids. Only new nodes are evaluated.
   Coordinates are checked for exact equality before reuse. Integer-count
   agreement, phase/magnitude limits, derivatives, Newton tolerances,
   boundary handling and parent/child conservation are unchanged. At the
   minimum three grids, 48+96+192 node evaluations become 192 unique nodes:
   42.9% fewer contour evaluations, without fewer contour samples.
2. **LU permutation representation.** `cmodes` uses the LU permutation
   vector, O(n) cycle parity rather than O(n²) inversion counting, and row
   selection rather than multiplying a dense permutation matrix by H'.
3. **MIMO cache lookup.** Replace per-node hexadecimal conversion and map
   dispatch with batched, bit-exact uint64 real/imaginary keys and bounded
   preallocated numeric pages. No rounded frequency keys, no persistent
   cache and no changes to cache capacity/eviction policy. Kernel snapshots
   still contain numerical kernels only, never evaluators or caches.
4. **Hybrid helpers.** Evaluate each reciprocal condition estimate once,
   retaining the same pivot threshold and unavailable-chart behavior.
5. **Build visibility and explicit packaging.** Each stage is timed. The
   full `release` path keeps its dependencies; a distinct `package` task
   exposes the existing packaging-only workflow without claiming validation.
   Both archives are built in a temporary staging area before destination
   replacement; a failed installer no longer first deletes a good ZIP.

The profiled baseline hybrid-kernel preparation made 286,786 calls to the
per-node key function, consuming 5.95 s of a profiled 14.87 s preparation.
Inclusive profiler times overlap and must not be added together. Warm,
unprofiled measurements are used for before/after speed comparisons.

## Verification criteria

- Same roots, counts, completeness decisions and scaled residual quality on
  independent polynomial, Lambert-W, distributed and aeroelastic fixtures.
- Same impulse/step/arbitrary-input behavior, singular terms, channel layout,
  error estimates, convergence requirements and reusable-kernel semantics.
- Cache eviction under small budgets must reproduce the generous-cache
  result; callbacks and evaluation counters agree, with no cross-call reuse.
- Both permutation parities and the trace-integral check must pass.
- A failed build stage records its duration and rethrows the original error.
  Unwritable timing logs must not prevent the operation or mask its failure.
- Full release validation dependencies stay present; package-only is explicit.
- ZIP/MLTBX payload, metadata and installation paths are checked without
  installing or publishing the temporary test artifacts.

## Reproduction

Run each baseline/optimized core benchmark in an otherwise idle MATLAB
session. Do not change the source while that benchmark is executing.

```matlab
addpath(fullfile(repo,'tools'))
benchmark_costs(repo,fullfile(scratch,'core'),'core');
benchmark_costs(repo,fullfile(scratch,'pipeline'),'pipeline');
benchmark_costs(repo,fullfile(scratch,'package'),'package');
compare_performance_runs(beforeCoreDir,afterCoreDir);
build_release(fullfile(scratch,'production-package'));
check_release_payload(fullfile(scratch,'production-package'));
```

Use distinct output directories. `core` performs a warm-up and three timed
repeats, plus a separately profiled kernel preparation. `pipeline` executes
tests/docs/studies and writes checkpoints after each stage; it does not
compile the PDF or publish anything. `package` reproduces the original
archive-stage sequence in a temporary output folder; it is separate from
the instrumented production packager `build_release(outputDir)`.

Do not replace accuracy tests with wall-clock thresholds. Timings vary with
warm-up, operating-system caching, thermal state and competing work.
Results from the measured runs and their limitations are summarized below.

## Measured results (2026-09-27)

### Numerical core

Median of three warm, unprofiled runs per implementation:

| Workload | Original [s] | Optimized [s] | Speedup |
|---|---:|---:|---:|
| `cmodes`, Goland at 150 m/s, four pieces | 3.364 | 1.794 | 1.88x |
| `croots(wing_delta)`, same spectral window | 2.644 | 1.469 | 1.80x |
| `ckernel`, hybrid wing, 2x2 force/twist channels | 7.946 | 3.217 | 2.47x |

The original and optimized roots matched with **zero measured discrepancy**.
Simulating the saved original/optimized kernels with the same two inputs
also gave zero discrepancy; resolved masks, error estimates, preparation
evaluation counts and cache-hit counts were identical. This is an A/B
regression result for these fixtures, not a proof of global numerical accuracy.
`compare_performance_runs` automates these checks; unlike wall time, they are
asserted. Existing independent-reference tests remain necessary.

After the cache change, the separately profiled kernel took 6.69 s versus
14.87 s originally. Matrix evaluation now accounts for about 6.08 s of that
6.69 s; the former per-node key-dispatch hotspot is gone. Do not compare
profiled times directly with the unprofiled medians above.

### Full-pipeline baseline: measured lower bound, not a completed release

| Completed stage | Original [s] |
|---|---:|
| Unit/regression tests (138 tests) | 131.79 |
| Executable Markdown and quickstarts (3 tests) | 329.58 |
| Delay study | 166.51 |
| Distributed-system study | 6.23 |
| Coupled-beam study | 28.04 |
| Aeroelastic section study | 209.87 |
| Time-response study | 243.10 |
| **Subtotal of completed stages** | **1115.12 (18 min 35 s)** |

The baseline run was deliberately stopped during the subsequent continuous
wing study, after more than 20 minutes of process wall time. Hybrid/matrix
mode studies, documentation figures and LaTeX compilation are **not** in this
subtotal. Therefore this is a lower bound, not a claimed complete-release
duration. The pipeline ran after the core/profile study in the same original
MATLAB process. Short additional test/package checks overlapped parts of
this diagnostic run, so these stage times are indicative costs, not a
controlled speedup comparison.

The executable documentation alone took 179 s for Markdown snippets and
149 s for standalone quickstarts. Time-response studies also include an
independent adaptive Bromwich inversion at selected times; they are not
merely drawing plots of previously saved results. Repeating these checks
costs real numerical work.

Warnings observed during the original long-running session included a full
JVM code cache (compiler disabled), graphics fallback and slow vector export.
These are possible secondary contributors to variability, **not isolated
causal measurements**. No JVM, graphics quality or solver settings were
weakened to obtain the core speedups.

### Archive-only validation

The instrumented production packager completed in **5.455 s** initially and
**4.942 s** on the final payload check, including staging, both archives and
final copies. This is comparable to the original
5.439 s; there is no demonstrated archive-compression speedup. The benefit
here is visibility, explicit workflow separation and safer replacement.

Both temporary archives were unpacked and their **201 payload files**
compared byte-for-byte against the working tree. UUID, toolbox name/version
and the two MATLAB installation paths were retained. The installer was not
installed; this check does not validate signatures or interactive installation.
Existing `dist/` artifacts were not replaced.

### Optimized qualification

All **149 distinct tests** passed: 146 unit/regression tests (144 in the
pipeline plus the two additional opaque-cache tests), and all three
documentation tests. The logging-failure fixture was also rerun after its
final additions. No optional-toolbox fixture was skipped in this environment.
The eight full studies below completed with all their existing assertions.

| Stage | Optimized [s] |
|---|---:|
| Unit/regression pipeline (144 tests) | 105.14 |
| Executable documentation (3 tests) | 194.71 |
| Delay study | 168.48 |
| Distributed-system study | 4.64 |
| Coupled-beam study | 21.13 |
| Aeroelastic section study | 83.66 |
| Time-response study | 167.92 |
| Continuous-wing study | 98.29 |
| Hybrid/implicit assembly comparison | 147.40 |
| Matrix-mode/determinant comparison | 11.90 |
| **Total measured pipeline stages** | **1003.27 (16 min 43 s)** |

This is **not** a complete `buildtool release` timing: PDF/figure regeneration
and archive creation are excluded from the table. Additional targeted tests
and packaging checks briefly overlapped this qualification run. Consequently,
do not use the pipeline table to claim a controlled whole-release speedup;
the controlled numerical-core comparison is given separately above.

Independent checks included the continuous-wing flutter point
U = 136.983977449 m/s and omega = 70.033012953 rad/s, finite-element convergence,
vacuum frequencies and time-domain inversion cross-checks. Hybrid versus
implicit transfer matrices agreed within 1.67e-13 relative error over the
sampled dry/Goland/tapered models. The full `cmodes` study retained its
subdivision/closed-form comparisons and detected the same determinant
underflow/overflow limitations at scales 1e-200 and 1e200.

Both normal and deliberately unwritable timing-log destinations were tested.
The latter produced warnings without suppressing an operation or changing
its original exception; the packager still produced valid, byte-checked
archives when its optional timing-report path was unavailable.

Static analysis found no errors in the changed numerical routines; remaining
messages are pre-existing compact-style/shared-loop-variable warnings and
the compatibility notice for `verLessThan`. `git diff --check` passed.
Qualification was performed on R2023b/macOS arm64, not every supported MATLAB
release or operating system. The existing PDF was packaged as-is; it was
not regenerated, installed or published during this audit.

## Remaining opportunities and limits

- Full studies intentionally redo parameter sweeps and independent inversion
  oracles. Incremental result reuse would require source/options/environment
  fingerprints and tests against stale results. It was not introduced here.
- The delay study did not improve materially (166.5 s originally, 168.5 s in
  the optimized run). It repeatedly calls the separate legacy `delay_roots`
  routine and exports vector figures. Its Newton/Padé sweeps were not
  rewritten: the speedups above apply to the identified core hot paths,
  not to every computation performed by the release.
- Documentation examples overlap scientific regressions, but execute in
  different contexts. Removing that coverage would be a behavioral change
  to release qualification, not a safe compression optimization.
- Cache lookup is now batched; it is not an asymptotically constant-time hash
  table. Very large node banks or unusual tiny-batch workloads deserve
  separate profiling. Allocations stay within the existing estimated cache
  allowance; MATLAB temporaries and callback allocations are not hard-capped.
- Further large-scale gains may come from sparse factorizations, analytic
  model batching or parallel parameter sweeps. Those require numerical,
  memory and reproducibility studies beyond the measured hot paths here.
- No manual/PDF content, numerical tolerances, contour refinement criteria,
  independent references or release-validation dependencies were removed.
