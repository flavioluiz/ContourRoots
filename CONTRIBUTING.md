# Contributing to ContourRoots

Thank you for your interest. Bug reports, examples and improvements are
welcome.

## Reporting a problem

Please include the MATLAB release (`version`), the ContourRoots version
(`contourroots_version`), a minimal script that reproduces the problem, and
the full output, including `info.status` and `info.unresolvedBoxes`.

## Development tasks

All tasks run from the repository folder with MATLAB's `buildtool`
(R2022b or later). `make <task>` does the same from a shell.

| Command | What it does |
|---|---|
| `buildtool test` | unit and regression tests (about a minute) |
| `buildtool docs` | runs every `matlab` block of the README and `docs/`, the quick-start scripts, and checks that code shown in tutorials matches the example files |
| `buildtool examples` | runs the full studies in `examples/` (potentially tens of minutes); results go to `output/` |
| `buildtool manual` | runs the studies (`examples`), regenerates `docs/assets`, and compiles the PDF manual (needs `latexmk`) |
| `buildtool release` | tests, docs, studies and manual, then builds `dist/ContourRoots.zip` and `dist/ContourRoots.mltbx` |
| `buildtool package` | packages the existing code/examples/docs/PDF only; does **not** run validation or regenerate results |
| `buildtool clean` | removes `dist/` and the manual build folder |

## Guidelines

- Keep changes to the numerical core (`matlab/complex_spectrum.m`,
  `matlab/private/`) separate from interface or documentation changes, and
  add a regression test for every behavior change.
- Every `matlab` code block in the documentation must run. Mark a block
  that must not run in tests (downloads, `savepath`) with `<!-- no-test -->`
  on the line before the fence.
- A code block preceded by `<!-- file: path -->` must be identical to that
  file.
- Use the words "exact", "complete" and "certified" carefully: see
  [Diagnostics and limits](docs/diagnostics_and_limits.md).
- Documentation is written in English.

## Releases

1. Update `contourroots_version.m`, `Contents.m`, `CITATION.cff` and
   `CHANGELOG.md`.
2. Run `buildtool release`.
3. Tag the commit (`vX.Y.Z`) and attach `dist/ContourRoots.zip` and
   `dist/ContourRoots.mltbx` to the GitHub release. Keep these exact file
   names: the README links to `/releases/latest/download/<name>`.

`release` intentionally retains all validation gates. If an already validated,
unchanged tree only needs its archives rebuilt, use `buildtool package`
(`make package` from a shell). This is **not** an incremental validation mode:
after code or documentation changes, run `release` again before publication.
Do not infer PDF freshness from successful packaging alone.
A full release can take tens of minutes: it recomputes the scientific studies
and executes the documentation, rather than only creating ZIP/MLTBX files.

Build stages print start/end times and write JSON reports below
`output/build_timings/`. The packager separately records its staging, ZIP,
installer and copy times in `dist/packaging_timings.json`. These distinguish
slow numerical studies from slow archive creation. `build_release(outputDir)`
can package into a scratch destination without touching existing releases;
archives are prepared successfully before existing destinations are replaced
(the two final replacements are not a single atomic transaction).

For reproducible performance measurements, see the source checkout's
`docs/development/performance_audit.md` and `tools/benchmark_costs.m`
(development notes/tools are intentionally omitted from installers).
`check_release_payload(outputDir)` verifies archive payload and installer
metadata without installing the add-on. Timing is observational, not a machine-dependent
pass/fail threshold. Correctness, work counts and cache isolation are tested
in `test_performance_invariants`.
