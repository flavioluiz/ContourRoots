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
| `buildtool examples` | runs the full studies in `examples/` (a few minutes); results go to `output/` |
| `buildtool manual` | runs the studies (`examples`), regenerates `docs/assets`, and compiles the PDF manual (needs `latexmk`) |
| `buildtool release` | tests, docs, studies and manual, then builds `dist/ContourRoots.zip` and `dist/ContourRoots.mltbx` |
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
