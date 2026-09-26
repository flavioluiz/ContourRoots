# ContourRoots: Repository Organization and Publication Plan

Status: proposed plan, not an implemented reorganization.

This document describes how to turn the current research repository into a
public, reproducible MATLAB toolbox. Existing function names and numerical
behavior should be preserved during the initial migration. Changes to the
algorithms should be reviewed and tested separately from file moves.

## 1. Purpose and positioning

**ContourRoots — A MATLAB Toolbox for Nonrational Poles and Zeros**

Suggested short description:

> Find complex roots of scalar characteristic equations and poles of scalar
> nonrational transfer functions in bounded regions, with contour-based
> numerical checks and worked examples from time-delay and distributed systems.

The first release should emphasize:

- A small interface: `characteristic_roots`, `transfer_poles`, and
  `complex_spectrum`.
- Symbolic expressions, numerical function handles, and explicit analytic
  numerator/denominator pairs, with documented support for SISO LTI models.
- A distinction between internal characteristic roots and observable
  input/output poles after numerical cancellation checks.
- Diagnostics, reproducible examples, and clearly stated numerical limits.
- Optional time-delay utilities, including Padé comparisons and critical-delay
  calculations for their supported model classes.

Do not claim a new root-finding principle, universal completeness, certified
root locations, or superiority over other solvers without supporting evidence.
Here, evaluating the original model without Padé or modal truncation does not
mean that the computed roots are exact symbolic values.

## 2. Proposed repository layout

```text
ContourRoots/
├── README.md                       # English landing page and quick start
├── LICENSE                         # Maintainer-approved code license
├── CITATION.cff                     # Verified authorship and release metadata
├── CHANGELOG.md
├── CONTRIBUTING.md
├── THIRD_PARTY_NOTICES.md
├── CONTOURROOTS_REPOSITORY_PLAN.md
├── setup_contourroots.m             # Add only the public MATLAB folder
├── Contents.m                      # MATLAB help entry point
├── Makefile                        # Optional developer/report commands
├── matlab/
│   ├── characteristic_roots.m
│   ├── transfer_poles.m
│   ├── complex_spectrum.m
│   ├── delay_roots.m               # Existing delay API, retained initially
│   ├── ...                        # Supported delay and Padé utilities
│   └── private/
│       ├── spectrum_model.m
│       ├── spectrum_count.m
│       ├── spectrum_newton.m
│       ├── spectrum_solve.m
│       └── spectrum_values.m
├── examples/
│   ├── quickstart/
│   ├── cancellations/
│   ├── time_delay/
│   ├── distributed_systems/
│   ├── coupled_beam/
│   └── models/                    # Reusable physical model helpers
├── tests/
│   ├── run_tests.m
│   ├── unit/
│   ├── regression/
│   └── examples/                  # Smoke tests of documented examples
├── benchmarks/
│   ├── README.md
│   ├── cases/
│   └── adapters/                  # Optional external-solver interfaces
├── docs/
│   ├── index.md
│   ├── getting_started.md
│   ├── model_inputs.md
│   ├── numerical_method.md
│   ├── diagnostics_and_limits.md
│   ├── troubleshooting.md
│   ├── api/
│   ├── tutorials/
│   ├── validation.md
│   ├── related_software.md
│   ├── assets/                    # Selected, reproducible published figures
│   └── reports/                   # Curated PDF manuals
├── report/                        # LaTeX/TikZ sources, not runtime code
├── references/
│   ├── README.md                  # Bibliography, DOI and source links
│   └── references.bib
├── tools/                         # Documentation/release build helpers
├── output/                        # Local generated results; normally ignored
└── .github/
    ├── workflows/
    └── ISSUE_TEMPLATE/
```

Keep public `.m` functions together in `matlab/` for the first release. This
preserves unqualified calls and MATLAB's `private/` lookup rules. Organize the
API conceptually through documentation rather than introducing a breaking
`+contourroots` namespace now.

`setup_contourroots` should locate itself with `mfilename('fullpath')` and add
only `matlab/` to the current session's path. It should not call `savepath`,
modify startup files, or recursively add tests and report folders. Examples
and test runners should add their own helper paths temporarily and restore
them on exit.

## 3. Migration of existing files

| Current files | Proposed destination and treatment |
|---|---|
| `complex_spectrum.m`, `characteristic_roots.m`, `transfer_poles.m` | `matlab/`; primary public API, unchanged names and signatures |
| `private/spectrum_*.m` | `matlab/private/`; move with their callers |
| `delay_roots.m`, `delay_root_count.m` | `matlab/`; document as a separate, existing delay-specific workflow |
| `critical_delays.m`, `rhp_root_bound.m`, `unstable_root_count.m` | `matlab/`; document exact model assumptions and exceptional cases |
| `pade_delay.m`, `pade_characteristic.m`, `pade_critical_delays.m`, `match_pade_roots.m` | `matlab/`; optional comparison utilities, not prerequisites for the general solver |
| `distributed_model.m`, `coupled_beam_model.m`, `coupled_beam_fem.m` | `examples/models/`; retain names, document the path change and example-helper status |
| `run_delay_study.m`, `run_pade_critical_study.m`, `run_pade_pitfalls.m` | `examples/time_delay/`; preserve reproducibility of the original studies |
| `run_nonrational_examples.m` | `examples/distributed_systems/` |
| `run_coupled_beam_study.m` | `examples/coupled_beam/` |
| `tests/run_tests.m`, `tests/test_complex_spectrum.m`, `tests/test_coupled_beam.m` | Retain a single test entry point; split the suite only after a passing relocation baseline |
| `docs/nonrational_toolbox.md`, `docs/coupled_beam.md` | English tutorial chapters; retain old paths as short navigation pages or clearly labeled Portuguese editions |
| `report/*.tex`, `report/*.tikz` | Retain as report sources; update input paths and consolidate bibliography references |
| `report/report_original.tex`, `report/PLANO_ESTRUTURA.md` | Historical material; keep outside the main user navigation and package |
| `output/` figures, tables, CSVs and MAT files | Separate curated documentation assets from disposable generated output |

Before moving anything, record the current commit, MATLAB/toolbox versions,
test results, and representative numerical outputs. After moving files, run
the same tests and studies and compare numerical results using tolerances,
not binary equality of PDFs or MAT files.

Audit all `pwd`, `addpath`, `fullfile`, LaTeX `input`/image paths, and Makefile
targets. Resolve resources relative to source locations, not the caller's
working directory. Verify operation from both the repository root and an
unrelated working directory, including a path containing spaces.

## 4. Main README specification

The English README should answer “what is it, can it solve my problem, and how
do I run it?” before discussing the research history.

Recommended order:

1. Name, one-sentence purpose, and one representative pole-map figure.
2. Supported problems and explicit non-goals.
3. Installation and dependency table with tested MATLAB versions.
4. A copy-and-paste quick start that does not require Symbolic Math Toolbox.
5. An explicit numerator/denominator example for transfer poles.
6. Optional symbolic input and links to LTI-model restrictions.
7. Meaning of the returned diagnostics and numerical completeness.
8. Links to tutorials, API reference, validation and related software.
9. How to run tests and reproduce selected figures.
10. License, citation, contributions and release information.

Suggested first numerical example, after installation:

```matlab
Delta = @(s) 1 + s + s.^2 + (2*s + 3).*exp(-s);
region = [-8 2 -20 20];
[r,info] = characteristic_roots(Delta,region, ...
    'AssumeAnalytic',true,'Plot',true);
disp(r)
disp(info.status)
```

Explain immediately that the analyticity declaration is justified here by
the expression, not inferred from an arbitrary function handle. State that
`info.complete` is a conditional numerical check, not a proof based on
interval arithmetic. Every README example must have an executable test.

Do not advertise arbitrary matrix-valued characteristic equations, general
MIMO models, automatic handling of branch cuts, or complete infinite spectra.
Do not imply that every `tf/ss` delay realization receives the same guarantees
as an explicit analytic characteristic function.

## 5. Complete, didactic documentation

### 5.1. Learning sequence

Organize the manual from a minimal example to physical modeling:

1. **Getting started:** installation, a first root search, plotting and reading
   `info`; choose a rectangle that does not pass through roots.
2. **What is a pole?** Characteristic modes versus transfer poles, observable
   modes, multiplicity, removable singularities and cancellations.
3. **Providing a model:** symbolic input, vectorized handles, analytic N/D
   pairs, polynomial coefficients and supported LTI adapters. Explain optional
   derivatives and initial seeds.
4. **How the algorithm works:** argument principle, winding counts,
   subdivision, contour refinement, Newton polishing and local cancellation
   checks. Include pseudocode and a small worked contour example.
5. **Diagnostics and failure cases:** unresolved cells, close roots, boundary
   roots, near cancellations, overflow, scaling, branch and accumulation
   points, and the exploratory single-handle pole mode.
6. **Time-delay systems:** derive the characteristic equation; explain the
   difference between the general solver and the legacy grid/Newton delay API.
7. **Padé and critical delays:** magnitude/phase conditions, crossing
   directions, stability intervals and limits of finite-order approximations.
   State the supported real-polynomial, single-delay assumptions explicitly.
8. **Distributed systems:** heat, wave, duct and beam examples, with boundary
   conditions, derivations and regularization of removable singularities.
9. **Coupled beam and oscillator:** schematic, PDE/ODE coupling, sign
   conventions, characteristic determinant, transfer functions, direct N/D
   implementation, parameter sweeps and independent finite-element validation.
10. **Validation and comparison:** reproducible cases, error metrics, known
    limitations, related packages and a justified choice of solver.

Each physical tutorial should contain: problem statement, assumptions and
units, schematic when useful, derivation, complete executable MATLAB code,
expected results, physical interpretation, numerical caveats and references.
Distinguish a window-dependent spectral abscissa from a global stability claim.

### 5.2. The direct beam example is a required tutorial

Preserve the complete example currently in section 5.1 of
`docs/coupled_beam.md`, not just a call to `coupled_beam_model`:

- Define the bare beam numerator and denominator explicitly.
- Assemble the coupled `N` and `D` with the mass, stiffness and damping.
- Include the entire local regularization function in the displayed script.
- Call both `transfer_poles` and `characteristic_roots`.
- Explain why the factors are analytic despite intermediate fourth roots.
- Show the expected five conjugate pairs and the decoupled cancellation case.
- Explain normalized parameters and how to use dimensional model helpers.

Create a standalone runnable `.m` example as the canonical code source. Use
marked code regions or a checked extraction process to keep the displayed
Markdown/LaTeX code identical to the executable example.

### 5.3. API reference and terminology

For each public function, document syntax, types, defaults, assumptions,
outputs, examples, errors, limitations and related functions. Give every
`info` field a precise meaning, particularly multiplicity, location radius,
residual, cancellation uncertainty, completeness and certification.

Audit the existing use of “exact”, “all”, “certified” and “stable”. In
particular, distinguish a mathematically valid half-plane root bound from
the reliability of the floating-point contour count used with that bound.
Apply the same terminology to help comments, README, tables and reports.

Use English as the authoritative publication language. Preserve existing
Portuguese material where useful, but label its language and version so that
readers do not mistake an older translation for the current API reference.

## 6. Examples, figures and reproducibility

Provide short examples separately from full parameter studies. A first-use
example should not compile LaTeX or launch every expensive sweep.

Each example should declare its parameters, region, tolerances, dependencies
and expected qualitative or quantitative result. Plot options should be
controllable for headless execution. Avoid `clear all`, persistent path
changes, and modifications to the user's workspace outside the example's
documented outputs.

Keep a small curated collection of figures and final reports in `docs/`.
Write regenerated CSVs, MAT files and temporary renders to `output/` or
`tmp/`. Preserve existing research results until the replacement organization
has been verified; do not delete them as part of a blind cleanup.

Record MATLAB version, optional toolbox availability, source commit,
parameters and solver settings with benchmark/study results. Keep schematics
as editable TikZ or SVG sources. Figure captions must state omitted conjugate
poles, different axis scales and the search window where relevant.

## 7. Tests and independent comparisons

Separate the following layers:

- **Core tests:** known polynomial roots, analytic transcendental examples,
  repeated roots, complex coefficients, conjugate symmetry where applicable,
  validation of inputs, contour boundaries and unresolved searches.
- **Pole tests:** total/partial cancellations, near pole-zero pairs, hidden
  modes and single-handle exploratory behavior.
- **Optional-adapter tests:** symbolic and Control System Toolbox inputs,
  including explicit skips when dependencies are unavailable.
- **Physics regressions:** static beam compliance, dimensional scaling,
  decoupled limits, distributed-model references and FEM convergence.
- **Delay regressions:** zero delay, known critical delays, crossing direction,
  Padé comparisons and applicability checks for half-plane bounds.
- **Documentation tests:** run all quick starts and the direct beam script;
  verify expected counts, residuals and statuses, not just lack of exceptions.

Add public benchmarks against GRPF and a contour-integration root solver;
use QPmR or TDS-CONTROL for compatible delay cases. Treat these as optional
development dependencies, never mandatory runtime dependencies. Review their
licenses before copying or distributing any external implementation.

Use matched domains and meaningful accuracy targets. Report missed and
spurious roots, multiplicities, residuals, reference-location errors,
function evaluations and elapsed time. Publish failures and parameter tuning
alongside successes. A small residual alone is not a location-error bound.
Do not rank methods solely by their default settings or by runtime alone.

Candidate reference projects:

- [GRPF](https://github.com/PioKow/GRPF)
- [Root Finding Cauchy Integration Method](https://www.mathworks.com/matlabcentral/fileexchange/68750-root-finding-cauchy-integration-method)
- [QPmR publication](https://doi.org/10.1007/978-3-319-01695-5_22)
- [TDS-CONTROL](https://twr.cs.kuleuven.be/research/software/delay-control/)
- [cxroots](https://github.com/rparini/cxroots)
- [RootsAndPoles.jl](https://github.com/fgasdia/RootsAndPoles.jl)

## 8. Dependencies, automation and distribution

Determine the minimum MATLAB release by testing, not by guessing from syntax.
Publish separate requirements for the numerical core, symbolic input, LTI
adapters, examples and document generation. Do not claim Octave compatibility
until the supported subset has been tested.

Use MATLAB scripts/functions as the primary entry points; retain Makefile
targets as convenience wrappers so Windows users are not required to install
`make`. Separate study generation from report compilation, allowing an
existing dataset to be used for documentation-only changes.

Proposed CI stages:

1. Run core tests and quick-start examples headlessly.
2. Run available optional-adapter tests and report skipped coverage.
3. Check documentation links and synchronized code snippets.
4. Run slower physical regressions in a separate job.
5. Build and visually inspect PDF reports for releases.

Check MATLAB CI licensing/runner requirements before promising a version
matrix. External solvers should not be downloaded and run implicitly by the
ordinary test command. Package a curated runtime/examples/docs subset for
File Exchange; consider a `.mltbx` installer only after the folder-based
installation is reliable.

## 9. License, references and public-release safety

Select a license with the maintainer's approval. MIT is a candidate for
original code, not an assumed decision. Verify authorship and provenance;
record third-party notices, bibliography and any separately licensed assets.
Do not invent author identifiers, affiliations, release DOIs or citations.

The current repository contains
`references/Curtain_Morris_2009_transfer_functions.pdf`. Local possession does
not establish permission to redistribute it. Before public publication:

- Confirm redistribution rights for the exact PDF version, or replace the
  public copy with bibliographic metadata, DOI and an authorized source link.
- Check repository history as well as the current tree: deleting the current
  file alone does not remove it from previous commits.
- If rights are unclear, prefer a reviewed, sanitized public export while
  preserving the local research repository. Any history rewrite or deletion
  requires explicit maintainer authorization and a recovery plan.

Audit tracked files for credentials, personal paths and nonredistributable
material. Do not push, publish a release, upload to File Exchange or modify
remote visibility without an explicit publication request.

## 10. Implementation phases and acceptance gates

### Phase A — Inventory and baseline

- [ ] Record the source commit, environment and passing test results.
- [ ] Inventory public functions, helpers, datasets and documentation links.
- [ ] Resolve licensing and the reference-PDF distribution decision.

Acceptance: a reproducible baseline and a documented public-release scope.

### Phase B — Layout and installation

- [ ] Move code and private helpers together; add setup and help entry points.
- [ ] Relocate examples and update all source-relative paths.
- [ ] Preserve public solver signatures and document helper-path changes.
- [ ] Pass the existing regression suite before algorithm changes.

Acceptance: clean installation and matching results from different directories.

### Phase C — User documentation

- [ ] Write the English README, input guide, API reference and limitations.
- [ ] Publish the complete executable direct N/D beam example.
- [ ] Translate and organize delay, distributed-system and beam tutorials.
- [ ] Synchronize code blocks and regenerate selected reports and figures.

Acceptance: a new user can reproduce a basic search and the beam example
without reading internal solver code or manually reconstructing missing steps.

### Phase D — Validation and packaging

- [ ] Add documentation smoke tests and optional-dependency coverage.
- [ ] Run at least one independent general-root-solver comparison and one
  compatible delay-solver comparison; record discrepancies and limitations.
- [ ] Prepare release metadata, license, citation and a clean-install archive.
- [ ] Check figures, links, package contents and MATLAB-version statements.

Acceptance: the candidate release is reproducible and contains no unsupported
accuracy, completeness, compatibility or novelty claims.

### Phase E — Maintainer-approved publication

- [ ] Review the final public tree and any history to be published.
- [ ] Choose the initial version and publication destinations.
- [ ] After explicit approval, publish the GitHub repository/release and
  File Exchange entry with consistent names and descriptions.
- [ ] Verify installation from the actual published artifact.

Defer GUI development, a breaking package namespace, MIMO nonlinear
eigenproblems, rigorous interval certification and broad performance claims
to separately scoped future work. The initial release should make the
existing capabilities dependable and understandable.
