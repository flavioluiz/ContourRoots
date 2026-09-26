# Time-response examples

Short scripts (run `setup_contourroots` first), explained in
[Tutorial 10](../../docs/tutorials/10_time_response.md):

| Script | Shows |
|---|---|
| `first_step.m` | `cstep` on a loop with a delay, and how to justify `SingularityBound` |
| `arbitrary_input.m` | `clsim` with a sampled input through a diffusion model |
| `delays_and_diffusion.m` | step responses of a delay and of diffusion against exact formulas |
| `unstable_response.m` | an unstable system: the inversion line right of the pole |

`run_time_response_study.m` is the full validation study (about a minute):
delayed feedback against a method-of-steps series, diffusion against
`erfc`, a heated rod against separation of variables, the coupled beam of
Tutorial 8 against finite elements, and the aeroelastic pitch response of
Tutorial 9 computed by two independent inversion methods. It writes CSV,
MAT and PNG files to `output/time_response/`. All responses start from rest.
