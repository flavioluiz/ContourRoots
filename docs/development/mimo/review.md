# Review of the MIMO proposal

Status: review of the planning documents, 2026-09-26, at release 0.5.0
(`e29c5be`). Reviewed: [implementation plan](../mimo_implementation_plan.md),
[API and contracts](api_and_contracts.md),
[numerical methods](numerical_methods.md),
[validation and acceptance](validation_and_acceptance.md).
Reviewed by Claude Code; the plan was drafted by Codex. Nothing here is
implemented.

## Verdict

The proposal is mathematically sound and unusually careful. It keeps the
three objects apart: characteristic values of $H$, transfer poles and
transmission zeros. It warns against the determinant shortcut, keeps the
scalar API unchanged, and makes the incomplete cases explicit. The
fixtures and references checked below are correct.

It is approved **as a plan** with four changes of emphasis (Section 3):

1. release a smaller first stage;
2. show explicitly what the matrix engine adds over the scalar
   determinant;
3. define the batched evaluation layout early;
4. consider MIMO zero-state responses, which are now cheap thanks to
   `ckernel`.

## 1. Independent checks performed

The fixtures were recomputed with MATLAB tools independent of the
proposed code. `pole`/`tzero` were applied to minimal realizations, and
`smithForm` and `det` were computed symbolically.

| Fixture | Plan statement | Independent result |
|---|---|---|
| F03 | `ones(2)/(s+1)`: pole −1 multiplicity 1; `eye(2)/(s+1)`: multiplicity 2 | minimal orders 1 and 2; no zeros ✓ |
| F04 | pole order 2, total multiplicity 3 at −1; residue rank insufficient | minimal order 3; residue `diag(0,1)` has rank 1 ✓ |
| F05 | poles −2 (×2), transmission zero 0, channel zero −1 not a transmission zero | Smith form of numerator `diag(1,s)`; `det G = s/(s+2)^2`; `tzero` = 0 ✓ |
| F06 | pole and zero both at −1, pole at −2; determinant hides the pair | `tzero` = −1, poles {−2, −1}, `det G = 1/(s+2)` ✓ |
| F07 | tall/wide: pole −1, zero 0 | ✓ both orientations |
| F09 | rank 1; pole −1, zero 2 | ✓ |
| F10 | modes −1, −2; transfer `1/(s+1)`; invariant zero −2 | `det P = s+2`; `tzero` of the nonminimal realization = −2; minimal transfer `1/(s+1)` ✓ |
| F11 | `[s+1 1;0 s+1]` geometric multiplicity 1 | `rank H(-1) = 1` ✓ |
| F12 | four zeros of `sin` in `[-0.4,10]×[-1,1]` | `croots`: 0, π, 2π, 3π ✓ |
| F13 | roots `-1+W_k(-e/2)` | `croots` agrees with `lambertw` to 1.5e-12; **see correction below** |

The algebraic statements also check out:
- `ord det G = m_z − m_p` for square full-rank matrices;
- `det P = det H · det G`, with the sign convention `[H −B; C D]`;
- the simple-pole residue `C v w* B / (w* H' v)`;
- the Beyn moments;
- the recovery of local orders from determinantal divisors,
  `α_k = d_k − d_{k−1} − h`;
- on an analytic region, any rank drop of `G` is a zero of every fixed
  `r×r` minor.

The four cited references resolve to the stated works:
- Beyn, arXiv:1003.1580;
- Brennan–Embree–Gugercin, arXiv:2012.14979;
- Güttel–Tisseur, Acta Numerica 2017;
- `nla-group/nep`.

The proposed public names (`cmodes`, `ctzeros`, `cmimo`, `cdyn`) collide
with nothing on the MATLAB path, and the executable documentation tests
skip `docs/development/`, so the non-runnable proposal snippets are safe.

## 2. Correction applied

**F13 branch set.** In the rectangle `[-6 1 -25 25]`, `croots` finds 8
roots. They come from the branches `k = −4, …, 3` of MATLAB's
`lambertw(k, −e/2)`; `k = 0` and `k = −1` give the rightmost conjugate
pair, `−1.1027 ± 1.5026i`. A reference built from the symmetric range
`k = −3…3` misses one root and would make a correct solver look
spurious. The fixture now states its branch set explicitly.

## 3. Recommendations

**3.1 Stage the delivery around `cmodes`.** P0 to P7 is a research
programme; P6 (local Smith–McMillan structure) is a numerical research
problem in its own right. The two motivating applications of this toolbox
are delay systems and the Theodorsen flutter model, and both are matrix
characteristic problems, `det(sI − A − B e^{−sT})` and `H(p,U)`. Suggested
releases:

- **first**: P1 and P2 only (`cdyn`, `cmodes`, channel operations), with the
  aeroelastic section (F15) as the physical acceptance case;
- **then**: P3 (transfer poles with simple multiplicities) and square
  transmission zeros (P4);
- **later**: P5 and P6, with the restricted-scope statements the plan
  already requires.

**3.2 Compare against the determinant route.** Today a user can already
write `croots(@(s) det(H(s)), region, 'AssumeAnalytic', true)` for a
small analytic $H$. The Tutorial 9 aeroelastic example does exactly this.
The plan should state, and test, what `cmodes` adds:
- counting by LU log-determinant without overflow;
- mode shapes (left and right vectors), residuals and conditioning;
- geometric multiplicity (F11);
- robustness when `det` is badly scaled.

Add one fixture where the scalar determinant fails or loses accuracy,
such as a larger or strongly scaled matrix. Also add one where both
routes must agree to about 1e-10. The future tutorial should *start*
from the determinant route readers already know, then show its limits.

**3.3 Define batched evaluation in P1.** Existing scalar handles are
vectorized over `s`. The proposed matrix evaluator takes one scalar node
per call. Contour counting needs thousands of nodes, and the Theodorsen
function calls `besselk` at each one. Per-call overhead will dominate for
small $n$. The plan defers batching; its layout (for example an
`n×n×N` page array) should still be fixed in the P1 contract so that it
can be added later without changing the API.

**3.4 MIMO zero-state responses are now cheap.** The plan excludes MIMO
time simulation. The zero-state response is, however, just superposition,
$y_i = \sum_j G_{ij} * u_j$. With `ckernel` (0.5.0), each channel's kernels
are prepared once. `clsim(M, U, t)` for a `cmimo`/`cdyn` model could be
a small, independent item, with no spectral theory involved. For
`cdyn`, `G_ij(s)` costs one linear solve per frequency node, shared by
all channels of one input. Suggest listing it as a separate, earlier
candidate rather than as out of scope.

**3.5 Minor points.**
- `cpoles(M)` means transfer poles, while `czeros(M)` without channel
  indices is an error. This asymmetry is deliberate: channel zeros and
  transmission zeros must not be confused. The reference pages should
  explain it in one sentence.
- The location target `1e-7` is looser than what the scalar engine
  delivers (typically 1e-12 on the F13 fixture). Keep `1e-7` as the
  acceptance floor for ill-conditioned cases, and add the tighter
  agreement of Recommendation 3.2 for well-conditioned fixtures.
- `ROADMAP.md` still lists "general MIMO nonlinear eigenvalue problems"
  as out of scope. P0 rightly defers rewording it until implementation
  starts. A pointer to this proposal has been added to the roadmap
  without changing that statement.
- The planning documents are dense rule lists, which suits developers.
  When the feature is implemented, the user documentation must follow the
  tutorial style of the rest of the toolbox, as Section 7 of the plan
  already requires. Begin with worked example F05 and the three
  spectral objects, and show a figure.
