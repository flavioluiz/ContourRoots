% ContourRoots: poles, zeros and roots of nonrational scalar functions.
% Version 0.8.0 27-Sep-2026
%
% Main functions
%   croots               - Roots of a scalar analytic function in a rectangle.
%   cpoles               - Poles of a transfer function (after cancellations).
%   czeros               - Zeros of a transfer function (after cancellations).
%   cpzmap               - Pole-zero map in a rectangle of the complex plane.
%   ndpair               - Transfer function given as numerator/denominator.
%
% Explicit matrix models and characteristic modes
%   cmodes               - Characteristic values of an analytic square H(s).
%   cmimo                - Matrix transfer from channels, constants or a handle.
%   cdyn                 - Structured model H*x=B*u, y=C*x+D*u.
%   ceval                - Matrix evaluation with shared block solves.
%   cchannel             - Scalar evaluator for one matrix channel.
%
% Zero-state time responses (SISO and explicit MIMO models)
%   cimpulse             - Ordinary impulse response and singular-term metadata.
%   cstep                - Unit-step response by inverse Laplace transform.
%   clsim                - Response to a sampled or function-handle input.
%   ckernel              - Prepare reusable step/ramp kernels for clsim.
%   cinvlaplace          - Inverse Laplace transform (FFT, de Hoog or quadrature).
%
% Solver and compatibility names
%   complex_spectrum     - Scalar spectral solver (cmodes has a separate engine).
%   characteristic_roots - Same as croots (name kept for compatibility).
%   transfer_poles       - Same as cpoles (name kept for compatibility).
%
% Time-delay systems  D(s) + N(s) exp(-s T) = 0  (folder delay/)
%   critical_delays      - Imaginary-axis crossings, delays and directions.
%   unstable_root_count  - Number of roots with Re(s) >= 0 (bound-based count).
%   rhp_root_bound       - Radius containing every right-half-plane root.
%   delay_roots          - Grid/Newton root search for one delay.
%   delay_root_count     - Argument-principle count in a rectangle.
%   pade_delay           - Diagonal [n/n] Padé approximation of exp(-sT).
%   pade_characteristic  - Characteristic polynomial with Padé.
%   pade_critical_delays - Critical delays predicted by a Padé model.
%   match_pade_roots     - Match exact roots to Padé roots.
%
% Utilities
%   contourroots         - Overview and version (type "contourroots").
%   contourroots_version - Installed version.
%
% Documentation: docs/index.md and docs/ContourRoots_manual.pdf in the
% installation folder.
