%% Response to an arbitrary input (the analogue of lsim)
% Diffusion through a slab, G(s) = exp(-0.7 sqrt(s)): an infinite-
% dimensional system whose only singularity is the branch cut on the
% negative real axis, so 'SingularityBound',0 is justified.
% Run setup_contourroots once per session before this script.

G = @(s) exp(-0.7*sqrt(s));
t = (0:0.02:6).';
u = sin(3*t).*exp(-0.4*t);            % input samples

% 'foh': the input is linear between samples ('zoh' holds each sample).
[y, t, info] = clsim(G, u, t, 'SingularityBound', 0, 'Interpolation', 'foh');
assert(info.converged, info.stopReason);

figure
plot(t, u, '--', t, y, 'LineWidth', 1.5), grid on
legend('input', 'output'), xlabel('t')
title('Diffusion: delayed and smoothed response')
