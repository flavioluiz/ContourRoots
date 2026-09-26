%% First step response of a nonrational transfer function
% Loop with a delay: x'(t) = -x(t) - 0.5 x(t-1) + u(t), i.e.
%     G(s) = 1/(s + 1 + 0.5 e^{-s}),
% which has infinitely many poles and no finite state-space model.
% Run setup_contourroots once per session before this script.

G = @(s) 1./(s + 1 + 0.5*exp(-s));
t = (0:0.02:8).';

% 'SingularityBound',0: G has no singularity with Re(s) > 0, because there
% |s+1| >= 1 > 0.5 >= |0.5 e^{-s}|. The inversion line is placed right of it.
[y, t, info] = cstep(G, t, 'SingularityBound', 0);
assert(info.converged, info.stopReason);

figure
plot(t, y, 'LineWidth', 1.5), grid on
yline(1/1.5, ':', 'static gain G(0)')
xlabel('t'), ylabel('step response')
title('Loop with a delay: step response computed from G(s) itself')
