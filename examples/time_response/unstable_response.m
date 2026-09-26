%% Step response of an unstable system
% G(s) = 1/(s - 0.3) has a pole at +0.3. The inversion line must lie to its
% right; for coefficient vectors ContourRoots finds this from the pole.
% (For a function handle you would pass 'SingularityBound',0.3.)
% Run setup_contourroots once per session before this script.

t = (0:0.05:10).';
[y, t, info] = cstep(ndpair(1, [1 -0.3]), t);
assert(info.converged, info.stopReason);

figure
plot(t, y, t, expm1(0.3*t)/0.3, 'k--', 'LineWidth', 1.5), grid on
legend('computed from G(s)', 'analytical (e^{0.3t}-1)/0.3')
xlabel('t'), ylabel('step response')
