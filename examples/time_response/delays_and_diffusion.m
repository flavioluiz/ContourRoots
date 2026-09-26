%% Two nonrational step responses with exact references
% A transport delay and a diffusion process, compared with their
% analytical step responses.
% Run setup_contourroots once per session before this script.

t = (0:0.02:4).';

% Delay: exp(-0.713 s)/(s+1). The delay time is not on the grid, so no
% sample falls exactly on the jump.
[delayed, ~, id] = cstep(@(s) exp(-0.713*s)./(s + 1), t, 'SingularityBound', 0);

% Diffusion: exp(-0.7 sqrt(s)); its step response is erfc(0.7/(2 sqrt(t))).
[diffusive, ~, ih] = cstep(@(s) exp(-0.7*sqrt(s)), t, 'SingularityBound', 0);
assert(id.converged && ih.converged);

reference = zeros(size(t));
reference(2:end) = erfc(0.7./(2*sqrt(t(2:end))));
figure
plot(t, delayed, t, diffusive, t, reference, 'k:', 'LineWidth', 1.5), grid on
legend('transport delay', 'diffusion', 'diffusion, analytical'), xlabel('t')
