%% Flutter of a wing section with exact Theodorsen aerodynamics
% Stability of a two-degree-of-freedom (plunge-pitch) wing section as the
% airspeed U grows, WITHOUT a rational (Pade-type, Jones or Roger)
% approximation of the unsteady aerodynamics. Run setup_contourroots once
% per session before this script. No additional toolbox is needed.
repo = fileparts(fileparts(which('croots')));
addpath(fullfile(repo, 'examples', 'models'))

%% 1. The model: a "characteristic polynomial" that is not a polynomial
% Structure M q'' + K q = aerodynamic loads, q = [plunge; pitch]. In the
% Laplace domain the loads depend on Theodorsen's function C(sb/U), built
% from Bessel functions. The modes are the roots of
%     Delta(s,U) = det( s^2 M + K - A(s,U) ),
% the analogue of det(sI - A) for a state-space model.
model = aeroelastic_section('nasa');             % NASA/TP-2015-218765 test case
Delta = @(s, U) aeroelastic_delta(s, U, model);  % analytic except on s <= 0 real

%% 2. Is the section stable at U = 170 ft/s?  Count the roots with Re(s) > 0
% Theodorsen's function has a branch cut on the NEGATIVE real axis, so the
% whole right half-plane is a valid search region.
rhp = [1e-3 60 -250 250];                        % [xmin xmax ymin ymax]
[p, info] = croots(@(s) Delta(s, 170), rhp, 'AssumeAnalytic', true);
fprintf('U = 170 ft/s: %d unstable roots (%s)\n', numel(p), info.status);

%% 3. Where are the modes?  Look above the branch cut
% The coefficients are real, so the roots come in conjugate pairs: the
% upper half-plane is enough.
upper = [-100 50 0.5 180];
modes = croots(@(s) Delta(s, 170), upper, 'AssumeAnalytic', true)

%% 4. Root locus in the airspeed
speeds = 0:20:240;
figure, hold on
for U = speeds
    r = croots(@(s) Delta(s, U), upper, 'AssumeAnalytic', true);
    plot(real(r), imag(r), 'o', 'MarkerFaceColor', [U/250 0.3 1-U/250], ...
        'MarkerEdgeColor', 'none')
end
xline(0, ':'), grid on, xlabel('Re(s) [1/s]'), ylabel('Im(s) [rad/s]')
title('Modes of the NASA section, U = 0 to 240 ft/s (blue to red)')

%% 5. Flutter speed: where the least-damped mode crosses Re(s) = 0
alpha = @(U) max(real(croots(@(s) Delta(s, U), upper, 'AssumeAnalytic', true)));
Uf = fzero(alpha, [150 200]);                    % NASA reports 173.26 ft/s
fprintf('Flutter speed: %.2f ft/s\n', Uf);
