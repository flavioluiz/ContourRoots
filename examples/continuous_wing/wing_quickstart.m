%% Flutter of a continuous wing, without discretizing the structure or the air
% The Goland wing: a uniform cantilever that bends and twists, in an
% incompressible flow. The structure is NOT discretized (no finite
% elements, no modal truncation) and Theodorsen's aerodynamics is NOT
% approximated (no rational fit): the modes are the zeros of an exact
% characteristic function Delta(s,U). Run setup_contourroots once per
% session before this script. No additional toolbox is needed.
repo = fileparts(fileparts(which('croots')));
addpath(fullfile(repo, 'examples', 'models'))

%% 1. The characteristic function of the continuous wing
% Inside the wing, bending w(y,s) and twist alpha(y,s) obey ODEs in the span
% coordinate y whose solution is an exact matrix exponential; clamping the
% root and freeing the tip gives a square matrix K(s,U). Its determinant
% plays the role of det(sI - A) for this infinite-dimensional system.
wing  = wing_model('goland');
Delta = @(s, U) wing_delta(s, U, wing);

%% 2. First check: in vacuum the beam frequencies are known exactly
dry = wing_model('dry');                          % no air, no bending-torsion coupling
r = croots(@(s) wing_delta(s, 0, dry), [-2 2 1 320], 'AssumeAnalytic', true);
e = dry.strips;  L = dry.L;
bending = [1.875104068711961 4.694091132974174].^2*sqrt(e.EI/e.mu)/L^2;
torsion = [1 3]*pi/(2*L)*sqrt(e.GJ/e.Ialpha);
exact = sort([bending torsion]).';
table(exact, sort(imag(r)), 'VariableNames', {'closed_form', 'croots'})

%% 3. Is the wing stable?  Count the roots with Re(s) > 0
% Theodorsen's function has a branch cut on the NEGATIVE real axis, so a
% rectangle in the right half-plane is a valid search region.
rhp = [1e-3 60 -400 400];
for U = [120 150]
    [p, info] = croots(@(s) Delta(s, U), rhp, 'AssumeAnalytic', true);
    fprintf('U = %3d m/s: %d unstable roots (%s)\n', U, numel(p), info.status);
end

%% 4. The flutter speed: where the least-damped mode crosses Re(s) = 0
window = [-60 35 1 380];                          % above the branch cut
alpha  = @(U) max(real(croots(@(s) Delta(s, U), window, 'AssumeAnalytic', true)));
Uf = fzero(alpha, [120 150]);
p  = croots(@(s) Delta(s, Uf), window, 'AssumeAnalytic', true);
[~, k] = max(real(p));
fprintf('Flutter: U = %.4f m/s, omega = %.4f rad/s\n', Uf, imag(p(k)));

%% 5. Time response to a tip force, below and above flutter
% Tip deflection in mm per kN of tip force, from the same continuous model:
% no modal ODE, no aerodynamic lag states. The inversion line must lie
% right of every singularity. SingularityBound = 5 is an ASSUMPTION,
% supported (not proved) by step 3: no unstable root at 120 m/s and one
% pair with Re(s) = 3.70 at 150 m/s, inside a finite rectangle.
% Tolerances are in the units of the output: 1e-3 mm (1 micrometre).
tip = @(U) @(s) 1e6*wing_transfer(s, U, wing, 1, 1);   % mm per kN
t = (0:0.005:0.5).';
opts = {'SingularityBound', 5, 'AbsTol', 1e-3, 'RelTol', 1e-3};
wBelow = cstep(tip(120), t, opts{:});                  % 1 kN step at the tip
wAbove = cstep(tip(150), t, opts{:});
figure, plot(t, wBelow, t, wAbove), grid on
legend('120 m/s (below flutter)', '150 m/s (above flutter)', 'Location', 'northwest')
xlabel('t [s]'), ylabel('tip deflection [mm]')
title('Continuous wing, exact Theodorsen aerodynamics: 1 kN tip step')
