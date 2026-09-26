%% Time-delay systems without Padé
% Stability of x'(t) = -x(t) - 2 x(t - T), whose characteristic equation
% is s + 1 + 2 e^{-sT} = 0, written as D(s) + N(s) e^{-sT} = 0 with
% D(s) = s + 1 and N(s) = 2 (coefficient vectors, descending powers).
% Run setup_contourroots once per session before this script.
N = 2;
D = [1 1];

%% 1. Critical delays: where roots cross the imaginary axis
crit = critical_delays(N, D, 10)        % all crossings for 0 <= T <= 10
% The first crossing is at T = 2*pi/(3*sqrt(3)) = 1.2092, frequency sqrt(3).
% CrossingSpeed = Re(ds/dT) > 0: the pair moves into the right half-plane.

%% 2. Number of unstable roots for a given delay
Z1 = unstable_root_count(N, D, 1.0)     % 0: stable
Z2 = unstable_root_count(N, D, 1.5)     % 2: one unstable pair

%% 3. Roots for one delay, with the general solver
T = 1.5;
F = @(s) s + 1 + 2*exp(-s*T);
[r,info] = croots(F, [-6 1 -30 30], 'AssumeAnalytic', true);
r(1:2)                                  % the unstable pair, Re(s) > 0

%% 4. What a first-order Padé model would predict
pade = pade_critical_delays(N, D, 10, 1)
% Padé [1/1] predicts the first critical delay at T = 2 instead of 1.2092:
% it would declare T = 1.5 stable. See docs/tutorials/06_pade_pitfalls.md.
