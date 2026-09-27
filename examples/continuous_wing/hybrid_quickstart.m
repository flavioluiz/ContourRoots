%% Hybrid ports and MIMO responses of the continuous wing
% A strip of the wing as a two-port with mixed (hybrid) variables, two
% strips connected through their interface, the same transfer from the
% implicit assembly, and the simultaneous response to a tip force and a
% tip torque. Nothing is discretized: every block is the exact solution of
% the beam equations with exact Theodorsen strip loads (Tutorial 11).
% Run setup_contourroots once per session before this script.
repo = fileparts(fileparts(which('croots')));
addpath(fullfile(repo, 'examples', 'models'))

%% 1. One strip as a two-port: displacements in on the left, loads in on the right
wing = wing_model('goland', 2);                   % two identical half-span strips
U = 120;  s = 3 + 70i;                            % airspeed, one complex frequency
[HA, ia] = wing_hybrid(s, U, wing.strips(1), wing.scale);
[HB, ib] = wing_hybrid(s, U, wing.strips(2), wing.scale);
assert(ia.available && ib.available)

%% 2. Connect them: solve the interface, do not multiply the matrices
[H, ij] = wing_hybrid_join(HA, HB);
assert(ij.available)
Hwhole = wing_hybrid(s, U, wing_model('goland').strips, wing.scale);  % one full-span strip
fprintf('joined halves vs one strip: %.1e\n', norm(H - Hwhole, 'fro')/norm(Hwhole, 'fro'));
Ghybrid = H(4:6, 4:6);                            % clamped root (qL = 0), tip loads in

%% 3. The same transfer from the implicit assembly: one solve, nine channels
M = wing_cdyn(U, wing);                           % K(s) x = B u,  y = C x
[Gimplicit, work] = ceval(M, s);
fprintf('hybrid vs implicit: %.1e, %d factorization for %d channels\n', ...
    norm(Ghybrid - Gimplicit, 'fro')/norm(Gimplicit, 'fro'), work.factorizations, numel(Gimplicit));

%% 4. Where the hybrid form breaks down: an artificial pole
% In vacuum, a half-span piece clamped at one end and free at the other
% resonates at 4 x 49.49 rad/s. There the element's hybrid matrix does not
% exist, although the whole wing has no mode there.
dry = wing_model('dry', 2);  e = dry.strips(1);
w0 = 1.875104068711961^2*sqrt(e.EI/e.mu)/e.length^2;       % 197.97 rad/s
[~, info] = wing_hybrid(1i*w0, 0, e, dry.scale);
G = ceval(wing_cdyn(0, dry), 1i*w0);
fprintf('at %.2f rad/s: hybrid available = %d, implicit |G11| = %.2e (finite)\n', ...
    w0, info.available, abs(G(1,1)));

%% 5. Tip force and tip torque together
% Inputs in kN and kN m, outputs in mm and mrad. The bank holds the kernels
% of all four channels; each clsim call is then a cheap convolution.
% SingularityBound = 5 is an assumption, supported by the root counts of
% Tutorial 11 (none unstable at 120 m/s), not a proof.
P = wing_response_model(U, wing, 'implicit', [1 3]);   % force/torque -> w/alpha
t = (0:0.005:0.4).';
loads = [1 + 0.2*sin(10*t), 0.3*cos(5*t)];
bank = ckernel(P, t, 'SingularityBound', 5, 'AbsTol', 1e-3, 'RelTol', 1e-3);
out = {'AbsTol', [0.01 0.01], 'RelTol', 2e-3};    % 10 micrometres, 10 microradians
[y, ~, info] = clsim(bank, loads, t, out{:});
yForce  = clsim(bank, [loads(:,1) 0*t], t, out{:});
yTorque = clsim(bank, [0*t loads(:,2)], t, out{:});
assert(info.converged && info.evaluations == 0)
fprintf('superposition residual: %.1e\n', max(abs(y - yForce - yTorque), [], 'all'));

figure, tiledlayout(1, 2)
names = {'tip deflection [mm]', 'tip twist [mrad]'};
for k = 1:2
    nexttile, plot(t, yForce(:,k), t, yTorque(:,k), t, y(:,k), 'k', 'LineWidth', 1.2)
    grid on, xlabel('t [s]'), ylabel(names{k})
end
legend('from the force', 'from the torque', 'both')
