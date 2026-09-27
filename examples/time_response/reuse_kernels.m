%% Several inputs, one preparation: a sine sweep through a delayed loop
% Run setup_contourroots once per session before this script.
% G is analytic in Re(s) > 0 because |s + 1| >= 1 > 0.5 there.

G = @(s) 1./(s + 1 + 0.5*exp(-s));
t = (0:0.02:20).';
K = ckernel(G, t, 'SingularityBound', 0);     % the inversions happen here

w = 0.5:0.25:2;                               % input frequencies (rad/s)
late = t > 10;                                % the transient has died out
measured = zeros(size(w));
for k = 1:numel(w)
    [y, ~, info] = clsim(K, sin(w(k)*t), t);  % no new inversion
    assert(info.converged && info.evaluations == 0)
    ab = [sin(w(k)*t(late)) cos(w(k)*t(late))] \ y(late);  % fit a sinusoid
    measured(k) = norm(ab);
end
predicted = abs(G(1i*w));
table(w.', measured.', predicted.', ...
      'VariableNames', {'omega', 'measured_gain', 'abs_G_jw'})
fprintf('%d transfer evaluations, all during preparation.\n', K.Info.evaluations);
