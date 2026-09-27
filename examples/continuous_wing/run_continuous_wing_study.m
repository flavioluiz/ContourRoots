function results = run_continuous_wing_study(outputDir)
%RUN_CONTINUOUS_WING_STUDY Flutter of the continuous Goland wing, validated.
%   RESULTS = RUN_CONTINUOUS_WING_STUDY() reproduces every number and figure
%   of Tutorial 11 and writes them to output/continuous_wing (about 3 min):
%     1. vacuum check: closed-form cantilever frequencies;
%     2. uniform strips: 1, 2, 4, 8 strips give the same characteristic function;
%     3. root locus 0-180 m/s, unstable-root counts and flutter speed;
%     4. independent check: finite elements + Hankel-form air loads;
%     5. a tapered wing, where the strips DO approximate the geometry;
%     6. time responses: vacuum pulse vs modal series, steps below/above
%        flutter, FFT vs adaptive quadrature.
%   The primary model never discretizes the structure or approximates
%   Theodorsen's function; the finite elements are only a reference.

    here = fileparts(mfilename('fullpath')); root = fileparts(fileparts(here));
    if exist('croots','file') ~= 2, run(fullfile(root,'setup_contourroots.m')); end
    oldPath = addpath(fullfile(root,'examples','models'));
    restorePath = onCleanup(@() path(oldPath));
    if nargin == 0, outputDir = fullfile(root,'output','continuous_wing'); end
    if ~isfolder(outputDir), mkdir(outputDir); end
    opts = {'AssumeAnalytic',true};
    results = struct;

    %% 1. Vacuum: closed-form frequencies of a clamped-free beam
    fprintf('1. Vacuum frequencies\n');
    dry = wing_model('dry'); e = dry.strips; L = dry.L;
    beta = [1.875104068711961 4.694091132974174 7.854757438237613];
    exact = sort([beta.^2*sqrt(e.EI/e.mu)/L^2, (1:2:21)*pi/(2*L)*sqrt(e.GJ/e.Ialpha)]).';
    exact = exact(exact < 900);
    [r, info] = croots(@(s) wing_delta(s,0,dry), [-2 2 1 900], opts{:});
    assert(info.complete && numel(r) == numel(exact));
    computed = sort(imag(r));
    results.vacuum = table(exact, computed, abs(computed-exact)./exact, ...
        'VariableNames', {'ClosedForm_rad_s','Contour_rad_s','RelativeError'});
    assert(max(results.vacuum.RelativeError) < 1e-10);

    %% 2. Uniform strips: the same function, whatever the number of strips
    fprintf('2. Strip invariance\n');
    z = [0.5, -3+60i, 2+140i, -10+300i];
    d1 = wing_delta(z, 130, wing_model('goland',1)); inv = zeros(4,2);
    for k = 1:4
        n = 2^(k-1); dn = wing_delta(z, 130, wing_model('goland',n));
        inv(k,:) = [n, max(abs(dn-d1)./abs(d1))];
    end
    results.strips = array2table(inv, 'VariableNames', {'Strips','MaxRelativeDifference'});
    assert(all(inv(:,2) < 1e-9));

    %% 3. Root locus, unstable-root counts and flutter speed
    fprintf('3. Root locus and flutter\n');
    wing = wing_model('goland');
    window = [-60 35 1 380]; rhp = [1e-3 60 -400 400];
    speeds = 0:10:180; nU = numel(speeds);
    poles = complex(NaN(nU,4)); unstable = NaN(nU,1); prev = [];
    for j = 1:nU
        U = speeds(j);
        [p, info] = croots(@(s) wing_delta(s,U,wing), window, opts{:});
        assert(info.complete && numel(p) == 4, 'Window must contain four modes.');
        if isempty(prev), [~,ix] = sort(imag(p)); p = p(ix); else, p = match(p, prev); end
        poles(j,:) = p.'; prev = p;
        if U > 0     % at U = 0 the undamped modes sit ON the imaginary axis
            [q, iq] = croots(@(s) wing_delta(s,U,wing), rhp, opts{:});
            assert(iq.complete); unstable(j) = numel(q);
        else
            unstable(j) = NaN;
        end
        fprintf('   U = %3d m/s: max Re = %8.4f, %d roots with Re > 0\n', U, max(real(p)), unstable(j));
    end
    results.locus = struct('speeds',speeds,'poles',poles,'unstable',unstable);
    alpha = @(U) max(real(croots(@(s) wing_delta(s,U,wing), window, opts{:})));
    Uf = fzero(alpha, [130 140], optimset('TolX',1e-9));
    p = croots(@(s) wing_delta(s,Uf,wing), window, opts{:}); [~,k] = max(real(p));
    neutral = neutral_point(@(s,U) wing_delta(s,U,wing), [Uf imag(p(k))]);
    results.flutter = struct('U',Uf,'omega',imag(p(k)),'neutralU',neutral(1), ...
        'neutralOmega',neutral(2));
    fprintf('   flutter: U = %.9f m/s, omega = %.9f rad/s\n', Uf, imag(p(k)));
    assert(abs(neutral(1)/Uf-1) < 1e-8);
    % Larger windows just below/above flutter: no other mode goes unstable.
    big = [];
    for f = [0.95 1.05]
        [pp, ib] = croots(@(s) wing_delta(s,f*Uf,wing), [-80 50 1 700], opts{:});
        assert(ib.complete);
        big(end+1,:) = [f*Uf, numel(pp), sum(real(pp) > 0)]; %#ok<AGROW>
    end
    results.largeWindow = array2table(big, 'VariableNames', {'U_m_s','ModesUpTo700','Unstable'});

    %% 4. Independent reference: finite elements + Hankel air loads
    fprintf('4. Finite-element reference\n');
    fe = [];
    for mesh = [8 8; 16 12; 32 20; 64 28; 128 36].'
        model = wing_fem(wing, mesh(1), mesh(2)); sc = sqrt(diag(model.K));
        D = @(s,U) det(wing_fem_matrix(s,U,model,'hankel')./(sc*sc.'));
        x = neutral_point(D, [Uf results.flutter.omega]);
        fe(end+1,:) = [model.nElements model.nModes x abs(x(1)/Uf-1)]; %#ok<AGROW>
        fprintf('   %2d elements, %2d modes: U = %.6f m/s\n', model.nElements, model.nModes, x(1));
    end
    results.fe = array2table(fe, 'VariableNames', ...
        {'Elements','Modes','U_m_s','Omega_rad_s','RelativeUError'});
    assert(fe(end,5) < 1e-4 && all(diff(fe(:,5)) < 0));

    %% 5. Tapered wing: strips now approximate the spanwise variation
    fprintf('5. Tapered wing\n');
    taper = [];
    guess = [Uf results.flutter.omega];
    for n = [2 4 8 16 32]
        wt = wing_model('tapered', n);
        guess = neutral_point(@(s,U) wing_delta(s,U,wt), guess);
        taper(end+1,:) = [n guess]; %#ok<AGROW>
    end
    results.taper = array2table(taper, 'VariableNames', {'Strips','U_m_s','Omega_rad_s'});

    %% 6. Time responses
    fprintf('6. Time responses\n');
    tol = {'AbsTol',1e-3,'RelTol',1e-3};          % outputs in mm: 1 micrometre
    % 6a. Vacuum, 1 kN half-sine tip pulse, against the modal series.
    t = (0:0.004:0.6).'; force = sin(pi*t/0.12).*(t <= 0.12);     % kN
    [wDry,~,iDry] = clsim(@(s) 1e6*wing_transfer(s,0,dry,1,1), force, t, ...
        'SingularityBound', 0, tol{:});
    assert(iDry.converged);
    series = 1e6*modal_series(dry, t, force, 150);
    results.timeDry = struct('t',t,'force',force,'contour',wDry,'series',series, ...
        'maxError',max(abs(wDry-series)),'maxResponse',max(abs(series)));
    fprintf('   vacuum pulse: max |error| = %.2e mm (response %.2f mm)\n', ...
        results.timeDry.maxError, results.timeDry.maxResponse);
    assert(results.timeDry.maxError < 1e-3*results.timeDry.maxResponse);
    % 6b. 1 kN tip step, below and above flutter (bound from step 3 counts).
    ts = (0:0.005:0.5).'; steps = zeros(numel(ts),2); speedsT = [120 150];
    for k = 1:2
        [q, iq] = croots(@(s) wing_delta(s,speedsT(k),wing), rhp, opts{:});
        assert(iq.complete && all(real(q) < 5));
        [y,~,iy] = cstep(@(s) 1e6*wing_transfer(s,speedsT(k),wing,1,1), ts, ...
            'SingularityBound', 5, tol{:});
        assert(iy.converged); steps(:,k) = y;
    end
    % 6c. The same step at 120 m/s by adaptive quadrature at a few times
    % (an independent inversion method). The de Hoog method is recorded
    % too: on this response, with many lightly damped beam modes, its
    % accelerated series settles on a slightly wrong plateau (known issue).
    tc = (0.1:0.1:0.5).'; ix = round(tc/0.005)+1;
    G120 = @(s) 1e6*wing_transfer(s,120,wing,1,1);
    [yq,~,iq] = cstep(G120, tc, 'SingularityBound', 5, 'Method', 'quadrature', tol{:});
    assert(iq.converged);
    yd = cstep(G120, tc, 'SingularityBound', 5, 'Method', 'dehoog', tol{:}, 'Warn', false);
    results.timeAero = struct('t',ts,'speeds',speedsT,'steps',steps, ...
        'checkTimes',tc,'quadrature',yq, ...
        'quadratureDifference',max(abs(yq-steps(ix,1))), ...
        'dehoogDifference',max(abs(yd-steps(ix,1))));
    assert(results.timeAero.quadratureDifference < 1e-3 + 1e-3*max(abs(steps(:,1))));
    fprintf('   FFT vs quadrature at 120 m/s: %.2e mm (de Hoog: %.2e mm)\n', ...
        results.timeAero.quadratureDifference, results.timeAero.dehoogDifference);

    %% Save tables and figures
    results.metadata = struct('matlab',version,'toolbox',contourroots_version);
    save(fullfile(outputDir,'continuous_wing_results.mat'),'results');
    for name = {'vacuum','strips','fe','taper','largeWindow'}
        writetable(results.(name{1}), fullfile(outputDir,[name{1} '.csv']));
    end
    figures(results, outputDir);
    fprintf('Results written to %s\n', outputDir);
end

function p = match(p, prev)
% Order the new roots like the previous ones (smallest total displacement).
    orders = perms(1:numel(p)); cost = zeros(size(orders,1),1);
    for k = 1:size(orders,1), cost(k) = sum(abs(p(orders(k,:))-prev)); end
    [~,k] = min(cost); p = p(orders(k,:));
end

function x = neutral_point(delta, x)
% Local Newton solve of Re/Im delta(i*omega,U) = 0 for x = [U omega].
% Used for the finite-element reference (Hankel loads exist only on the
% imaginary axis) and as a cross-check; it is not a regional search.
    x = x(:);
    for it = 1:40
        f = delta(1i*x(2),x(1)); F = [real(f); imag(f)]; J = zeros(2);
        for j = 1:2
            h = 1e-5*x(j); xp = x; xm = x; xp(j) = xp(j)+h; xm(j) = xm(j)-h;
            df = (delta(1i*xp(2),xp(1))-delta(1i*xm(2),xm(1)))/(2*h);
            J(:,j) = [real(df); imag(df)];
        end
        step = J\F; x = x - step;
        if norm(step./x) < 1e-12, break; end
    end
    assert(it < 40, 'Neutral-point iteration did not converge.');
    x = x.';
end

function y = modal_series(wing, t, force, nModes)
% Tip deflection of the vacuum cantilever from its exact bending modes
% (independent of the inversion): ramp response summed over FOH slopes.
    e = wing.strips(1); L = wing.L; beta = zeros(nModes,1);
    for j = 1:nModes
        c = (j-0.5)*pi;
        beta(j) = fzero(@(x) cos(x)+1/cosh(x), [max(0.01,c-0.7) c+0.7]);
    end
    omega = beta.^2*sqrt(e.EI/e.mu)/L^2;
    ramp = (4/(e.mu*L))*sum((t-sin(t*omega.')./omega.')./(omega.'.^2), 2);
    slopes = diff(force)/(t(2)-t(1)); w = [slopes(1); diff(slopes); 0];
    y = conv(w, ramp); y = y(1:numel(t));
end

function figures(r, outputDir)
    blue = [0 .35 .65]; orange = [.85 .35 .05]; green = [.35 .55 .25];
    % Root locus and spectral abscissa
    f = figure('Visible','off','Position',[100 100 1100 430]);
    subplot(1,2,1), hold on, box on, grid on
    U = r.locus.speeds; P = r.locus.poles; cmap = parula(numel(U));
    for j = 1:numel(U)
        plot(real(P(j,:)), imag(P(j,:)), 'o', 'MarkerFaceColor', cmap(j,:), ...
            'MarkerEdgeColor', 'none', 'MarkerSize', 6)
    end
    for k = 1:4, plot(real(P(:,k)), imag(P(:,k)), '-', 'Color', [.6 .6 .6]), end
    names = {'1st bending','1st torsion','2nd torsion','2nd bending/torsion'};
    [~, order] = sort(imag(P(1,:)));                  % still-air modes, U = 0
    for k = 1:4
        text(real(P(1,order(k)))+1.5, imag(P(1,order(k))), names{k}, 'FontSize', 9)
    end
    xline(0, 'k:'), xlabel('Re(s) [1/s]'), ylabel('Im(s) [rad/s]')
    colormap(parula), cb = colorbar; cb.Label.String = 'U [m/s]'; clim([U(1) U(end)])
    title('Modes of the continuous wing, U = 0 to 180 m/s')
    subplot(1,2,2), hold on, box on, grid on
    yyaxis left, plot(U, max(real(P),[],2), 'o-', 'Color', blue, 'MarkerFaceColor', blue)
    ylabel('largest Re(s) [1/s]'), yline(0, 'k:')
    yyaxis right, stairs(U, r.locus.unstable, '-', 'Color', orange, 'LineWidth', 1.5)
    ylabel('roots with Re(s) > 0'), ylim([-0.2 3])
    xline(r.flutter.U, '--', sprintf('U_f = %.2f m/s', r.flutter.U))
    xlabel('U [m/s]'), title('Stability: abscissa and unstable-root count')
    exportgraphics(f, fullfile(outputDir,'wing_root_locus.png'), 'Resolution', 150); close(f)
    % Finite-element convergence toward the continuous result
    f = figure('Visible','off','Position',[100 100 1100 400]);
    subplot(1,2,1)
    loglog(r.fe.Elements, r.fe.RelativeUError, 'o-', 'Color', blue, ...
        'MarkerFaceColor', blue, 'LineWidth', 1.2), grid on
    xlabel('finite elements (retained modes: 8 to 36)'), ylabel('|U_f^{FE}/U_f - 1|')
    title('Discretized model converges to the continuous one')
    subplot(1,2,2), plot(r.taper.Strips, r.taper.U_m_s, 's-', 'Color', green, ...
        'MarkerFaceColor', green, 'LineWidth', 1.2), grid on
    set(gca, 'XScale', 'log'), xticks(r.taper.Strips)
    xlabel('strips'), ylabel('U_f [m/s]'), title('Tapered wing: strips approximate the taper')
    exportgraphics(f, fullfile(outputDir,'wing_convergence.png'), 'Resolution', 150); close(f)
    % Time responses
    f = figure('Visible','off','Position',[100 100 1100 400]);
    subplot(1,2,1), d = r.timeDry;
    plot(d.t, d.contour, '-', 'Color', blue, 'LineWidth', 1.5), hold on
    plot(d.t(1:6:end), d.series(1:6:end), 'o', 'Color', orange), grid on
    legend('ContourRoots (clsim)', '150-mode analytic series', 'Location', 'southwest')
    xlabel('t [s]'), ylabel('tip deflection [mm]'), title('Vacuum: 1 kN half-sine pulse')
    subplot(1,2,2), a = r.timeAero;
    plot(a.t, a.steps(:,1), '-', 'Color', blue, 'LineWidth', 1.5), hold on
    plot(a.t, a.steps(:,2), '-', 'Color', orange, 'LineWidth', 1.5)
    plot(a.checkTimes, a.quadrature, 'ko'), grid on
    legend('120 m/s', '150 m/s', 'quadrature check', 'Location', 'northwest')
    xlabel('t [s]'), ylabel('tip deflection [mm]'), title('Air flow: 1 kN tip step')
    exportgraphics(f, fullfile(outputDir,'wing_time_response.png'), 'Resolution', 150); close(f)
end
