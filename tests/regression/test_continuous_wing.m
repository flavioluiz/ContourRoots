function tests = test_continuous_wing
%TEST_CONTINUOUS_WING Continuous Goland wing: structure, air loads, flutter, time.
%   Reference values come from closed-form beam theory, the independent
%   Hankel-form air loads and the study run_continuous_wing_study.
    tests = functiontests(localfunctions);
end

function setupOnce(tc)
    root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    tc.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,'matlab')));
    tc.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,'examples','models')));
end

function testStaticCompliance(tc)
    % s = 0 in vacuum: the classical cantilever flexibility, and reciprocity.
    dry = wing_model('dry'); e = dry.strips; L = dry.L;
    exact = [L^3/(3*e.EI) L^2/(2*e.EI) 0; L^2/(2*e.EI) L/e.EI 0; 0 0 L/e.GJ];
    G = zeros(3);
    for o = 1:3, for i = 1:3, G(o,i) = wing_transfer(0,0,dry,o,i); end, end
    tc.verifyLessThan(norm(G-exact,'fro')/norm(exact,'fro'), 1e-12);
end

function testVacuumFrequencies(tc)
    dry = wing_model('dry'); e = dry.strips; L = dry.L;
    exact = sort([[1.875104068711961 4.694091132974174].^2*sqrt(e.EI/e.mu)/L^2, ...
                  [1 3]*pi/(2*L)*sqrt(e.GJ/e.Ialpha)]);
    [r, info] = croots(@(s) wing_delta(s,0,dry), [-2 2 1 320], 'AssumeAnalytic', true);
    tc.verifyTrue(info.complete);
    tc.verifyEqual(numel(r), 4);
    tc.verifyLessThan(max(abs(sort(imag(r)).'-exact)./exact), 1e-10);
    tc.verifyLessThan(max(abs(real(r))), 1e-8);
end

function testUniformStripsAndSubstepsAreExact(tc)
    % Splitting a uniform strip changes the bookkeeping, not the model.
    z = [0.5, -3+60i, 2+140i];
    d1 = wing_delta(z, 130, wing_model('goland',1));
    for n = [2 4]
        dn = wing_delta(z, 130, wing_model('goland',n));
        tc.verifyLessThan(max(abs(dn-d1)./abs(d1)), 1e-9);
    end
    w = wing_model('goland'); s = 5+2e4i; g = zeros(1,2);
    for k = 1:2          % transfer at high frequency with 12 and 24 pieces
        K = wing_matrix(s,120,w,12*k); f = zeros(size(K,1),1);
        f(end-2) = 1/w.scale(4); x = K\f; g(k) = w.scale(1)*x(end-5);
    end
    tc.verifyLessThan(abs(diff(g))/abs(g(1)), 1e-10);
    tc.verifyLessThan(abs(wing_transfer(s,120,w,1,1)-g(1))/abs(g(1)), 1e-10);
end

function testAnalyticAndConjugate(tc)
    w = wing_model('goland'); s = -3+60i; h = 1e-4; U = 130;
    d = wing_delta(s,U,w);
    tc.verifyLessThan(abs(wing_delta(conj(s),U,w)-conj(d))/abs(d), 1e-11);
    dx = (wing_delta(s+h,U,w)-wing_delta(s-h,U,w))/(2*h);
    dy = (wing_delta(s+1i*h,U,w)-wing_delta(s-1i*h,U,w))/(2i*h);
    tc.verifyLessThan(abs(dx-dy)/abs(dx), 1e-6);      % Cauchy-Riemann
end

function testAirLoadsAgainstHankelForm(tc)
    e = wing_model('goland').strips;
    for U = [40 140 220]
        for omega = [5 70 300]
            A = theodorsen_loads(1i*omega,U,e.b,e.a,e.rho);
            B = theodorsen_loads_hankel(omega,U,e.b,e.a,e.rho);
            tc.verifyLessThan(norm(A-B,'fro')/norm(B,'fro'), 1e-12);
        end
    end
    tc.verifyEqual(theodorsen_loads(1+2i,100,e.b,e.a,0), zeros(2));   % vacuum
    tc.verifyError(@() theodorsen_loads(-1,100,e.b,e.a,e.rho), 'theodorsen:BranchCut');
end

function testStabilityCountsAndFlutter(tc)
    w = wing_model('goland'); rhp = [1e-3 60 -400 400]; opts = {'AssumeAnalytic',true};
    [p, info] = croots(@(s) wing_delta(s,130,w), rhp, opts{:});
    tc.verifyTrue(info.complete); tc.verifyEmpty(p);
    [p, info] = croots(@(s) wing_delta(s,145,w), rhp, opts{:});
    tc.verifyTrue(info.complete); tc.verifyEqual(numel(p), 2);
    window = [-60 35 1 380];
    alpha = @(U) max(real(croots(@(s) wing_delta(s,U,w), window, opts{:})));
    Uf = fzero(alpha, [130 140], optimset('TolX',1e-8));
    tc.verifyLessThan(abs(Uf-136.983977449), 1e-5);   % study value
end

function testFiniteElementReference(tc)
    % Independent discretize-then-truncate model with Hankel-form loads.
    w = wing_model('goland'); fem = wing_fem(w, 32, 20); sc = sqrt(diag(fem.K));
    D = @(s,U) det(wing_fem_matrix(s,U,fem,'hankel')./(sc*sc.'));
    x = [137; 70];
    for it = 1:30                          % Newton on Re/Im D(i*omega,U) = 0
        f = D(1i*x(2),x(1)); J = zeros(2);
        for j = 1:2
            h = 1e-5*x(j); dx = zeros(2,1); dx(j) = h;
            v = (D(1i*(x(2)+dx(2)),x(1)+dx(1))-D(1i*(x(2)-dx(2)),x(1)-dx(1)))/(2*h);
            J(:,j) = [real(v); imag(v)];
        end
        step = J\[real(f); imag(f)]; x = x-step;
        if norm(step./x) < 1e-12, break; end
    end
    tc.verifyLessThan(abs(x(1)/136.983977449-1), 2e-4);   % 32 elements: 1.5e-4
    tc.verifyGreaterThan(x(1), 136.983977449);            % converges from above
end

function testVacuumStepAgainstModalSeries(tc)
    dry = wing_model('dry'); e = dry.strips; L = dry.L; t = (0:0.004:0.2).';
    [y, ~, info] = cstep(@(s) 1e6*wing_transfer(s,0,dry,1,1), t, ...
        'SingularityBound', 0, 'AbsTol', 1e-3, 'RelTol', 1e-3);   % mm per kN
    beta = zeros(150,1);
    for j = 1:150
        c = (j-0.5)*pi; beta(j) = fzero(@(x) cos(x)+1/cosh(x), [max(0.01,c-0.7) c+0.7]);
    end
    omega = beta.^2*sqrt(e.EI/e.mu)/L^2;
    ref = 1e6*4/(e.mu*L)*sum((1-cos(t*omega.'))./(omega.'.^2), 2);
    tc.verifyTrue(info.converged);
    tc.verifyLessThan(max(abs(y-ref)), 1e-3 + 1e-3*max(abs(ref)));
end

function testModelValidation(tc)
    tc.verifyError(@() wing_model('unknown'), 'wing:Model');
    tc.verifyError(@() wing_model('goland',0), 'MATLAB:expectedPositive');
    t = wing_model('tapered',4);                 % strips follow the taper
    tc.verifyTrue(all(diff([t.strips.b]) < 0));
    tc.verifyLessThan(abs(sum([t.strips.length])-t.L), 1e-12);
end

function testDeHoogNoFalseConvergence(tc)
    % Regression (0.6.0): with a 1.4x degree ladder, de Hoog settled on a
    % plateau that omitted lightly damped high modes (8.0987 mm instead of
    % 8.0553 mm at t = 0.3 s) and still reported convergence. It must now
    % either agree with quadrature within tolerance or report unresolved.
    w = wing_model('goland'); G = @(s) 1e6*wing_transfer(s,120,w,1,1);
    opts = {'SingularityBound',5,'AbsTol',1e-3,'RelTol',1e-3,'Warn',false};
    [q,~,iq] = cstep(G, 0.3, opts{:}, 'Method', 'quadrature');
    [d,~,id] = cstep(G, 0.3, opts{:}, 'Method', 'dehoog');
    tc.verifyTrue(iq.converged);
    if id.converged
        tc.verifyLessThanOrEqual(abs(d-q), 2*(1e-3+1e-3*abs(q)));
    end
end
