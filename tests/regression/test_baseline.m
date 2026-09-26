function tests = test_baseline
%TEST_BASELINE Compare current results with the research-repository baseline.
%   baseline_57af5bf.mat was produced by the research code at commit 57af5bf,
%   before the reorganization into ContourRoots. Values must agree within
%   numerical tolerances (not bitwise).
    tests = functiontests(localfunctions);
end

function setupOnce(tc)
    here = fileparts(mfilename('fullpath'));
    repo = fileparts(fileparts(here));
    import matlab.unittest.fixtures.PathFixture
    tc.applyFixture(PathFixture(fullfile(repo,'matlab')));
    tc.applyFixture(PathFixture(fullfile(repo,'matlab','delay')));
    tc.applyFixture(PathFixture(fullfile(repo,'examples','models')));
    data = load(fullfile(here,'baseline_57af5bf.mat'));
    tc.TestData.b = data.b;
end

function testCriticalDelays(tc)
    b = tc.TestData.b;
    compareTables(critical_delays(2,[1 1],10), b.p1critical);
    compareTables(critical_delays(3,[1 .8 4],10), b.p3critical);
    compareTables(pade_critical_delays(2,[1 1],10,[1 2 4 8]), b.p1pade);
    compareTables(pade_critical_delays(3,[1 .8 4],10,[1 2 4 8]), b.p3pade);
end

function testBoundsAndCounts(tc)
    b = tc.TestData.b;
    assert(abs(rhp_root_bound(2,[1 1])-b.p1bound) < 1e-12);
    assert(abs(rhp_root_bound(3,[1 .8 4])-b.p3bound) < 1e-12);
    Z = arrayfun(@(t) unstable_root_count(3,[1 .8 4],t), b.p3ZDelays);
    assert(isequal(Z, b.p3Z));
end

function testDelayRoots(tc)
    b = tc.TestData.b;
    [r,info] = delay_roots(2,[1 1],1,[-8 2],[-20 20],36,60,'MaxRefinements',2);
    assert(info.argumentPrincipleCount == b.delayRootsCount);
    sameSet(r, b.delayRoots, 1e-8);
end

function testGeneralSolver(tc)
    b = tc.TestData.b;
    Delta = @(s) 1+s+s.^2+(2*s+3).*exp(-s);
    [r,info] = characteristic_roots(Delta,[-8 2 -20 20],'AssumeAnalytic',true);
    assert(strcmp(info.status,b.userStatus));
    assert(isequal(info.multiplicity,b.userMultiplicity));
    sameSet(r, b.userRoots, 1e-9);
    kinds = fieldnames(b.distributed);
    for j = 1:numel(kinds)
        ref = b.distributed.(kinds{j});
        [G,meta] = distributed_model(kinds{j},'Sensor',.5);
        [p,ii] = transfer_poles(G,ref.region,'AssumeAnalytic',true, ...
            'Singularities',meta.singularities);
        assert(strcmp(ii.status,ref.status), 'Status changed: %s', kinds{j});
        sameSet(p, ref.poles, 1e-9);
        sameSet(ii.cancelledLocations, ref.cancelled, 1e-9);
    end
end

function testCoupledBeam(tc)
    b = tc.TestData.b;
    [D,G] = coupled_beam_model();
    [r,info] = characteristic_roots(D,[-30 1 -130 130],'AssumeAnalytic',true);
    assert(strcmp(info.status,b.beamStatus));
    sameSet(r, b.beamRoots, 1e-9);
    p = transfer_poles(G.TipForceToTip,[-30 1 -130 130],'AssumeAnalytic',true);
    sameSet(p, b.beamTipPoles, 1e-9);
end

function sameSet(a,b,tol)
    assert(numel(a)==numel(b), 'Different number of values (%d vs %d).', numel(a), numel(b));
    for z = b(:).'
        assert(min(abs(a-z)) <= tol*(1+abs(z)), 'Value %s moved.', num2str(z));
    end
end

function compareTables(a,b)
    assert(height(a)==height(b), 'Different number of rows.');
    names = intersect(a.Properties.VariableNames, {'Frequency','Delay','Order','Branch'});
    for k = 1:numel(names)
        x = a.(names{k}); y = b.(names{k});
        assert(max(abs(x-y)./(1+abs(y))) < 1e-9, 'Column %s changed.', names{k});
    end
end
