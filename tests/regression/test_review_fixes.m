function tests = test_review_fixes
%TEST_REVIEW_FIXES Regression tests for the problems found in the v0.1.0 review.
    tests = functiontests(localfunctions);
end

function setupOnce(tc)
    repo = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    import matlab.unittest.fixtures.PathFixture
    tc.applyFixture(PathFixture(fullfile(repo,'matlab')));
    tc.applyFixture(PathFixture(fullfile(repo,'matlab','delay')));
    tc.applyFixture(PathFixture(fullfile(repo,'tools')));
end

function testScalarOnlyHandleIsNotTreatedAsConstant(tc)
    % A non-vectorized handle returns a scalar for vector input.
    F = @(s) prod([s-1, s+1]);
    [r,info] = croots(F,[-2 2 -1 1],'AssumeAnalytic',true);
    tc.verifyTrue(info.complete);
    tc.verifyEqual(sort(real(r)),[-1;1],'AbsTol',1e-8);
    % Same roots as the vectorized form, also for poles of a pair.
    G1 = ndpair(@(s) prod([s-0.5, s+0.5]), @(s) prod([s-0.3i, s+0.3i, s+2]));
    G2 = ndpair(@(s) s.^2-0.25, @(s) (s.^2+0.09).*(s+2));
    p1 = cpoles(G1,[-3 1 -1 1],'AssumeAnalytic',true);
    p2 = cpoles(G2,[-3 1 -1 1],'AssumeAnalytic',true);
    tc.verifyEqual(numel(p1),3); tc.verifyEqual(numel(p2),3);
    tc.verifyLessThan(max(arrayfun(@(z) min(abs(p1-z)),p2)),1e-8);
    % A genuinely constant handle has no roots.
    [r,info] = croots(@(s) 3,[-1 1 -1 1],'AssumeAnalytic',true);
    tc.verifyEmpty(r); tc.verifyTrue(info.complete);
end

function testCriticalDelaysScaleInvariant(tc)
    a = critical_delays(2,[1 1],10);
    for c = [1e-8 1e8]
        b = critical_delays(2*c,[c c],10);
        tc.verifyEqual(b.Delay,a.Delay,'RelTol',1e-12);
        tc.verifyEqual(b.Frequency,a.Frequency,'RelTol',1e-12);
    end
    [~,info] = critical_delays(2e-12,[1e-12 1e-12],10);
    tc.verifyFalse(info.zeroRootForAllDelays);          % F(0) = 3e-12, not 0
    [~,info] = critical_delays(1,[1 -1],10);
    tc.verifyTrue(info.zeroRootForAllDelays);           % F(0) = 0 exactly
end

function testCriticalDelaysTimeScaleInvariant(tc)
    % Time scaled by 1e6: frequencies divided and delays multiplied by 1e6.
    a = critical_delays(2,[1 1],10);
    b = critical_delays(2e-6,[1 1e-6],1e7);
    tc.verifyEqual(b.Delay,1e6*a.Delay,'RelTol',1e-10);
    tc.verifyEqual(b.Frequency,1e-6*a.Frequency,'RelTol',1e-10);
    tc.verifyEqual(b.Delay(1),2*pi/(3*sqrt(3))*1e6,'RelTol',1e-12);
end

function testDegenerateAndPersistentCases(tc)
    [crit,info] = tc.verifyWarning(@() critical_delays(1,1,10),'critical_delays:Degenerate');
    tc.verifyEmpty(crit); tc.verifyTrue(info.degenerate);
    % D = (s+1)(s^2+1), N = 2(s^2+1): s = +-i is a root for every delay.
    [crit,info] = tc.verifyWarning(@() critical_delays([2 0 2],[1 1 1 1],10), ...
        'critical_delays:PersistentImaginaryRoots');
    tc.verifyEqual(info.persistentFrequencies,1,'AbsTol',1e-8);
    tc.verifyEqual(crit.Frequency,sqrt(3)*ones(3,1),'RelTol',1e-12);   % s+1+2e^{-sT}
end

function testUnstableCountNeverFakesZero(tc)
    % Small delays: resolved and exact.
    tc.verifyEqual(unstable_root_count(2,[1 1],1),0);
    tc.verifyEqual(unstable_root_count(2,[1 1],100),56);
    % Large delay with roots crowding the imaginary axis: inconclusive, not 0.
    [Z,info] = unstable_root_count(2,[1 1],1e5,[],2^14);
    tc.verifyTrue(isnan(Z)); tc.verifyFalse(info.resolved);
    tc.verifyWarning(@() unstable_root_count(2,[1 1],1e5,[],2^14), ...
        'unstable_root_count:Inconclusive');
    % A root exactly on the imaginary axis (critical delay of P1).
    Tc = 2*pi/(3*sqrt(3));
    [Z,~] = unstable_root_count(2,[1 1],Tc,[],2^14);
    tc.verifyTrue(isnan(Z) || Z == 0);   % boundary root: never a wrong positive
    % delay_root_count: inconclusive instead of max(0,...).
    [c,~,i] = delay_root_count(2,[1 1],1,[-1 1],[-sqrt(3) 2],64,2^12);   % no root on edges
    tc.verifyTrue(i.resolved); tc.verifyGreaterThanOrEqual(c,0);
end

function testPadeCharacteristicDegree(tc)
    for T = [0.1 1 10]
        for n = [1 4 8 12 16]
            p = pade_characteristic(2,[1 1],T,n);
            tc.verifyEqual(numel(p)-1, 1+n, sprintf('T=%g, n=%d',T,n));
        end
    end
end

function testPadeCrossingsDoNotDisappear(tc)
    a = pade_critical_delays(2,[1 1],10,[1 2 4 8]);
    b = pade_critical_delays(2,[1 1],1e6,[1 2 4 8]);
    tc.verifyEqual(height(a),7);   % 1 + 1 + 2 + 3 events up to T = 10
    for k = 1:height(a)
        match = b(b.Order==a.Order(k),:);
        tc.verifyLessThan(min(abs(match.Delay-a.Delay(k))),1e-9*(1+a.Delay(k)));
    end
    % Order n has at most ceil((n*pi-theta)/(2*pi)) crossings for all T.
    tc.verifyEqual(sum(b.Order==8),4);
    % Known values: [1/1] gives T = 2, [2/2] gives sqrt(5)-1.
    tc.verifyEqual(b.Delay(b.Order==1),2,'AbsTol',1e-12);
    tc.verifyEqual(b.Delay(b.Order==2),sqrt(5)-1,'AbsTol',1e-12);
end

function testPartialCancellationMetadata(tc)
    % (s+0.2)/(s+0.2)^2: one pole remains at -0.2 (partial cancellation).
    [p,info] = cpoles(ndpair([1 .2],[1 .4 .04]),[-1 1 -1 1]);
    tc.verifyEqual(real(p),-0.2,'AbsTol',1e-8);
    tc.verifyEqual(info.multiplicity,1);
    tc.verifyEqual(info.cancellationOrders,1);
    tc.verifyFalse(info.cancellationComplete);
    % Total cancellation is flagged as complete.
    [~,info] = cpoles(ndpair(@(s) sinh(s/2),@(s) sinh(s)),[-1 1 -10 10],'AssumeAnalytic',true);
    tc.verifyTrue(all(info.cancellationComplete));
end

function testPoleWarningMentionsNdpair(tc)
    w = warning('off','all'); restore = onCleanup(@() warning(w)); %#ok<NASGU>
    warning('on','ContourRoots:Exploratory');
    lastwarn('');
    cpoles(@(s) (s-0.4)./(s+0.2),[-1 1 -1 1],'GridSize',[9 9]);
    [msg,id] = lastwarn;
    tc.verifyEqual(id,'ContourRoots:Exploratory');
    tc.verifySubstring(msg,'ndpair');
    tc.verifySubstring(msg,'1./G');
end

function testDocRunnerKeepsUserFigures(tc)
    f = figure('Visible','off');
    tc.addTeardown(@() close(f(isvalid(f))));
    check_doc_snippets({});
    tc.verifyTrue(isvalid(f));
end
