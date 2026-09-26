function tests = test_delay_tools
%TEST_DELAY_TOOLS Regression tests of the time-delay module (matlab/delay).
    tests = functiontests(localfunctions);
end

function setupOnce(tc)
    repo = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    import matlab.unittest.fixtures.PathFixture
    tc.applyFixture(PathFixture(fullfile(repo,'matlab')));
    tc.applyFixture(PathFixture(fullfile(repo,'matlab','delay')));
end

function testPadeFirstOrder(~)
    [num,den] = pade_delay(2,1);
    assert(norm(num-[-1 1]) < 1e-13);
    assert(norm(den-[1 1]) < 1e-13);
end

function testZeroDelayCollapsesToPolynomial(~)
    % At T=0, the quasi-polynomial collapses exactly to D+N.
    [r,info] = delay_roots(2,[1 1],0,[-5 2],[-3 3],20,20, ...
        'PadeSeedOrders',1,'MaxRefinements',1,'VerifyCount',true);
    assert(numel(r)==1 && abs(r+3)<1e-8);
    assert(info.countMatches && max(info.residuals)<1e-8);
end

function testConjugateSymmetryAndResiduals(~)
    [r,info] = delay_roots(2,[1 1],1,[-8 2],[-20 20],36,60, ...
        'MaxRefinements',2,'VerifyCount',true);
    assert(info.countMatches);
    assert(max(info.residuals)<1e-6);
    for k=1:numel(r)
        assert(min(abs(r-conj(r(k))))<2e-5);
    end
end

function testClosedFormCriticalDelayP1(~)
    [crit,cinfo] = critical_delays(2,[1 1],2);
    expectedT = 2*pi/(3*sqrt(3));
    assert(~cinfo.zeroRootForAllDelays && ~isempty(crit));
    assert(abs(crit.Delay(1)-expectedT)<1e-11);
    assert(abs(crit.Frequency(1)-sqrt(3))<1e-11);
    assert(crit.CrossingSpeed(1)>0 && crit.Residual(1)<1e-10);
end

function testDelayIndependentStabilityP2(~)
    % |2+i*omega| can never equal one: P2 is stable for every delay.
    assert(isempty(critical_delays(1,[1 2],5)));
end

function testCrossingResidualsP3(~)
    crit = critical_delays(3,[1 0.8 4],3);
    assert(~isempty(crit));
    assert(max(crit.Residual)<1e-9);
end

function testPadeKeepsFrequencyShiftsDelay(~)
    expectedT = 2*pi/(3*sqrt(3));
    padeCrit = pade_critical_delays(2,[1 1],3,[1 2 4]);
    p1 = padeCrit(padeCrit.Order==1,:);
    assert(height(p1)==1 && abs(p1.Frequency-sqrt(3))<1e-11);
    assert(abs(p1.Delay-2)<1e-9);          % [1/1] prediction versus exact 1.2092...
    p4 = padeCrit(padeCrit.Order==4,:);
    assert(abs(p4.Delay(1)-expectedT)<2e-5);
    assert(max(padeCrit.Residual)<1e-8);
end

function testBoundBasedUnstableCounts(~)
    assert(abs(rhp_root_bound(2,[1 1])-3)<1e-12);
    assert(unstable_root_count(2,[1 1],1.0)==0);
    assert(unstable_root_count(2,[1 1],1.5)==2);
    assert(unstable_root_count(3,[1 .8 4],1.0)==2);
    assert(unstable_root_count(3,[1 .8 4],2.75)==0);   % P3 stability window
    assert(unstable_root_count(3,[1 .8 4],3.0)==2);
    assert(unstable_root_count(1,[1 2],7.3)==0);
end
