function tests = test_api
%TEST_API Tests of the public MATLAB-style interface.
    tests = functiontests(localfunctions);
end

function setupOnce(tc)
    repo = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    import matlab.unittest.fixtures.PathFixture
    tc.applyFixture(PathFixture(fullfile(repo,'matlab')));
    tc.applyFixture(PathFixture(fullfile(repo,'matlab','delay')));
end

function setup(tc)
    % Figures created by plotting tests are closed after each test.
    tc.addTeardown(@() close(findall(groot,'Type','figure','Tag','cr_test')));
end

function testCrootsMatchesCharacteristicRoots(tc)
    F = @(s) 1+s+s.^2+(2*s+3).*exp(-s);
    [r,info] = croots(F,[-8 2 -20 20],'AssumeAnalytic',true);
    [r0,info0] = characteristic_roots(F,[-8 2 -20 20],'AssumeAnalytic',true);
    tc.verifyEqual(r,r0);
    tc.verifyEqual(info.status,info0.status);
    tc.verifyTrue(info.complete);
    tc.verifyEqual(numel(r),7);
end

function testCrootsPolynomialLikeRoots(tc)
    [r,info] = croots([1 0 -1],[-2 2 -1 1]);
    tc.verifyTrue(info.complete);
    tc.verifyEqual(sort(real(r)),[-1;1],'AbsTol',1e-8);
    [r,info] = croots([1 -1],[-2 2 -1 1]);        % single root at s=1
    tc.verifyTrue(info.complete);
    tc.verifyEqual(r,1,'AbsTol',1e-8);
end

function testCpolesWithCancellation(tc)
    G = ndpair(@(s) sinh(s/2), @(s) sinh(s));
    [p,info] = cpoles(G,[-1 1 -10 10],'AssumeAnalytic',true);
    tc.verifyTrue(info.complete);
    expected = 1i*pi*[-3 -1 1 3];
    tc.verifyEqual(numel(p),4);
    for z = expected
        tc.verifyLessThan(min(abs(p-z)),1e-6);
    end
    tc.verifyEqual(numel(info.cancelledLocations),3);   % 0 and +-2*pi*i
end

function testCzerosAfterCancellation(tc)
    G = ndpair(@(s) (s+1).*exp(-s), [1 3 2]);
    [z,info] = czeros(G,[-3 1 -3 3],'AssumeAnalytic',true);
    tc.verifyTrue(info.complete);
    tc.verifyEmpty(z);
    [p,info] = cpoles(G,[-3 1 -3 3],'AssumeAnalytic',true);
    tc.verifyTrue(info.complete);
    tc.verifyEqual(p,-2,'AbsTol',1e-8);
end

function testPolynomialPairIsAnalyticWithoutAssertion(tc)
    G = ndpair([1 2],[1 4 3]);        % (s+2)/((s+1)(s+3))
    [p,info] = cpoles(G,[-4 1 -1 1]);
    tc.verifyTrue(info.complete);
    tc.verifyEqual(sort(real(p)),[-3;-1],'AbsTol',1e-8);
    z = czeros(G,[-4 1 -1 1]);
    tc.verifyEqual(z,-2,'AbsTol',1e-8);
end

function testExploratoryWarnsWithoutInfo(tc)
    F = @(s) s.^2+1;
    tc.verifyWarning(@() croots(F,[-2 2 -2 2]),'ContourRoots:Exploratory');
    tc.verifyWarningFree(@() croots(F,[-2 2 -2 2],'Warn',false));
    [~,info] = croots(F,[-2 2 -2 2]);          % no warning when info is requested
    tc.verifyEqual(info.status,'exploratory');
    tc.verifyFalse(info.complete);
end

function testIncompleteWarns(tc)
    F = @(s) s;                                % root on the boundary
    tc.verifyWarning(@() croots(F,[-1 0 -1 1],'AssumeAnalytic',true), ...
        'ContourRoots:Incomplete');
end

function testInputErrors(tc)
    tc.verifyError(@() croots(@(s) s),'ContourRoots:Region');
    tc.verifyError(@() croots(@(s) s,[0 1]),'ContourRoots:Region');
    tc.verifyError(@() croots(@(s) s,[-1 1 -1 1],'Mode','poles'),'ContourRoots:Mode');
    tc.verifyError(@() ndpair(@(s) s),'ContourRoots:ndpair');
    tc.verifyError(@() ndpair('s',1),'ContourRoots:ndpair');
    tc.verifyError(@() ndpair([0 0],1),'ContourRoots:ndpair');
end

function testCpzmapOutputsAndPlot(tc)
    G = ndpair(@(s) s+0.5, @(s) (s+0.5).*(s+2).*(s.^2+1));
    [p,z,ip,iz] = cpzmap(G,[-3 1 -2 2],'AssumeAnalytic',true);
    tc.verifyTrue(ip.complete && iz.complete);
    tc.verifyEmpty(z);                          % s=-0.5 cancels
    tc.verifyEqual(numel(p),3);
    f = figure('Visible','off','Tag','cr_test');
    ax = axes(f);
    cpzmap(ax,G,[-3 1 -2 2],'AssumeAnalytic',true);
    lines = findobj(ax,'Type','line','DisplayName','Poles');
    tc.verifyEqual(numel(lines(1).XData),3);
    tc.verifyNotEmpty(findobj(ax,'Type','line','DisplayName','Cancelled'));
end

function testVersionAndOverview(tc)
    tc.verifyMatches(contourroots_version(),'^\d+\.\d+\.\d+$');
    out = evalc('contourroots');
    tc.verifySubstring(out,'croots');
end
