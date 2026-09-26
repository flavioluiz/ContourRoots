function tests=test_response_api
    tests=functiontests(localfunctions);
end
function setupOnce(tc)
    root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
    tc.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,'matlab')));
end
function testShapeAndConstants(tc)
    t=0:.1:1;
    [y,tt,i]=cstep(2,t); tc.verifySize(y,[11 1]); tc.verifyEqual(tt,t.');
    tc.verifyEqual(y,2*ones(11,1)); tc.verifyTrue(i.converged);
    [g,~,i]=cimpulse(2,t); tc.verifyEqual(g,zeros(11,1));
    tc.verifyEqual(i.singularTerms.weight,2); tc.verifyEqual(i.singularTerms.time,0);
    [y,~,i]=clsim(0,sin(t),t); tc.verifyEqual(y,zeros(11,1)); tc.verifyTrue(i.converged);
end
function testContracts(tc)
    f=@(s) 1./(s+1); t=0:.1:1;
    tc.verifyError(@() cstep(f,t),'ContourRoots:ResponseDomain');
    tc.verifyError(@() cimpulse(f,t,'SingularityBound',0),'ContourRoots:ResponseImpulse');
    tc.verifyError(@() cstep(f,t,'Abscissa',-1,'SingularityBound',0),'ContourRoots:ResponseDomain');
    tc.verifyError(@() cstep(ndpair(1,[1 -1]),t,'AssumeStable',true),'ContourRoots:ResponseDomain');
    tc.verifyError(@() cstep(f,t,'Abscissa',0),'ContourRoots:ResponseDomain');
    tc.verifyError(@() cstep(ndpair([1 0],1),t),'ContourRoots:ResponseModel');
end
function testTimesAndOptions(tc)
    G=ndpair(1,[1 1]);
    tc.verifyError(@() cstep(G,[0 .1 .4]),'ContourRoots:ResponseTime');
    tc.verifyError(@() cstep(G,[0 0]),'ContourRoots:ResponseTime');
    tc.verifyError(@() clsim(G,[1 1],[1 2]),'ContourRoots:ResponseTime');
    tc.verifyError(@() cstep(G,[0 1],'Typo',1),'ContourRoots:ResponseOption');
    tc.verifyError(@() cstep(G,[0 1],'Warn',NaN),'ContourRoots:ResponseOption');
    tc.verifyError(@() cstep(G,[0 1],'Feedthrough',2),'ContourRoots:ResponseModel');
end
function testScalarOnlyHandle(tc)
    f=@(s) 1/prod([s+1 s+2]); t=(0:.1:2).';
    [g,~,i]=cimpulse(f,t,'SingularityBound',0,'RegularImpulse',true,'InitialValue',0,'Method','dehoog');
    tc.verifyTrue(i.converged); tc.verifyLessThan(max(abs(g-(exp(-t)-exp(-2*t)))),1e-6);
end
function testSymmetryAndOverflow(tc)
    tc.verifyError(@() cstep(@(s) 1i./(s+1),[0 .1 1],'Method','dehoog','SingularityBound',0),'ContourRoots:ResponseSymmetry');
    tc.verifyError(@() cstep(@(s) NaN*s,[0 1],'SingularityBound',0),'ContourRoots:ResponseEvaluation');
    tc.verifyError(@() cstep(@(s) 1e-15i./(s+1),[0 .1 .2],'SingularityBound',0),'ContourRoots:ResponseSymmetry');
end
function testBudgetAndUnknownOrigin(tc)
    [y,~,i]=cstep(ndpair(1,[1 1]),[0 .1 1],'Method','dehoog','MaxPoints',8,'Warn',false);
    tc.verifyFalse(i.converged); tc.verifyTrue(any(isnan(y)));
    [g,~,i]=cinvlaplace(@(s) 1./sqrt(s),[0 .1 1],'Method','dehoog','SingularityBound',0,'Warn',false);
    tc.verifyTrue(isnan(g(1))); tc.verifyFalse(i.resolvedMask(1)); tc.verifyTrue(all(i.resolvedMask(2:end)));
end
function testQuadratureFailures(tc)
    [y,~,i]=cstep(ndpair(1,[1 1]),[0 .1 1],'Method','quadrature', ...
        'MaxPoints',8,'Warn',false);
    tc.verifyFalse(i.converged); tc.verifyTrue(all(isnan(y(2:end))));
    tc.verifyError(@() cstep(@(s) 1i./(s+1),[0 .1 1], ...
        'Method','quadrature','SingularityBound',0),'ContourRoots:ResponseSymmetry');
end
function testPlotsDoNotChangeUserFigures(tc)
    f=figure('Visible','off'); tc.addTeardown(@() close(f)); ax=axes(f);
    [~,~,i]=cstep(2,[0 1],'Parent',ax,'Plot',true);
    tc.verifyTrue(isgraphics(f)); tc.verifyFalse(ishold(ax)); tc.verifyTrue(i.converged);
    hold(ax,'on'); cimpulse(2,[0 1],'Parent',ax); tc.verifyTrue(ishold(ax));
end
function testSymbolicAdapters(tc)
    tc.assumeTrue(~isempty(ver('symbolic')),'Symbolic Math Toolbox not available.');
    syms s z
    t=(0:.1:1).';
    [y,~,i]=cstep(1/(s+1),t,'SingularityBound',0,'Method','dehoog');
    tc.verifyTrue(i.converged); tc.verifyLessThan(max(abs(y-(1-exp(-t)))),1e-6);
    tc.verifyError(@() cstep(ndpair(s+1,z+2),t,'SingularityBound',0),'ContourRoots:ResponseModel');
end
