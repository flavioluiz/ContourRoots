function tests=test_kernel_reuse
    tests=functiontests(localfunctions);
end
function setupOnce(tc)
    root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
    tc.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,'matlab')));
end
function testNoEvaluationAndFrozenClosure(tc)
    t=(0:.05:2).'; calls=0; rate=1; enabled=true;
    [K,setup]=ckernel(@model,t,'SingularityBound',0);
    setupCalls=calls; tc.verifyGreaterThan(setupCalls,0);
    rate=2; enabled=false;
    for u={ones(size(t)),sin(t),@(q) .5+sin(q)}
        [~,~,i]=clsim(K,u{1},t);
        tc.verifyTrue(i.converged&&i.kernelReused);
        tc.verifyEqual(i.evaluations,0);
        tc.verifyEqual(i.kernelPreparationEvaluations,setup.evaluations);
        tc.verifyEqual(calls,setupCalls);
    end
    [y,~,i]=clsim(K,ones(size(t)),t);
    tc.verifyLessThan(max(abs(y-(1-exp(-t)))),2e-6);
    tc.verifyEqual(i.stepKernelInfo.evaluations,0);
    enabled=true; K2=ckernel(@model,t,'SingularityBound',0);
    y2=clsim(K2,ones(size(t)),t);
    tc.verifyLessThan(max(abs(y2-(1-exp(-2*t))/2)),2e-6);
    function v=model(s)
        assert(enabled,'Transfer must not be evaluated on reuse.');
        calls=calls+1; v=1./(s+rate);
    end
end
function testFreshAgreementAllMethodsAndHolds(tc)
    t=(0:.1:1).'; G=ndpair([.4 1.4],[1 1]);
    for method={'fft','dehoog','quadrature'}
        for hold={'foh','zoh'}
            K=ckernel(G,t,'Method',method{1},'Interpolation',hold{1});
            for u={sin(t),1+cos(2*t)}
                [a,~,ia]=clsim(K,u{1},t);
                [b,~,ib]=clsim(G,u{1},t,'Method',method{1},'Interpolation',hold{1});
                tc.verifyTrue(ia.converged&&ib.converged);
                tc.verifyLessThan(max(abs(a-b)),2e-6);
                tc.verifyFalse(ib.kernelReused); tc.verifyEqual(ia.evaluations,0);
            end
        end
    end
end
function testNonrationalAndUnstable(tc)
    t=(0:.05:2).';
    plants={@(s) exp(-.7*sqrt(s)),@(s) 1./(s+1+.5*exp(-s)),ndpair(1,[1 -.3])};
    bounds=[0 0 .3]; u=1+sin(t);
    for k=1:numel(plants)
        K=ckernel(plants{k},t,'SingularityBound',bounds(k));
        [a,~,ia]=clsim(K,u,t);
        [b,~,ib]=clsim(plants{k},u,t,'SingularityBound',bounds(k));
        tc.verifyTrue(ia.converged&&ib.converged);
        tc.verifyLessThan(max(abs(a-b)),5e-6);
    end
end
function testConstantsAndZero(tc)
    t=(0:.1:1).'; u=sin(t);
    for gain=[0 2]
        K=ckernel(gain,t); [y,~,i]=clsim(K,u,t);
        tc.verifyEqual(y,gain*u); tc.verifyTrue(i.converged);
        tc.verifyEqual(K.Info.evaluations,0); tc.verifyEqual(i.evaluations,0);
    end
end
function testRechecksOutputAccuracy(tc)
    t=(0:.05:2).'; K=ckernel(ndpair(1,[1 1]),t);
    [~,~,i]=clsim(K,1e6*sin(20*t),t,'AbsTol',1e-20,'RelTol',1e-16,'Warn',false);
    tc.verifyFalse(i.converged); tc.verifyTrue(any(~i.resolvedMask));
    tc.verifyGreaterThan(max(i.errorEstimate),0); tc.verifyEqual(i.evaluations,0);
    tc.verifyWarning(@() unresolved(K,t),'ContourRoots:ResponseUnresolved');
end
function unresolved(K,t)
    [~,~,~]=clsim(K,1e6*sin(20*t),t,'AbsTol',1e-20,'RelTol',1e-16);
end
function testGridOptionsAndFailedPreparation(tc)
    t=(0:.1:1).'; G=ndpair(1,[1 1]); K=ckernel(G,t);
    tc.verifyError(@() clsim(K,ones(size(t)),t*2),'ContourRoots:KernelGrid');
    tc.verifyError(@() clsim(K,ones(size(t)),t,'Interpolation','zoh'),'ContourRoots:KernelInterpolation');
    for key={'Method','SingularityBound','Feedthrough','MaxRefinements'}
        tc.verifyError(@() clsim(K,ones(size(t)),t,key{1},1),'ContourRoots:KernelOption');
    end
    tc.verifyError(@() clsim(K,ones(size(t)),t,'AbsTol'),'ContourRoots:KernelOption');
    tc.verifyError(@() ckernel(G,t,'MaxPoints',8,'Warn',false),'ContourRoots:KernelUnresolved');
    tc.verifyError(@() ckernel(G,t,'Plot',true),'ContourRoots:KernelOption');
    tc.verifyError(@() ckernel(G,t,'MaxMemoryMB',.001),'ContourRoots:ResponseBudget');
    tc.verifyError(@() clsim(K,ones(size(t)),t,'MaxPoints',8),'ContourRoots:ResponseBudget');
    [y,tt,i]=clsim(K,ones(size(t)),t.','Interpolation','FOH');
    tc.verifyEqual(tt,t); tc.verifySize(y,size(t)); tc.verifyTrue(i.converged);
end
function testReadOnlyAndSerialization(tc)
    t=(0:.1:1).'; K=ckernel(ndpair(1,[1 1]),t);
    info=K.Info; info.evaluations=-1;
    tc.verifyGreaterThan(K.Info.evaluations,0);
    failed=false;
    try, K.Time=t*2; catch, failed=true; end
    tc.verifyTrue(failed); tc.verifyEqual(K.Time,t);
    folder=tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
    file=fullfile(folder.Folder,'kernel.mat'); save(file,'K'); loaded=load(file,'K');
    u=sin(t); [a,~,ia]=clsim(K,u,t); [b,~,ib]=clsim(loaded.K,u,t);
    tc.verifyEqual(a,b); tc.verifyEqual(ia,ib);
end
function testPlotting(tc)
    t=(0:.1:1).'; K=ckernel(2,t,'Interpolation','zoh');
    fig=figure('Visible','off'); tc.addTeardown(@() close(fig)); ax=axes(fig);
    clsim(K,t,t,'Parent',ax);
    tc.verifyEqual(numel(findobj(ax,'Type','line')),2); tc.verifyFalse(ishold(ax));
end
function testOptionalExternalDelay(tc)
    tc.assumeTrue(~isempty(ver('control')),'Control System Toolbox not available.');
    t=(0:.05:2).'; sys=tf([1 2],[1 1],'InputDelay',.35); u=1+t;
    for hold={'foh','zoh'}
        K=ckernel(sys,t,'Interpolation',hold{1});
        [a,~,ia]=clsim(K,u,t); [b,~,ib]=clsim(sys,u,t,'Interpolation',hold{1});
        tc.verifyTrue(ia.converged&&ib.converged);
        tc.verifyLessThan(max(abs(a-b)),2e-6);
        tc.verifyEqual(K.Info.delay,.35); tc.verifyEqual(K.Info.feedthrough,1);
    end
end
function testExactDelayedReference(tc)
    % Independent reference: the delay-free lsim response shifted by the
    % delay (a whole number of samples). lsim itself is not used with the
    % delay because its FOH ramps the initial jump u(0) over one sample.
    tc.assumeTrue(~isempty(ver('control')),'Control System Toolbox not available.');
    t=(0:.01:5).'; n=30; s0=tf([2 1],[1 .6 4]); sys=s0; sys.InputDelay=n*.01;
    inputs={1+.5*t,sin(3*t),double(t>1.234)};
    for hold={'foh','zoh'}
        K=ckernel(sys,t,'Interpolation',hold{1});
        for k=1:numel(inputs)
            y0=lsim(s0,inputs{k},t,hold{1}); exact=[zeros(n,1);y0(1:end-n)];
            [y,~,i]=clsim(K,inputs{k},t);
            tc.verifyTrue(i.converged);
            tc.verifyLessThan(max(abs(y-exact)),5e-8);
        end
    end
end
function testSuperposition(tc)
    t=(0:.02:4).'; K=ckernel(@(s) 1./(s+1+.5*exp(-s)),t,'SingularityBound',0);
    u1=sin(2*t); u2=1+cos(5*t);
    y1=clsim(K,u1,t); y2=clsim(K,u2,t); y3=clsim(K,3*u1-2*u2,t,'Warn',false);
    tc.verifyLessThan(max(abs(y3-(3*y1-2*y2))),1e-12);
end
