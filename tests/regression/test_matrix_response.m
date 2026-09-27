function tests=test_matrix_response, tests=functiontests(localfunctions); end
function setupOnce(tc)
    root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
    tc.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,'matlab')));
end
function testStepImpulseAndDirectTerms(tc)
    M=cmimo({ndpair(1,[1 1]),2;0,ndpair(2,[1 2])}); t=(0:.1:1).';
    [y,~,i]=cstep(M,t); tc.verifySize(y,[11 2 2]); tc.verifyTrue(i.converged);
    tc.verifyLessThan(max(abs(y(:,1,1)-(1-exp(-t)))),2e-6);
    tc.verifyEqual(y(:,1,2),2*ones(size(t))); tc.verifyEqual(y(:,2,1),zeros(size(t)));
    [h,~,ii]=cimpulse(M,t); tc.verifyTrue(ii.converged);
    tc.verifyLessThan(max(abs(h(:,1,1)-exp(-t))),2e-6);
    tc.verifyEqual(ii.singularTerms,struct('output',1,'input',2,'time',0,'order',0,'weight',2));
end
function testSuperpositionAndReuse(tc)
    t=(0:.1:1).'; U=[1+sin(t),cos(t)]; M=cmimo({ndpair(1,[1 1]),ndpair(2,[1 2]);0,1});
    for hold={'foh','zoh'}
        [K,ki]=ckernel(M,t,'Interpolation',hold{1});
        [y,~,i]=clsim(K,U,t); [fresh,~,j]=clsim(M,U,t,'Interpolation',hold{1});
        ref=clsim(ndpair(1,[1 1]),U(:,1),t,'Interpolation',hold{1})+clsim(ndpair(2,[1 2]),U(:,2),t,'Interpolation',hold{1});
        tc.verifyLessThan(max(abs(y(:,1)-ref)),3e-6); tc.verifyEqual(y(:,2),U(:,2));
        tc.verifyLessThan(max(abs(y(:)-fresh(:))),3e-6); tc.verifyTrue(i.converged&&j.converged);
        tc.verifyEqual(i.evaluations,0); tc.verifyEqual(i.kernelPreparationEvaluations,ki.evaluations);
    end
end
function testSharedEvaluationAndFrozenBank(tc)
    t=(0:.1:1).'; calls=0; enabled=true; rate=1;
    M=cdyn(@H,eye(2),eye(2),zeros(2),'Dimensions',[2 2 2]);
    [K,info]=ckernel(M,t,'SingularityBound',0); old=calls;
    tc.verifyEqual(calls,info.factorizations); tc.verifyGreaterThan(info.cacheHits,0);
    enabled=false; rate=2;
    [a,~,ia]=clsim(K,[ones(size(t)) zeros(size(t))],t);
    [b,~,ib]=clsim(K,[zeros(size(t)) ones(size(t))],t);
    tc.verifyEqual(calls,old); tc.verifyTrue(ia.converged&&ib.converged);
    tc.verifyLessThan(max(abs(a(:,1)-(1-exp(-t)))),2e-6);
    tc.verifyLessThan(max(abs(b(:,2)-(1-exp(-2*t))/2)),2e-6);
    tc.verifyEqual(ia.factorizations,0); tc.verifyEqual(ib.evaluations,0);
    folder=tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture); file=fullfile(folder.Folder,'bank.mat');
    save(file,'K'); loaded=load(file); tc.verifyEqual(clsim(loaded.K,[ones(size(t)) zeros(size(t))],t),a);
    function v=H(s), assert(enabled); calls=calls+1; v=diag([s+rate s+2*rate]); end
end
function testCancellationAndContracts(tc)
    t=(0:.1:1).'; M=cmimo({ndpair(1,[1 1]),ndpair(1,[1 1])}); K=ckernel(M,t);
    [y,~,i]=clsim(K,[ones(size(t)) -ones(size(t))],t,'AbsTol',1e-18,'RelTol',1e-16,'Warn',false);
    tc.verifyLessThan(max(abs(y)),1e-12); tc.verifyFalse(i.converged); tc.verifyGreaterThan(max(i.errorEstimate),0);
    tc.verifyError(@() clsim(K,zeros(numel(t),2),t*2),'ContourRoots:KernelGrid');
    tc.verifyError(@() clsim(K,zeros(numel(t),2),t,'Method','fft'),'ContourRoots:KernelOption');
    tc.verifyError(@() clsim(K,zeros(numel(t),2),t,'Interpolation','zoh'),'ContourRoots:KernelInterpolation');
    tc.verifyError(@() cstep(M,t,'MaxMemoryMB',.001),'ContourRoots:MatrixBudget');
    tc.verifyError(@() cstep(cmimo(@(s) 1/(s+1),'Size',[1 1]),t),'ContourRoots:ResponseDomain');
end
function testDelaysAndOpaqueMetadata(tc)
    t=(0:.05:1).'; M=cmimo([2 3],'Delay',[.2 .4]);
    [y,~,i]=cstep(M,t); tc.verifyTrue(i.converged);
    tc.verifyEqual(y(:,1,1),2*(t>=.2)); tc.verifyEqual(y(:,1,2),3*(t>=.4));
    [~,~,i]=cimpulse(M,t); tc.verifyEqual([i.singularTerms.time],[.2 .4]);
    K=ckernel(M,t); U=[1+t 1+2*t]; y=clsim(K,U,t);
    ref=2*(1+t-.2).*(t>=.2)+3*(1+2*(t-.4)).*(t>=.4);
    tc.verifyLessThan(max(abs(y-ref)),1e-12);
    A=cmimo(@(s) [1/(s+1) 2],'Size',[1 2],'Feedthrough',[0 2],'InitialValue',[1 0]);
    [y,~,i]=cimpulse(A,t,'SingularityBound',0,'RegularImpulse',true);
    tc.verifyTrue(i.converged); tc.verifyLessThan(max(abs(y(:,1,1)-exp(-t))),2e-6);
end
function testBranchUnstableAndSingleton(tc)
    t=(0:.1:1).'; M=cmimo({@(s) exp(-.7*sqrt(s)),ndpair(1,[1 -.3])});
    K=ckernel(M,t,'SingularityBound',.3); [y,~,i]=clsim(K,ones(numel(t),2),t);
    ref=clsim(@(s) exp(-.7*sqrt(s)),ones(size(t)),t,'SingularityBound',.3)+clsim(ndpair(1,[1 -.3]),ones(size(t)),t);
    tc.verifyTrue(i.converged); tc.verifyLessThan(max(abs(y-ref)),5e-6);
    [a,~,ia]=cstep(cmimo({ndpair(1,[1 1])}),t);
    tc.verifyEqual(a,cstep(ndpair(1,[1 1]),t)); tc.verifyTrue(ia.converged);
end
function testMethodsHandlesAndFailureContracts(tc)
    t=(0:.2:1).'; M=cmimo({ndpair(1,[1 1]);2});
    for method={'fft','dehoog','quadrature'}
        [K,~]=ckernel(M,t,'Method',method{1});
        [a,~,i]=clsim(K,@(q) 1+q,t,'AbsTol',[1e-6 1e-8]);
        [b,~,j]=clsim(M,1+t,t,'Method',method{1},'AbsTol',[1e-6 1e-8]);
        tc.verifyTrue(i.converged&&j.converged); tc.verifyLessThan(max(abs(a(:)-b(:))),2e-6);
    end
    tc.verifyError(@() ckernel(M,t,'MaxPoints',8),'ContourRoots:KernelUnresolved');
    tc.verifyError(@() cstep(M,t,'MaxEvaluations',1),'ContourRoots:MatrixBudget');
    tc.verifyError(@() cstep(cmimo(@(s) 1/(s+1),'Size',[1 1],'DomainCheck',@(s) false),t, ...
        'SingularityBound',0),'ContourRoots:MatrixDomain');
    failed=false;
    try, clsim(M,@(q) ones(1,numel(q)),t); catch, failed=true; end
    tc.verifyTrue(failed);
    K=ckernel(M,t); failed=false;
    try, K.Time=t*2; catch, failed=true; end
    tc.verifyTrue(failed);
end
function testOptionalMIMOOracle(tc)
    tc.assumeTrue(~isempty(ver('control')));
    t=(0:.05:2).'; sys=ss([-1 .2;0 -2],eye(2),[1 2;3 1],[.1 0;0 .2]);
    M=cmimo(sys); U=[sin(t) cos(2*t)];
    for hold={'foh','zoh'}
        [y,~,info]=clsim(M,U,t,'Interpolation',hold{1});
        ref=lsim(sys,U,t,hold{1}); tc.verifyTrue(info.converged);
        tc.verifyLessThan(max(abs(y(:)-ref(:))),3e-6);
    end
    delayed=tf({1,2},{[1 1],[1 2]},'IODelay',[.2 .4]); A=cmimo(delayed);
    K=ckernel(A,t); [y,~,i]=clsim(K,ones(numel(t),2),t);
    ref=(1-exp(-(t-.2))).*(t>=.2)+(1-exp(-2*(t-.4))).*(t>=.4);
    tc.verifyTrue(i.converged); tc.verifyLessThan(max(abs(y-ref)),2e-6);
    gain=tf(2,1,'InputDelay',.1); A=cmimo(gain,'Delay',.2);
    [v,~,iv]=cstep(A,[0;.3;.6]); tc.verifyTrue(iv.converged); tc.verifyEqual(v,[0;2;2]);
    K=ckernel(A,[0;.3;.6]); tc.verifyEqual(clsim(K,ones(3,1),[0;.3;.6]),[0;2;2]);
end
function testPlotAndPreparationBudgets(tc)
    t=(0:.1:1).'; M=cmimo([1 2;3 4]);
    fig=figure('Visible','off'); tc.addTeardown(@() close(fig)); ax=axes(fig);
    cstep(M,t,'Parent',ax); tc.verifyEqual(numel(findobj(ax,'Type','line')),4);
    K=ckernel(M,t); clsim(K,[t t],t,'Parent',ax);
    tc.verifyEqual(numel(findobj(ax,'Type','line')),2); tc.verifyFalse(ishold(ax));
    tc.verifyError(@() ckernel(M,t,'Plot',true),'ContourRoots:KernelOption');
    tc.verifyError(@() clsim(K,[t t],t,'MaxMemoryMB',1e-6),'ContourRoots:MatrixBudget');
end

function testChannelDelaysDefaultBudgetAgainstLsim(tc)
    % A 2x2 rational model with a different delay per channel, simulated
    % with the default (unlimited) shared budget; only the requested channel
    % is evaluated at each node. ZOH lsim is exact for delays = k*dt.
    tc.assumeTrue(~isempty(ver('control')),'Control System Toolbox not available.');
    s=tf('s'); t=(0:0.01:4).';
    G=[1/(s+1), 2/(s^2+0.4*s+4); (s+3)/(s+2), 0.5/(s+0.5)];
    G.IODelay=[0 0.3; 0.5 0.1]; u=[sin(2*t)+1, double(t>1.005)];
    [y,~,info]=clsim(cmimo(G),u,t,'Interpolation','zoh');
    tc.verifyTrue(info.converged);
    tc.verifyLessThan(max(abs(y-lsim(G,u,t,'zoh')),[],'all'),1e-7);
end
