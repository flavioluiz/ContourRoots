function tests=test_shared_grid, tests=functiontests(localfunctions); end
function setupOnce(tc)
    root=fileparts(fileparts(fileparts(mfilename('fullpath')))); tc.TestData.root=root;
    tc.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,'matlab')));
end
function testScalarAndCommonDiagnostics(tc)
    t=(0:.1:.5).'; G=ndpair(1,[1 1]); M=cmimo({G});
    for method={'fft','dehoog'}
        for hold={'foh','zoh'}
            args={'Method',method{1},'Interpolation',hold{1},'AbsTol',1e-7,'RelTol',1e-5};
            A=ckernel(G,t,args{:}); [B,b]=ckernel(M,t,args{:},'SharedGrid',true);
            tc.verifyTrue(b.sharedGrid&&b.converged); U=1+sin(2*t);
            tc.verifyEqual(clsim(A,U,t),clsim(B,U,t),'AbsTol',2e-6);
        end
    end
end
function testBlockSolveWorstChannelAndReuse(tc)
    t=(0:.1:.6).'; seen=[]; enabled=true;
    B=[1 .5;-.2 1;1 -1]; C=[1 0 0;0 .2 1;2 1 -.3]; D=[.1 0;0 .2;-.1 .3];
    M=cdyn(@H,B,C,D,'Dimensions',[3 3 2],'Feedthrough',D,'Delay',[0 .1;.2 0;.15 .3]);
    U=[1+sin(t) cos(2*t)]; opts={'SingularityBound',0,'AbsTol',[1e-6;2e-7;3e-7],'RelTol',1e-5};
    for method={'fft','dehoog'}
        [A,a]=ckernel(M,t,opts{:},'Method',method{1},'SharedGrid',false); seen=[];
        [K,b]=ckernel(M,t,opts{:},'Method',method{1},'SharedGrid',true);
        tc.verifyEqual(b.factorizations,numel(seen)); tc.verifyEqual(b.factorizations,numel(unique(seen)));
        tc.verifyEqual(b.rhsColumns,2*b.factorizations); tc.verifyTrue(b.converged&&a.converged);
        tc.verifyEqual(b.AbsTol,[1e-6;2e-7;3e-7]);
        enabled=false;
        [y,~,info]=clsim(K,U,t,'AbsTol',1e-4); z=clsim(A,U,t,'AbsTol',1e-4);
        tc.verifyTrue(info.converged); tc.verifyEqual(info.factorizations,0);
        tc.verifyEqual(y,z,'AbsTol',5e-6); enabled=true;
        % Every active channel and both orders used the same refinement grids.
        for c=2:numel(b.channelInfo)
            tc.verifyEqual(b.channelInfo{c}.step.points,b.channelInfo{1}.step.points);
            tc.verifyEqual(b.channelInfo{c}.ramp.points,b.channelInfo{1}.ramp.points);
        end
    end
    function K=H(s), assert(enabled); seen(end+1)=s; K=[s+1 .2 0;0 s+2 .1;0 0 s+3]; end
end
function testChannelsBoundsAndExactTerms(tc)
    t=(0:.1:.6).'; M=cmimo({ndpair(1,[1 -.2]),2;0,ndpair([1 3],[1 2])}, ...
        'Delay',[.1 .2;.3 .15]); U=[1+t cos(t)];
    for method={'fft','dehoog'}
        for hold={'foh','zoh'}
            args={'Method',method{1},'Interpolation',hold{1},'AbsTol',[1e-7 2e-7],'RelTol',1e-5};
            A=ckernel(M,t,args{:}); K=ckernel(M,t,args{:},'SharedGrid',true);
            tc.verifyEqual(clsim(K,U,t),clsim(A,U,t),'AbsTol',4e-6);
            tc.verifyEqual(K.Info.channelInfo{1}.delay,.1);
            tc.verifyEqual(K.Info.channelInfo{4}.feedthrough,1);
        end
    end
    [K,k]=ckernel(cmimo([2 0;0 3],'Delay',[.2 0;0 1]),t,'SharedGrid',true);
    tc.verifyEqual(k.evaluations,0); tc.verifyEqual(clsim(K,[1+t t],t),[2*(1+t-.2).*(t>=.2) 0*t],'AbsTol',1e-14);
end
function testBatchedBudgetsAndFailures(tc)
    t=(0:.1:.4).'; weights=[1 2;3 4];
    A=cmimo(@(s) weights/(s+1),'Size',[2 2]);
    B=cmimo(@(s) weights.*reshape(1./(s+1),1,1,[]),'Size',[2 2],'Batched',true);
    opts={'SingularityBound',0,'AbsTol',1e-6,'RelTol',1e-4,'SharedGrid',true};
    [a,ia]=ckernel(A,t,opts{:}); [b,ib]=ckernel(B,t,opts{:});
    tc.verifyEqual(clsim(a,[1+t t],t),clsim(b,[1+t t],t),'AbsTol',1e-12);
    tc.verifyEqual(ia.evaluations,ib.evaluations); tc.verifyLessThan(ib.factorEvaluations,ia.factorEvaluations);
    [limited,li]=ckernel(A,t,opts{:},'MaxEvaluations',ia.evaluations);
    tc.verifyEqual(li.evaluations,ia.evaluations); tc.verifyEqual(limited.Info.converged,true);
    tc.verifyError(@() ckernel(A,t,opts{:},'MaxEvaluations',ia.evaluations-1),'ContourRoots:MatrixBudget');
    tc.verifyError(@() ckernel(A,t,opts{:},'MaxPoints',8),'ContourRoots:KernelUnresolved');
    tc.verifyError(@() ckernel(A,t,opts{:},'MaxMemoryMB',.01),'ContourRoots:MatrixBudget');
    tc.verifyError(@() ckernel(A,t,opts{:},'Method','quadrature'),'ContourRoots:KernelOption');
    tc.verifyError(@() ckernel(A,t,opts{:},'SharedGrid',2),'ContourRoots:KernelOption');
    tc.verifyError(@() ckernel(A,t,opts{:},'Abscissa',0),'ContourRoots:ResponseDomain');
    bad=cmimo(@(s) weights/(s+1),'Size',[2 2],'DomainCheck',@(s) false);
    tc.verifyError(@() ckernel(bad,t,opts{:}),'ContourRoots:MatrixDomain');
    complexModel=cmimo(@(s) 1i/(s+1),'Size',[1 1]);
    tc.verifyError(@() ckernel(complexModel,t,opts{:}),'ContourRoots:ResponseSymmetry');
end
function testOldSnapshotsAndNewSerialization(tc)
    saved=load(fullfile(tc.TestData.root,'tests','fixtures','kernel_v080.mat'));
    t=saved.t; U=[ones(size(t)) t];
    tc.verifyEqual(clsim(saved.scalar,U(:,1),t),1-exp(-t),'AbsTol',2e-6);
    tc.verifyEqual(clsim(saved.matrix,U,t),[1-exp(-t)+2*t t-(1-exp(-2*t))/2],'AbsTol',2e-6);
    M=cmimo({ndpair(1,[1 1]),2;0,ndpair(2,[1 2])});
    K=ckernel(M,t,'SharedGrid',true);
    tmp=tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture); file=fullfile(tmp.Folder,'shared.mat');
    save(file,'K'); loaded=load(file);
    tc.verifyEqual(clsim(loaded.K,U,t),clsim(K,U,t));
    tc.verifyEqual(clsim(loaded.K,U,t),clsim(saved.matrix,U,t),'AbsTol',2e-6);
end
function testOptionalLTIDelays(tc)
    tc.assumeTrue(~isempty(ver('control')));
    t=(0:.1:.6).'; G=tf({1,2},{[1 1],[1 2]},'IODelay',[.2 .3]);
    M=cmimo(G); args={'AbsTol',1e-7,'RelTol',1e-5}; U=[1+t sin(t)];
    A=ckernel(M,t,args{:}); B=ckernel(M,t,args{:},'SharedGrid',true);
    tc.verifyEqual(clsim(A,U,t),clsim(B,U,t),'AbsTol',2e-6);
end
function testOscillatoryWorstChannelAndEviction(tc)
    t=(0:.1:1).'; M=cmimo({ndpair(1,[1 1]);ndpair(100,[1 .4 100])});
    w=sqrt(100-.2^2); ref=[1-exp(-t),1-exp(-.2*t).*(cos(w*t)+.2/w*sin(w*t))];
    for method={'fft','dehoog'}
        [K,i]=ckernel(M,t,'SharedGrid',true,'Method',method{1},'AbsTol',[1e-6 1e-7],'RelTol',1e-5);
        [y,~,out]=clsim(K,ones(size(t)),t,'AbsTol',1e-5);
        tc.verifyTrue(i.converged&&out.converged); tc.verifyEqual(y,ref,'AbsTol',3e-6);
    end
    A=cmimo(@(s) [1 2;3 4]/(s+1),'Size',[2 2]);
    opts={'SharedGrid',true,'SingularityBound',0,'AbsTol',1e-6,'RelTol',1e-4};
    [large,a]=ckernel(A,t,opts{:}); [small,b]=ckernel(A,t,opts{:},'MaxMemoryMB',8);
    tc.verifyEqual(clsim(large,[1+t t],t),clsim(small,[1+t t],t),'AbsTol',1e-12);
    tc.verifyGreaterThan(b.evaluations,a.evaluations);
end
