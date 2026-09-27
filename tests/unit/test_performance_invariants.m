function tests=test_performance_invariants, tests=functiontests(localfunctions); end
function setupOnce(tc)
    root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
    tc.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,'matlab')));
    tc.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,'tools')));
end
function testMatrixContourReusesOnlyNestedNodes(tc)
    seen=[]; derivativeCalls=0;
    [r,i]=cmodes(@H,[-1.7 2.3 -.4 .9],'AssumeAnalytic',true,'Derivative',@D);
    tc.verifyEmpty(r); tc.verifyTrue(i.complete&&i.traceCheck.ok);
    tc.verifyEqual(i.factorizations,i.history(1).samples);
    tc.verifyEqual(numel(seen),1+i.factorizations); % reference plus contour
    tc.verifyEqual(derivativeCalls,i.factorizations);
    tc.verifyEqual(i.evaluations,numel(seen)+derivativeCalls);
    tc.verifyEqual(i.contourReuses,144); % 48 and 96 retained nodes
    tc.verifyEqual(numel(unique(seen(2:end))),i.factorizations);
    function A=H(s), seen(end+1)=s; A=[1 2;0 3]; end
    function A=D(~), derivativeCalls=derivativeCalls+1; A=zeros(2); end
end
function testScalarContourReusesWithoutGlobalCache(tc)
    seen=[]; value=1;
    [r,i]=croots(@F,[-1.7 2.3 -.4 .9],'AssumeAnalytic',true,'ContourPoints',12);
    tc.verifyTrue(i.complete); tc.verifyEmpty(r);
    tc.verifyEqual(numel(seen),192); tc.verifyEqual(numel(unique(seen)),192);
    seen=[]; value=2;
    [~,i]=croots(@F,[-1.7 2.3 -.4 .9],'AssumeAnalytic',true,'ContourPoints',12);
    tc.verifyTrue(i.complete); tc.verifyEqual(numel(seen),192);
    function y=F(s), seen=[seen s(:).']; y=value+zeros(size(s)); end
end
function testOddEvenPermutationSignsAndTrace(tc)
    % Dense row permutations, odd/even parity; analytic complex mixing.
    permutations={[2 1 3],[2 3 1],[3 2 1]};
    for k=1:numel(permutations)
        P=eye(3); P=P(permutations{k},:); R=[1 .1i 0;0 1 .3;0 0 1];
        H=@(s) P*diag([s+1 s+2 s+3])*R;
        [r,i]=cmodes(H,[-3.5 -.5 -1 1],'AssumeAnalytic',true,'Derivative',@(s) P*R);
        tc.verifyTrue(i.complete&&i.traceCheck.ok); tc.verifyEqual(r,[-1;-2;-3],'AbsTol',1e-10);
    end
end
function testCacheEvictionRectangularAndFreshModels(tc)
    rate=1; calls=0;
    M=cmimo(@G,'Size',[3 2]); t=(0:.1:1).';
    [large,a]=ckernel(M,t,'SingularityBound',0); tc.verifyEqual(calls,a.evaluations);
    calls=0;
    [small,b]=ckernel(M,t,'SingularityBound',0,'MaxMemoryMB',16); tc.verifyEqual(calls,b.evaluations);
    U=[1+t sin(t)]; [y,~,iy]=clsim(large,U,t); [z,~,iz]=clsim(small,U,t);
    tc.verifyTrue(iy.converged&&iz.converged); tc.verifyEqual(y,z,'AbsTol',1e-10);
    tc.verifyGreaterThan(b.evaluations,a.evaluations); tc.verifyGreaterThan(a.cacheHits,0);
    old=calls; clsim(small,[t t],t); tc.verifyEqual(calls,old);
    rate=2; fresh=ckernel(M,t,'SingularityBound',0);
    tc.verifyGreaterThan(max(abs(clsim(fresh,U,t)-y),[],'all'),1e-3);
    function v=G(s)
        calls=calls+1; v=[1/(s+rate) 1/(s+2*rate);0 1/(s+3*rate);2/(s+4*rate) 0];
    end
end
function testBuildTimerDoesNotSwallowFailures(tc)
    tmp=tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture); called=0;
    info=time_build_step('success',@job,tmp.Folder);
    tc.verifyEqual(called,1); tc.verifyEqual(info.status,'passed');
    tc.verifyError(@() time_build_step('failure',@bad,tmp.Folder),'test:ExpectedFailure');
    log=jsondecode(fileread(fullfile(tmp.Folder,'failure.json')));
    tc.verifyEqual(log.status,'failed'); tc.verifyEqual(log.errorIdentifier,'test:ExpectedFailure');
    occupied=fullfile(tmp.Folder,'not-a-folder'); fid=fopen(occupied,'w'); fclose(fid);
    tc.verifyWarning(@() time_build_step('unwritable',@job,occupied),'time_build_step:Log');
    tc.verifyEqual(called,2); % logging failure did not skip the operation
    tc.verifyError(@() time_build_step('unwritableFailure',@bad,occupied),'test:ExpectedFailure');
    function job(), called=called+1; end
    function bad(), error('test:ExpectedFailure','Expected failure'); end
end
function testReleaseKeepsAllValidationDependencies(tc)
    root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
    old=pwd; cleanup=onCleanup(@() cd(old)); cd(root); plan=buildfile;
    tc.verifyEqual(sort(string(plan('release').Dependencies)),sort(["test" "docs" "manual"]));
    tc.verifyEqual(string(plan('manual').Dependencies),"examples");
    tc.verifyEmpty(plan('package').Dependencies);
end
