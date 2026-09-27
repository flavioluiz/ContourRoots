function tests=test_matrix_modes, tests=functiontests(localfunctions); end
function setupOnce(tc)
    root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
    tc.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,'matlab')));
end
function testSimpleAndDeterminant(tc)
    H=@(s) [s+1 .2*exp(-s);0 s+2]; box=[-3 0 -1 1];
    [r,i]=cmodes(H,box,'AssumeAnalytic',true,'Derivative',@(s) [1 -.2*exp(-s);0 1]);
    [d,j]=croots(@(s) det(H(s)),box,'AssumeAnalytic',true);
    tc.verifyTrue(i.complete&&j.complete); tc.verifyEqual(r,d,'AbsTol',1e-10);
    tc.verifyTrue(i.traceCheck.ok); tc.verifyEqual(i.count,2);
    tc.verifyLessThan(max(i.residuals),1e-10);
    for k=1:numel(r)
        tc.verifyLessThan(norm(H(r(k))*i.rightVectors{k}),1e-10);
        tc.verifyLessThan(norm(i.leftVectors{k}'*H(r(k))),1e-10);
    end
end
function testScalingAndPivotParity(tc)
    for magnitude=[1e-200 1e200]
        H=@(s) magnitude*[0 s+2;s+1 0];
        [r,i]=cmodes(H,[-3 0 -1 1],'AssumeAnalytic',true);
        tc.verifyTrue(i.complete); tc.verifyEqual(r,[-1;-2],'AbsTol',1e-10);
        tc.verifyEqual(i.count,2); tc.verifyLessThan(max(i.residuals),1e-10);
        tc.verifyTrue(det(H(.5))==0||~isfinite(det(H(.5))));
    end
    H=@(s) diag([1e-120*(s+1),1e120*(s+2)]);
    [r,i]=cmodes(H,[-3 0 -1 1],'AssumeAnalytic',true);
    tc.verifyTrue(i.complete); tc.verifyEqual(r,[-1;-2],'AbsTol',1e-10);
end
function testMoreModesThanDimension(tc)
    [r,i]=cmodes(@(s) sinh(s),[-1 1 -10 10],'AssumeAnalytic',true);
    tc.verifyTrue(i.complete); tc.verifyEqual(i.count,7);
    tc.verifyEqual(sort(imag(r)),(-3:3)'*pi,'AbsTol',1e-10);
end
function testRepeatedAndDefective(tc)
    [r,i]=cmodes(@(s) (s+1)*eye(2),[-2 0 -1 1],'AssumeAnalytic',true);
    tc.verifyTrue(i.complete); tc.verifyEqual(r,-1,'AbsTol',1e-12);
    tc.verifyEqual(i.multiplicity,2); tc.verifyEqual(i.nullity,2);
    [r,i]=cmodes(@(s) [s+1 1;0 s+1],[-2 0 -1 1],'AssumeAnalytic',true);
    tc.verifyFalse(i.complete); tc.verifyTrue(i.countComplete);
    tc.verifyEqual(i.count,2); tc.verifyEqual(r,-1,'AbsTol',1e-7);
    tc.verifyEqual(i.localCounts,2); tc.verifyEqual(i.nullity,1);
    tc.verifyTrue(isnan(i.multiplicity)); tc.verifyNotEmpty(i.clusters);
end
function testCloseClusterNotRepeated(tc)
    H=@(s) diag([s+1,s+1+1e-9]);
    [~,i]=cmodes(H,[-2 0 -1 1],'AssumeAnalytic',true);
    tc.verifyTrue(i.countComplete); tc.verifyEqual(i.count,2);
    tc.verifyFalse(i.multiplicityComplete); tc.verifyNotEmpty(i.clusters);
end
function testStructuredHiddenModeAndGuard(tc)
    calls=0;
    M=cdyn(@H,[1;0],[1 0],0,'Dimensions',[2 1 1],'DomainCheck',@(s) real(s)<=0);
    [r,i]=cmodes(M,[-3 0 -1 1],'AssumeAnalytic',true);
    tc.verifyTrue(i.complete); tc.verifyEqual(r,[-1;-2],'AbsTol',1e-10);
    tc.verifyEqual(calls,i.evaluations);
    tc.verifyError(@() cmodes(M,[-3 1 -1 1],'AssumeAnalytic',true),'ContourRoots:MatrixDomain');
    function A=H(s), calls=calls+1; A=diag([s+1 s+2]); end
end
function testBoundariesBudgetsAndContracts(tc)
    [~,i]=cmodes(@(s) s,[-1 0 -1 1],'AssumeAnalytic',true);
    tc.verifyFalse(i.complete); tc.verifyNotEmpty(i.unresolvedBoxes);
    [~,i]=cmodes(@(s) s,[-1 1 -1 1],'AssumeAnalytic',true,'MaxEvaluations',10);
    tc.verifyFalse(i.complete); tc.verifyLessThanOrEqual(i.evaluations,10);
    [~,i]=cmodes(@(s) diag([s+1 s+2]),[-3 0 -1 1],'AssumeAnalytic',true,'MaxCells',1);
    tc.verifyFalse(i.complete); tc.verifyLessThanOrEqual(i.cells,1);
    [~,i]=cmodes(@(s) eye(30),[-1 1 -1 1],'AssumeAnalytic',true,'MaxMemoryMB',.01);
    tc.verifyFalse(i.complete); tc.verifySubstring(i.stopReason,'MaxMemoryMB');
    [~,i]=cmodes(zeros(2),[-1 1 -1 1],'AssumeAnalytic',true);
    tc.verifyFalse(i.complete);
    [r,i]=cmodes(eye(2),[-1 1 -1 1],'AssumeAnalytic',true);
    tc.verifyTrue(i.complete); tc.verifyEmpty(r); tc.verifyEqual(i.count,0);
    tc.verifyError(@() cmodes(@(s) s,[-1 1 -1 1]),'ContourRoots:MatrixAnalytic');
    tc.verifyError(@() cmodes(cmimo(eye(2)),[-1 1 -1 1],'AssumeAnalytic',true),'ContourRoots:MatrixRepresentation');
    tc.verifyError(@() cmodes(@(s) ones(2,3),[-1 1 -1 1],'AssumeAnalytic',true),'ContourRoots:MatrixShape');
end
function testQuadraticMixingAndRankSensitivity(tc)
    M=diag([2 1]); C=[.3 .02;.02 .2]; K=[3 -1;-1 2];
    oracle=eig([zeros(2) eye(2);-M\K -M\C]);
    L=[1 2i;-.3 1]; R=[1 .4;-.2i 2];
    H=@(s) L*(s^2*M+s*C+K)*R;
    before=rng;
    for tolerance=[1e-9 1e-8 1e-7]
        [r,i]=cmodes(H,[-1 .2 -3 3],'AssumeAnalytic',true,'RankTolerance',tolerance);
        tc.verifyTrue(i.complete); verifySpectrum(tc,r,oracle,1e-10);
        tc.verifyLessThan(max(i.unscaledResiduals),1e-10);
        tc.verifyLessThan(max(i.unscaledLeftResiduals),1e-10);
    end
    [r,i]=cmodes(H,[-1 .2 -3 3],'AssumeAnalytic',true,'Scaling','none');
    tc.verifyTrue(i.complete); verifySpectrum(tc,r,oracle,1e-10); tc.verifyEqual(rng,before);
end
function testDelayLambertReference(tc)
    % Independently generated in MATLAB Symbolic Math Toolbox from exact
    % -1+lambertw(k,-exp(sym(1))/2), k=-4:3, converted to double. No optional
    % toolbox is required to RUN this regression. Plus the hidden root -2.
    rr=[-3.7120063323966379;-3.3440066442276852;-2.7506884347870275;-1.1026594768180491];
    ii=[-20.287461160139785;-13.970937813330176;-7.628391593321969;-1.5025802096948044];
    oracle=[rr+1i*ii;rr-1i*ii;-2];
    H=@(s) diag([s+1+.5*exp(-s),s+2]);
    [r,i]=cmodes(H,[-6 1 -25 25],'AssumeAnalytic',true);
    tc.verifyTrue(i.complete); tc.verifyEqual(i.count,9); verifySpectrum(tc,r,oracle,1e-9);
end
function testGuardedDerivativeAndConservation(tc)
    box=[-.7 .9 -.8 .6]; visited=[];
    [r,i]=cmodes(@H,box,'AssumeAnalytic',true);
    tc.verifyTrue(i.complete); verifySpectrum(tc,r,[-.3;.4],1e-10);
    tc.verifyTrue(all(real(visited)>=box(1)&real(visited)<=box(2)));
    tc.verifyTrue(all(imag(visited)>=box(3)&imag(visited)<=box(4)));
    [~,j]=cmodes(@(s) diag([s+.3,s-.4]),box,'AssumeAnalytic',true,'MaxDepth',0);
    tc.verifyTrue(j.countComplete); tc.verifyFalse(j.complete); tc.verifyNotEmpty(j.unresolvedBoxes);
    [~,j]=cmodes(@(s) s,[-1 1 -1 1],'AssumeAnalytic',true,'Derivative',@(s) 0,'ContourRefinements',3);
    tc.verifyFalse(j.complete); tc.verifyFalse(j.traceCheck.ok);
    tc.verifyError(@() cmodes(@bad,box,'AssumeAnalytic',true),'test:CallbackFailure');
    function A=H(s)
        assert(real(s)>=box(1)&&real(s)<=box(2)&&imag(s)>=box(3)&&imag(s)<=box(4));
        visited(end+1)=s; A=diag([s+.3,s-.4]);
    end
    function A=bad(~), error('test:CallbackFailure','User evaluator failure'); A=[]; end %#ok<UNRCH>
end
function verifySpectrum(tc,r,oracle,tol)
    tc.assertEqual(numel(r),numel(oracle));
    % One-to-one matching: remove each matched oracle, so duplicates cannot
    % mask a missing distinct mode in these well-separated fixtures.
    for z=r(:).'
        [distance,k]=min(abs(z-oracle)); tc.verifyLessThan(distance/(1+abs(oracle(k))),tol);
        oracle(k)=[];
    end
end
