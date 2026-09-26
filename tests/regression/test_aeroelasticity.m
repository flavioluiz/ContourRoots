function tests=test_aeroelasticity
    tests=functiontests(localfunctions);
end
function setupOnce(tc)
    root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
    import matlab.unittest.fixtures.PathFixture
    tc.applyFixture(PathFixture(fullfile(root,'matlab')));
    tc.applyFixture(PathFixture(fullfile(root,'examples','models')));
end
function testBesselAgainstHankel(tc)
    k=logspace(-5,2,100);
    h1=besselh(1,2,k); h0=besselh(0,2,k);
    tc.verifyLessThan(max(abs(theodorsen_laplace(1i*k)-h1./(h1+1i*h0))),2e-13);
    tc.verifyEqual(theodorsen_laplace(0),1);
    tc.verifyLessThan(abs(theodorsen_laplace(1e6)-.5),2e-7);
end
function testBranchAndConjugacy(tc)
    z=[-.4+.2i .2+1i 3+4i];
    tc.verifyEqual(theodorsen_laplace(conj(z)),conj(theodorsen_laplace(z)),'AbsTol',1e-13);
    tc.verifyError(@() theodorsen_laplace(-1),'theodorsen:BranchCut');
    m=aeroelastic_section;
    tc.verifyError(@() aeroelastic_delta(1i,100,m,'pk'),'aeroelastic:Nonanalytic');
end
function testLoadsIndependently(tc)
    m=aeroelastic_section('dlr'); b=m.b; a=m.a; rho=m.rho; U=175;
    for s=[-3+55i 4+75i]
        C=theodorsen_laplace(s*b/U); Q=zeros(2);
        for j=1:2
            q=zeros(2,1); q(j)=1; h=q(1); alpha=q(2);
            downwash=s*h+U*alpha+b*(.5-a)*s*alpha;
            lift=pi*rho*b^2*(s^2*h+U*s*alpha-a*b*s^2*alpha) ...
                +2*pi*rho*U*b*C*downwash;
            moment=pi*rho*b^2*(a*b*s^2*h-U*b*(.5-a)*s*alpha ...
                -b^2*(1/8+a^2)*s^2*alpha)+2*pi*rho*U*b^2*(.5+a)*C*downwash;
            Q(:,j)=[-lift;moment];
        end
        actual=s^2*m.M+m.K-aeroelastic_matrix(s,U,m);
        tc.verifyLessThan(norm(actual-Q,'fro')/norm(Q,'fro'),2e-14);
    end
end
function testNASAFlutter(tc)
    m=aeroelastic_section('nasa'); f=aeroelastic_flutter(m);
    tc.verifyLessThan(abs(f.U-m.reference.U),m.reference.velocityTolerance);
    tc.verifyLessThan(abs(f.k-m.reference.k),m.reference.frequencyTolerance);
    tc.verifyGreaterThan(f.crossingSpeed,0);
    tc.verifyLessThan(f.residual,1e-9);
    for factor=[.98 1.02]
        [r,i]=croots(@(s) aeroelastic_delta(s,f.U*factor,m), ...
            [-100 50 .5 180],'AssumeAnalytic',true);
        tc.verifyTrue(i.complete); tc.verifyEqual(i.count,2);
        tc.verifyEqual(sign(max(real(r))),sign(factor-1));
    end
end
function testDLRFlutter(tc)
    m=aeroelastic_section('dlr'); f=aeroelastic_flutter(m);
    tc.verifyLessThan(abs(f.U-m.reference.U),m.reference.velocityTolerance);
    pk=aeroelastic_flutter(m,'pk');
    tc.verifyLessThan(abs(f.U-pk.U),1e-6);
end
function testZeroSpeedAddedMass(tc)
    m=aeroelastic_section; b=m.b; a=m.a;
    Ma=pi*m.rho*b^2*[1 -a*b;-a*b (1/8+a^2)*b^2];
    reference=sort(sqrt(eig(m.K,m.M+Ma)));
    [r,i]=croots(@(s) aeroelastic_delta(s,0,m),[-1 1 1 160],'AssumeAnalytic',true);
    tc.verifyTrue(i.complete);
    tc.verifyEqual(sort(imag(r)),reference,'AbsTol',1e-8);
    dry=m; dry.rho=0;
    tc.verifyEqual(aeroelastic_matrix(2i,0,dry),-4*m.M+m.K,'AbsTol',1e-9);
end
function testContoursAndModalResidual(tc)
    m=aeroelastic_section; U=180;
    [r,i]=croots(@(s) aeroelastic_delta(s,U,m),[-100 50 .5 180],'AssumeAnalytic',true);
    [rr,ii]=croots(@(s) aeroelastic_delta(s,U,m),[-140 70 .1 230], ...
        'AssumeAnalytic',true,'ContourPoints',64,'RootTolerance',1e-8);
    tc.verifyTrue(i.complete && ii.complete); tc.verifyEqual(numel(rr),numel(r));
    for s=r(:).'
        tc.verifyLessThan(min(abs(rr-s)),1e-7);
        D=aeroelastic_matrix(s,U,m); sigma=svd(D);
        tc.verifyLessThan(sigma(end)/sigma(1),1e-9);
    end
    [pk,res]=aeroelastic_pk_roots(U,m,r);
    tc.verifyLessThan(max(res),1e-8);
    tc.verifyGreaterThan(max(abs(pk-r)),1e-3);
end
function testJonesStateSpace(tc)
    m=aeroelastic_section; U=180;
    r=aeroelastic_rational_roots(U,m,'jones'); upper=r(imag(r)>1e-6);
    tc.verifyEqual(numel(r),6); tc.verifyEqual(numel(upper),2);
    [rc,info]=croots(@(s) aeroelastic_delta(s,U,m,'jones'),[-100 50 20 150],'AssumeAnalytic',true);
    tc.verifyTrue(info.complete); tc.verifyEqual(info.count,2);
    for s=upper(:).'
        tc.verifyLessThan(min(abs(rc-s)),1e-7);
        tc.verifyLessThan(abs(aeroelastic_delta(s,U,m,'jones')),1e-10);
    end
end
function testQuasiSteadyRealization(tc)
    m=aeroelastic_section;
    r=aeroelastic_rational_roots(150,m,'quasisteady');
    tc.verifyEqual(numel(r),4);
    for s=r(:).'
        tc.verifyLessThan(abs(aeroelastic_delta(s,150,m,'quasisteady')),1e-10);
    end
    m=aeroelastic_section('dlr');
    tc.verifyGreaterThan(max(real(aeroelastic_rational_roots(.01,m,'quasisteady'))),0);
end
function testIndependentHarmonicFlutter(tc)
    for kind={'nasa','dlr'}
        m=aeroelastic_section(kind{1}); f=aeroelastic_flutter(m);
        independent=aeroelastic_harmonic_flutter(m);
        tc.verifyEqual(independent.U,f.U,'AbsTol',1e-6);
        tc.verifyEqual(independent.k,f.k,'AbsTol',1e-8);
    end
end
function testRightHalfPlaneStabilityCount(tc)
    % The branch cut is on the NEGATIVE real axis: a right-half-plane window
    % is admissible and counts every unstable root, including real ones.
    m=aeroelastic_section; rhp=[1e-3 60 -250 250];
    [r,i]=croots(@(s) aeroelastic_delta(s,170,m),rhp,'AssumeAnalytic',true);
    tc.verifyTrue(i.complete); tc.verifyEmpty(r);
    [r,i]=croots(@(s) aeroelastic_delta(s,174,m),rhp,'AssumeAnalytic',true);
    tc.verifyTrue(i.complete); tc.verifyEqual(numel(r),2);   % crashed in 0.2.0
    Ud=sqrt(m.K(2,2)/(2*pi*m.rho*m.b^2*(.5+m.a)));
    [r,i]=croots(@(s) aeroelastic_delta(s,1.05*Ud,m),[1e-3 200 -300 300], ...
        'AssumeAnalytic',true,'ContourRefinements',10);
    tc.verifyTrue(i.complete);
    tc.verifyTrue(any(abs(imag(r))<1e-8 & real(r)>0));      % divergence: real root
end
function testFlutterBySpectralAbscissa(tc)
    m=aeroelastic_section;
    alpha=@(U) max(real(croots(@(s) aeroelastic_delta(s,U,m),[-100 50 .5 180],'AssumeAnalytic',true)));
    tc.verifyEqual(fzero(alpha,[150 200]),aeroelastic_flutter(m).U,'AbsTol',1e-5);
end
