function tests=test_time_response
    tests=functiontests(localfunctions);
end
function setupOnce(tc)
    root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
    tc.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,'matlab')));
    tc.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,'examples','models')));
end
function testFirstOrderAllMethods(tc)
    t=(0:.05:3).'; G=ndpair(1,[1 1]);
    for method={'fft','dehoog','quadrature'}
        [g,~,ig]=cimpulse(G,t,'Method',method{1});
        [y,~,iy]=cstep(G,t,'Method',method{1});
        tc.verifyTrue(ig.converged&&iy.converged);
        tc.verifyLessThan(max(abs(g-exp(-t))),2e-6);
        tc.verifyLessThan(max(abs(y-(1-exp(-t)))),2e-6);
    end
end
function testUnstableAndIntegrator(tc)
    t=(0:.1:5).';
    for method={'fft','dehoog','quadrature'}
        [g,~,i]=cimpulse(ndpair(1,[1 -.3]),t,'Method',method{1});
        tc.verifyTrue(i.converged); tc.verifyLessThan(max(abs(g-exp(.3*t))),5e-5);
        [y,~,i]=cstep(ndpair(1,[1 0]),t,'Method',method{1});
        tc.verifyTrue(i.converged); tc.verifyLessThan(max(abs(y-t)),2e-6);
    end
end
function testOscillator(tc)
    t=(0:.05:5).'; z=.02; w=5; wd=w*sqrt(1-z^2);
    reference=exp(-z*w*t).*sin(wd*t)/wd;
    for method={'fft','dehoog'}
        [g,~,i]=cimpulse(ndpair(1,[1 2*z*w w^2]),t,'Method',method{1});
        tc.verifyTrue(i.converged); tc.verifyLessThan(max(abs(g-reference)),2e-5);
    end
end
function testDiffusion(tc)
    t=(0:.05:3).'; a=.7;
    exact=zeros(size(t)); exact(2:end)=erfc(a./(2*sqrt(t(2:end))));
    for method={'fft','dehoog','quadrature'}
        [y,~,i]=cstep(@(s) exp(-a*sqrt(s)),t,'Method',method{1},'SingularityBound',0);
        tc.verifyTrue(i.converged); tc.verifyLessThan(max(abs(y-exact)),3e-6);
    end
end
function testSingularImpulse(tc)
    t=logspace(-3,1,20).';
    [g,~,i]=cimpulse(@(s) 1./sqrt(s),t,'Method','dehoog','SingularityBound',0,'RegularImpulse',true);
    tc.verifyTrue(i.converged); tc.verifyLessThan(max(abs(g.*sqrt(pi*t)-1)),1e-6);
    [y,~,i]=cstep(@(s) 1./sqrt(s),[0;t],'Method','dehoog','SingularityBound',0);
    tc.verifyTrue(i.converged); tc.verifyLessThan(max(abs(y-2*sqrt([0;t]/pi))),1e-6);
end
function testZOHAnalyticalRecursion(tc)
    t=(0:.05:3).'; u=sin(2*t); dt=t(2); reference=zeros(size(t));
    for k=2:numel(t), reference(k)=exp(-dt)*reference(k-1)+(1-exp(-dt))*u(k-1); end
    [y,~,i]=clsim(ndpair(1,[1 1]),u,t,'Interpolation','zoh');
    tc.verifyTrue(i.converged); tc.verifyLessThan(max(abs(y-reference)),2e-6);
end
function testFOHRampAndFeedthrough(tc)
    t=(0:.05:3).'; reference=t-1+exp(-t);
    for method={'fft','dehoog','quadrature'}
        [y,~,i]=clsim(ndpair(1,[1 1]),t,t,'Method',method{1});
        tc.verifyTrue(i.converged); tc.verifyLessThan(max(abs(y-reference)),2e-6);
        [y,~,i]=clsim(ndpair([2 3],[1 1]),ones(size(t)),t,'Method',method{1});
        tc.verifyTrue(i.converged); tc.verifyLessThan(max(abs(y-(3-exp(-t)))),2e-6);
    end
end
function testLinearityAndNoFutureLeakage(tc)
    t=(0:.1:3).'; u=sin(t); v=cos(t); G=ndpair(1,[1 1]);
    yu=clsim(G,u,t); yv=clsim(G,v,t); y=clsim(G,2*u-3*v,t);
    tc.verifyLessThan(max(abs(y-(2*yu-3*yv))),3e-6);
    changed=u; changed(t>1.5)=100;
    yc=clsim(G,changed,t);
    tc.verifyLessThan(max(abs(yc(t<=1.5)-yu(t<=1.5))),3e-6);
end
function testDelayAgainstAnalyticalStep(tc)
    t=(0:.05:3).'; T=.713;
    exact=(t>=T).*(1-exp(-max(0,t-T)));
    [y,~,i]=cstep(@(s) exp(-T*s)./(s+1),t,'SingularityBound',0,'Method','fft');
    tc.verifyTrue(i.converged); tc.verifyLessThan(max(abs(y-exact)),5e-6);
end
function testLTIExternalDelay(tc)
    tc.assumeTrue(~isempty(ver('control')),'Control System Toolbox not available.');
    t=(0:.05:3).'; T=.713; G=tf(1,[1 1],'InputDelay',T);
    [g,~,i]=cimpulse(G,t); ref=(t>=T).*exp(-max(0,t-T));
    tc.verifyTrue(i.converged); tc.verifyEqual(g,ref,'AbsTol',2e-6);
    pure=tf(2,1,'InputDelay',T); u=t;
    [y,~,i]=clsim(pure,u,t); tc.verifyTrue(i.converged);
    tc.verifyEqual(y,2*max(0,t-T),'AbsTol',1e-12);
    [g,~,i]=cimpulse(pure,t); tc.verifyEqual(g,zeros(size(t)));
    tc.verifyEqual(i.singularTerms.time,T); tc.verifyEqual(i.singularTerms.weight,2);
end
function testAgainstMATLABHolds(tc)
    tc.assumeTrue(~isempty(ver('control')),'Control System Toolbox not available.');
    t=(0:.04:4).'; u=sin(3*t)+.2*cos(5*t); G=tf([1 2],[1 .8 3]);
    for hold={'foh','zoh'}
        [y,~,i]=clsim(G,u,t,'Interpolation',hold{1});
        ref=lsim(G,u,t,[],hold{1});
        tc.verifyTrue(i.converged); tc.verifyLessThan(max(abs(y-ref)),5e-6);
    end
end
function testRefinementBudgetReported(tc)
    [~,~,i]=cstep(ndpair(1,[1 .02 100]),(0:.1:4).','MaxRefinements',1,'Warn',false);
    tc.verifyFalse(i.converged); tc.verifyFalse(i.certified);
    tc.verifyNotEmpty(i.history);
end
function testNumericalScale(tc)
    t=(0:.1:2).';
    [reference,~,ir]=cstep(ndpair(1,[1 1]),t,'AbsTol',1e-8);
    [small,~,is]=cstep(ndpair(1e-12,[1 1]),t,'AbsTol',1e-20);
    tc.verifyTrue(ir.converged&&is.converged);
    tc.verifyLessThan(max(abs(small/1e-12-reference)),1e-6);
end
function testDuhamelQuadrature(tc)
    errors=zeros(2,1);
    for j=1:2
        dt=.02/2^(j-1); t=(0:dt:2).'; u=sin(t); g=exp(-t);
        raw=dt*conv(g,u); direct=raw(1:numel(t))-dt/2*(g(1)*u+u(1)*g);
        y=clsim(ndpair(1,[1 1]),u,t);
        errors(j)=max(abs(y-direct));
    end
    tc.verifyLessThan(errors(2),errors(1)/3);
end
function testHighFrequencyLateTransient(tc)
    t=(.6:.05:1).'; w=75; z=.02; wd=w*sqrt(1-z^2);
    ref=(1-exp(-z*w*t).*(cos(wd*t)+z/sqrt(1-z^2)*sin(wd*t)))/w^2;
    [y,~,i]=cstep(ndpair(1,[1 2*z*w w^2]),t,'Method','quadrature', ...
        'SingularityBound',5,'AbsTol',1e-10,'RelTol',1e-5);
    tc.verifyTrue(i.converged); tc.verifyLessThan(max(abs(y-ref)),2e-9);
    tc.verifyGreaterThan(max(i.history(:,2)),w);
end
function testDeHoogOscillatoryGuard(tc)
    t=(.6:.05:1).'; w=75; z=.02; wd=w*sqrt(1-z^2);
    exact=(1-exp(-z*w*t).*(cos(wd*t)+z/sqrt(1-z^2)*sin(wd*t)))/w^2;
    [y,~,i]=cstep(ndpair(1,[1 2*z*w w^2]),t,'Method','dehoog', ...
        'SingularityBound',5,'AbsTol',1e-10,'RelTol',1e-5,'Warn',false);
    % A difficult double-precision case may be unresolved, never falsely
    % accepted on the low-degree plateau before the mode enters the band.
    tc.verifyLessThanOrEqual(abs(y(i.resolvedMask)-exact(i.resolvedMask)), ...
        1e-10+1e-5*abs(exact(i.resolvedMask)));
    tc.verifyGreaterThanOrEqual(max(i.history(:,2)),40);
end
function testNegativeExplicitImpulseLine(tc)
    t=(.05:.05:1).';
    [g,~,i]=cimpulse(@(s) 1./(s+20),t,'Abscissa',-10, ...
        'RegularImpulse',true,'InitialValue',1,'Method','quadrature');
    tc.verifyTrue(i.converged); tc.verifyLessThan(max(abs(g-exp(-20*t))),2e-6);
end
