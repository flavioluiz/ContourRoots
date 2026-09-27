function tests=test_matrix_modes_wing, tests=functiontests(localfunctions); end
function setupOnce(tc)
    repo=fileparts(fileparts(fileparts(mfilename('fullpath'))));
    tc.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(repo,'matlab')));
    tc.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(repo,'examples','models')));
end
function testDryClosedFormAndSubdivision(tc)
    wing=wing_model('dry'); e=wing.strips; w=1.875104068711961^2*sqrt(e.EI/e.mu)/wing.L^2;
    for pieces=[1 2 4]
        [r,i]=cmodes(wing_cdyn(0,wing,pieces),[-2 2 40 60],'AssumeAnalytic',true);
        tc.verifyTrue(i.complete); tc.verifyEqual(r,1i*w,'AbsTol',1e-8);
        tc.verifyLessThan(max(i.residuals),1e-10);
    end
end
function testTheodorsenAndDeterminant(tc)
    wing=wing_model('goland'); box=[1 10 60 80];
    [d,j]=croots(@(s) wing_delta(s,150,wing),box,'AssumeAnalytic',true);
    for pieces=[1 2]
        [r,i]=cmodes(wing_cdyn(150,wing,pieces),box,'AssumeAnalytic',true);
        tc.verifyTrue(i.complete&&j.complete); tc.verifyEqual(r,d,'AbsTol',1e-8);
    end
end
function testFlutterModeShape(tc)
    % The null vector of K at the flutter root, propagated exactly along the
    % span, satisfies the clamped and free-tip conditions and does not
    % depend on the number of exact pieces used to build K.
    wing=wing_model('goland'); Uf=136.983977449; ratios=zeros(1,2);
    for pieces=[1 2]
        [l,i]=cmodes(wing_cdyn(Uf,wing,pieces),[-2 2 60 80],'AssumeAnalytic',true);
        tc.verifyTrue(i.complete); tc.verifyLessThan(abs(real(l)),1e-6);
        [~,~,z]=wing_mode_shape(l,Uf,wing,i.rightVectors{1},[0 wing.L]);
        tc.verifyLessThan(norm(z(1:3,1))/norm(z(:,1)),1e-12);          % clamped root
        tc.verifyLessThan(norm(z(4:6,end))/norm(z(4:6,1)),1e-10);      % free tip
        ratios(pieces)=wing.strips.b*z(3,end)/z(1,end);
    end
    tc.verifyLessThan(abs(ratios(1)-ratios(2))/abs(ratios(1)),1e-8);
    tc.verifyLessThan(abs(abs(ratios(1))-1.034),5e-3);
    tc.verifyLessThan(abs(angle(ratios(1))*180/pi+62.2),0.5);
end
