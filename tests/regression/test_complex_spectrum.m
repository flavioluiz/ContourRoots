function tests = test_complex_spectrum
%TEST_COMPLEX_SPECTRUM Regression tests of the general solver: analytic
%   spectra, cancellations and failure diagnostics.
    tests = functiontests(localfunctions);
end

function setupOnce(tc)
    repo = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    import matlab.unittest.fixtures.PathFixture
    tc.applyFixture(PathFixture(fullfile(repo,'matlab')));
    tc.applyFixture(PathFixture(fullfile(repo,'matlab','delay')));
    tc.applyFixture(PathFixture(fullfile(repo,'examples','models')));
end

function testMultipleRootAndScaleInvariance(~)
    base={'AssumeAnalytic',true};
    [r,i]=characteristic_roots(@(s) (s-.3-.2i).^3,[-1 1 -1 1],base{:});
    assert(i.complete && numel(r)==1 && i.multiplicity==3 && abs(r-.3-.2i)<1e-6);

    % Same finite roots under a huge change of function scale.
    known=[-.3+.2i; .7-.1i; -.6-.8i];
    for scale=[1e-100 1 1e100]
        f=@(s) scale*(s-known(1)).*(s-known(2)).*(s-known(3));
        [r,i]=characteristic_roots(f,[-1 1 -1 1],base{:});
        check(r,i,known,1e-6);
    end
end

function testBoundaryRootsAndBudgetsAreNotComplete(~)
    base={'AssumeAnalytic',true};
    % Outer-boundary zeros and deliberately exhausted budgets must not pass.
    [~,i]=characteristic_roots(@(s) s,[-1 0 -1 1],base{:});
    assert(~i.complete);
    [~,i]=characteristic_roots(@(s) (s-.3).*(s+.3),[-1 1 -1 1],base{:},'MaxCells',1);
    assert(~i.complete);
    [r,i]=characteristic_roots(@(s) exp(s),[-2 2 -3 3],base{:});
    assert(i.complete && isempty(r));
end

function testHeatModelsAndCancellation(~)
    base={'AssumeAnalytic',true};
    % Heat models: compare with closed-form spectra and sensor cancellation.
    kinds={'heat_neumann','heat_dirichlet','heat_mixed'};
    expected={-(2*pi*(0:1)).^2, -(pi*[1 3]).^2, -(pi*((0:2)+.5)).^2};
    for j=1:3
        G=distributed_model(kinds{j},'Sensor',.5);
        [r,i]=transfer_poles(G,[-100 2 -3 3],base{:});
        check(r,i,expected{j},2e-6);
    end
    % At x=L, Dirichlet input/output is identically one: all modes cancel.
    G=distributed_model('heat_dirichlet','Sensor',1);
    [r,i]=transfer_poles(G,[-100 2 -3 3],base{:});
    assert(i.complete && isempty(r) && numel(i.cancelledLocations)==3);
end

function testPartialCancellationDuctAndWave(~)
    base={'AssumeAnalytic',true};
    % Partial cancellation leaves the correct pole order.
    pair=struct('Numerator',@(s) (s+.2).^2,'Denominator',@(s) (s+.2).^3);
    [r,i]=transfer_poles(pair,[-1 1 -1 1],base{:});
    check(r,i,-.2,1e-6); assert(i.multiplicity==1 && i.cancellationOrders==2);

    G=distributed_model('duct_constant','Reflectance',.5);
    [r,i]=transfer_poles(G,[-2 1 -11 11],base{:});
    check(r,i,log(.5)/2+1i*pi*(-3:3),1e-6);
    G=distributed_model('wave','Damping',0);
    [r,i]=transfer_poles(G,[-1 1 -20 20],base{:});
    k=[-6:-1 1:6]; k=k(mod(k,4)~=0);
    check(r,i,1i*pi*k,2e-6);
    assert(all(abs(i.cancelledLocations)>0)); % no spurious pole at zero
end

function testOpaqueHandleAndSingularities(~)
    base={'AssumeAnalytic',true};
    G=distributed_model('wave','Damping',0);
    % Opaque meromorphic handle with Z-P=0: never claim completeness.
    [r,i]=transfer_poles(@(s) (s-.4)./(s+.2),[-1 1 -1 1],'GridSize',[9 9]);
    assert(~i.complete && strcmp(i.status,'exploratory') && min(abs(r+.2))<1e-6);
    failed=false;
    try
        transfer_poles(G,[-2 2 -3 3],base{:},'Singularities',0);
    catch err
        failed=strcmp(err.identifier,'complex_spectrum:SingularityInRegion');
    end
    assert(failed);
end

function testSymbolicInput(tc)
    tc.assumeTrue(license('test','Symbolic_Toolbox')==1, 'Symbolic_Toolbox not available.');
    syms s;
    Delta=1+s+s^2+(2*s+3)*exp(-s);
    [r,i]=characteristic_roots(Delta,[-8 2 -20 20]);
    [old,oldInfo]=delay_roots([2 3],[1 1 1],1,[-8 2],[-20 20],30,45);
    assert(oldInfo.countMatches); check(r,i,old,2e-6);
    [r,i]=transfer_poles(cosh(s)/sinh(s),[-1 1 -10 10]);
    check(r,i,1i*pi*(-3:3),1e-6);
    [r,i]=transfer_poles(sinh(s/2)/sinh(s),[-1 1 -10 10]);
    check(r,i,1i*pi*[-3 -1 1 3],1e-6);
    failed=false;
    try, characteristic_roots(log(s),[-2 2 -2 2]);
    catch err, failed=strcmp(err.identifier,'complex_spectrum:AnalyticContract'); end
    assert(failed);
end

function testControlSystemInput(tc)
    tc.assumeTrue(license('test','Control_Toolbox')==1, 'Control_Toolbox not available.');
    s=tf('s');
    [r,i]=transfer_poles(exp(-2*s)/(s+1),[-3 1 -3 3]);
    check(r,i,-1,1e-6);
    [r,i]=characteristic_roots((s+2)/(s+1),[-3 1 -3 3]);
    check(r,i,-2,1e-6);
end

function check(r,info,expected,tol)
    assert(info.complete,'Contour/cancellation checks unresolved.');
    assert(numel(r)==numel(expected),'Wrong number of distinct roots.');
    for z=expected(:).', assert(min(abs(r-z))<tol*(1+abs(z))); end
end
