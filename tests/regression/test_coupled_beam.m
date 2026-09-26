function tests = test_coupled_beam
%TEST_COUPLED_BEAM Static compliance, scaling, limits, spectra and an
%   independent finite-element check of the coupled beam example.
    tests = functiontests(localfunctions);
end

function setupOnce(tc)
    repo = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    import matlab.unittest.fixtures.PathFixture
    tc.applyFixture(PathFixture(fullfile(repo,'matlab')));
    tc.applyFixture(PathFixture(fullfile(repo,'examples','models')));
end

function testSpectrumAndFiniteElementConvergence(~)
    base={'AssumeAnalytic',true}; region=[-30 1 -130 130];
    [D,G,meta]=coupled_beam_model();
    assert(abs(G.TipForceToTip.Numerator(0)/D(0)-1/3)<1e-12);
    assert(abs(G.MassForceToMass.Numerator(0)/D(0)-(1/3+1/meta.parameters.Stiffness))<1e-12);
    [r,i]=characteristic_roots(D,region,base{:}); assert(i.complete && numel(r)==10);
    assert(max(real(r))<0 && max(meta.normalizedResidual(r))<1e-8);
    [p,j]=transfer_poles(G.TipForceToTip,region,base{:}); assert(j.complete && numel(p)==10);

    % Independent FE convergence, including all poles in this window.
    errors=zeros(1,3); meshes=[8 16 32];
    for a=1:3
        pf=coupled_beam_fem(meta.parameters,meshes(a));
        errors(a)=max(arrayfun(@(z) min(abs(pf-z)),r));
    end
    assert(all(diff(errors)<0) && errors(end)<.002);
end

function testDecoupledOscillatorHiddenModes(~)
    base={'AssumeAnalytic',true};
    % Decoupled oscillator: two internal zero modes disappear from tip TF.
    [D,G]=coupled_beam_model('Stiffness',0,'Damping',0);
    [r,i]=characteristic_roots(D,[-1 1 -30 30],base{:}); assert(i.complete);
    q=find(abs(r)<1e-5); assert(isscalar(q) && i.multiplicity(q)==2);
    [p,j]=transfer_poles(G.TipForceToTip,[-1 1 -30 30],base{:});
    assert(j.complete && numel(p)==4 && all(abs(p)>1));
end

function testBareBeamLimit(~)
    base={'AssumeAnalytic',true};
    % No boundary attachments reproduces the first bare-beam eigenfrequency.
    [D,G]=coupled_beam_model('Topology','grounded-tip','Mass',0,'Stiffness',0,'Damping',0);
    [r,i]=characteristic_roots(D,[-1 1 -5 5],base{:});
    beta=fzero(@(b) cos(b)+1/cosh(b),[1 2]);
    assert(i.complete && numel(r)==2 && max(abs(abs(imag(r))-beta^2))<1e-7);
    assert(abs(G.TipForceToTip.Numerator(0)/D(0)-1/3)<1e-12);
end

function testDimensionalScaling(~)
    base={'AssumeAnalytic',true};
    % Dimensional rescaling leaves nondimensional characteristic unchanged.
    [Dn,~,mn]=coupled_beam_model('Mass',.2,'Stiffness',2,'Damping',.1);
    L=2; EI=9; rhoA=3; tau=sqrt(rhoA*L^4/EI);
    [Dp,~,mp]=coupled_beam_model('Length',L,'FlexuralRigidity',EI, ...
        'MassPerLength',rhoA,'Mass',.2*rhoA*L, ...
        'Stiffness',2*EI/L^3,'Damping',.1*EI*tau/L^3);
    z=[0 -1+2i -.1+25i];
    assert(max(abs(Dn(z)-Dp(z/tau))./(1+abs(Dn(z))))<1e-12);
    assert(abs(mn.massRatio-mp.massRatio)<1e-12);
end
