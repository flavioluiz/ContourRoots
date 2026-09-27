function wing = wing_model(kind,nStrips)
%WING_MODEL Continuous cantilever wing: bending, torsion and strip aerodynamics.
%   WING = WING_MODEL('goland') returns the classical Goland wing (SI units):
%   a straight, uniform, clamped-free wing that bends (Euler-Bernoulli) and
%   twists (Saint-Venant), with inertial coupling between the two (the
%   centre of mass lies aft of the elastic axis). Each spanwise strip carries
%   the exact two-dimensional Theodorsen air loads (THEODORSEN_LOADS).
%
%   WING = WING_MODEL(KIND,N) splits the span into N strips of constant
%   properties. KIND is
%     'goland'   uniform wing: every N gives the SAME model (the strips are
%                a bookkeeping device, not a discretization);
%     'dry'      the same structure in vacuum (rho = 0) and without
%                inertial coupling (S = 0): the classical cantilever
%                frequencies are known in closed form;
%     'tapered'  a synthetic wing whose chord, stiffness and inertia vary
%                along the span; the N strips approximate that variation
%                by piecewise-constant properties (converges as N grows).
%
%   Fields: L (span, m), strips (struct array, one per strip) and scale
%   (fixed scaling of the state [w; w'; alpha; V; M; T], see WING_MATRIX).
%   Strip fields: length, EI, GJ (N m^2), mu (kg/m), S (kg), Ialpha (kg m),
%   b (semichord, m), a (elastic axis aft of midchord, in semichords),
%   rho (air density, kg/m^3).
%
%   Parameters: M. Goland, The flutter of a uniform cantilever wing,
%   J. Appl. Mech. 12 (1945) A197-A208, as tabulated by A. A. Cal (PhD thesis,
%   City University London, 1992, Table 2.4), with h positive downward.
%
%   See also WING_MATRIX, WING_DELTA, WING_TRANSFER, THEODORSEN_LOADS.
    if nargin<1, kind='goland'; end
    if nargin<2, nStrips=1; end
    validateattributes(nStrips,{'numeric'},{'scalar','integer','positive'});
    L=6.096;                                   % span, 20 ft
    strip=struct('length',L/nStrips, ...
        'EI',9.773e6,'GJ',9.876e5, ...         % bending and torsion stiffness
        'mu',35.717,'S',35.717*0.183, ...      % mass/span; static moment (c.g. 0.183 m aft of e.a.)
        'Ialpha',8.642, ...                    % pitch inertia/span about the e.a.
        'b',0.9144,'a',-0.34,'rho',1.225);     % semichord; e.a. at 33% chord; sea-level air
    switch lower(kind)
        case 'goland'
        case 'dry', strip.rho=0; strip.S=0;
        case 'tapered'
        otherwise, error('wing:Model','Use ''goland'', ''dry'' or ''tapered''.');
    end
    wing=struct('name',lower(kind),'L',L,'strips',repmat(strip,1,nStrips));
    if strcmpi(kind,'tapered')
        for j=1:nStrips
            c=1-0.3*(j-0.5)/nStrips;           % chord ratio at the strip midpoint
            wing.strips(j).b=strip.b*c;
            wing.strips(j).EI=strip.EI*c^3;    wing.strips(j).GJ=strip.GJ*c^3;
            wing.strips(j).mu=strip.mu*c^2;    wing.strips(j).S=strip.S*c^3;
            wing.strips(j).Ialpha=strip.Ialpha*c^4;
        end
    end
    % Fixed scales of [w; w'; alpha; V; M; T]; they never depend on s.
    wing.scale=[L;1;1;strip.EI/L^2;strip.EI/L;strip.GJ/L];
    check(wing);
end

function check(wing)
    for e=wing.strips
        for name={'length','EI','GJ','mu','Ialpha','b'}
            validateattributes(e.(name{1}),{'numeric'},{'scalar','positive','finite'});
        end
        if e.mu*e.Ialpha-e.S^2<=0
            error('wing:Mass','The section inertia matrix must be positive definite.');
        end
    end
    if abs(sum([wing.strips.length])-wing.L)>1e-12*wing.L
        error('wing:Length','Strip lengths must add up to the span.');
    end
end
