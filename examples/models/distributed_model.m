function [G,meta] = distributed_model(kind,varargin)
%DISTRIBUTED_MODEL Nonrational benchmarks from Curtain & Morris (2009).
%   [G,META]=DISTRIBUTED_MODEL(KIND,...) returns analytic N/D handles for
%   COMPLEX_SPECTRUM. Evaluators regularize removable singularities at zero.
%   Models are scalar input/output maps, not full PDE spectra.
%
%   KIND: heat_neumann, heat_dirichlet, heat_mixed (normalized L=alpha=K0=1)
%         wave (normalized L=c=1, b=1 on [0,1/2], damping by output feedback)
%         duct_constant (normalized L=c=rho=1)
%         duct_radiation (dimensional parameters from Fig. 5)
%         beam (normalized L=EI=I=1, Kelvin-Voigt damping)
%   Options: Sensor (.5), Damping (.02), Reflectance (.5).
%   For beam, exclude META.singularities = -1/Damping from the closed region.
    p=inputParser;
    addParameter(p,'Sensor',.5,@(x) isscalar(x) && isreal(x) && x>0 && x<=1);
    addParameter(p,'Damping',.02,@(x) isscalar(x) && isreal(x) && isfinite(x) && x>=0);
    addParameter(p,'Reflectance',.5,@(x) isscalar(x) && isreal(x) && x>0 && x<1);
    parse(p,varargin{:}); o=p.Results;
    x=o.Sensor; e=o.Damping;
    meta=struct('kind',kind,'singularities',[],'reference', ...
        'Curtain & Morris, Automatica 45 (2009), 1101-1116, doi:10.1016/j.automatica.2009.01.008', ...
        'parameters',o);
    switch lower(kind)
        case 'heat_neumann'
            N=@(s) cosh(x*sqrt(s)); D=@(s) s.*sinhc_sqrt(s);
        case 'heat_dirichlet'
            N=@(s) x*sinhc_sqrt(x^2*s); D=@sinhc_sqrt;
        case 'heat_mixed'
            N=@(s) x*sinhc_sqrt(x^2*s); D=@(s) cosh(sqrt(s));
        case 'wave'
            N=@(s) s.*wave_n(s);
            D=@(s) sinhc(s)+e*N(s);
        case 'duct_constant'
            a=o.Reflectance;
            N=@(s) exp(-x*s).*(1+a*exp(2*(x-1)*s));
            D=@(s) 1-a*exp(-2*s);
        case 'duct_radiation'
            % Article Fig. 5: length [m], sound speed [m/s], density [kg/m^3].
            L=3.54; c=341; rho=1.20; a=.101; x0=x*L;
            R2=rho*c/(pi*a^2); R1=.504*R2;
            C=5.44*a^3/(rho*c^2); M=.1952*rho/a;
            zN=pi*a^2*[R1*R2*M*C (R1+R2)*M 0];
            zD=[R1*M*C M+R1*R2*C R1+R2];
            aN=zN-rho*c*zD; aD=zN+rho*c*zD;
            % Clear the reflectance denominator in BOTH factors.
            N=@(s) rho*c*exp(-x0*s/c).*(polyval(aD,s)+polyval(aN,s).*exp(2*(x0-L)*s/c));
            D=@(s) polyval(aD,s)-polyval(aN,s).*exp(-2*L*s/c);
            meta.parameters.Length=L; meta.parameters.SoundSpeed=c;
            meta.parameters.Density=rho; meta.parameters.Radius=a;
        case 'beam'
            % Correct shear-force / tip-velocity model, section 4.2.
            q=@(s) -s.^2./(1+e*s);
            N=@(s) s.*beam_n(q(s));
            D=@(s) (1+e*s).*beam_d(q(s));
            if e>0, meta.singularities=-1/e; end
        otherwise
            error('distributed_model:Kind','Unknown distributed model: %s',kind);
    end
    G=struct('Numerator',N,'Denominator',D);
end

function y=sinhc_sqrt(q)
% Entire function sum(q^k/(2k+1)!), independent of sqrt branch.
    y=ones(size(q)); small=abs(q)<1e-4;
    if any(small(:))
        v=q(small); y(small)=1+v/6+v.^2/120+v.^3/5040+v.^4/362880;
    end
    v=sqrt(q(~small)); y(~small)=sinh(v)./v;
end
function y=sinhc(s)
    y=sinhc_sqrt(s.^2);
end
function y=wave_n(s)
% Remove the fourth-order zero of the wave numerator analytically.
    y=zeros(size(s)); small=abs(s)<.5;
    v=s(small);
    for k=2:12
        a=.5/factorial(2*k-1)+(2^(1-2*k)-1.5)/factorial(2*k);
        y(small)=y(small)+a*v.^(2*k-4);
    end
    v=s(~small);
    y(~small)=(.5*v.*sinh(v)+2*cosh(v/2)-3*cosh(v/2).^2+1)./v.^4;
end
function y=beam_d(q)
% cosh(m)cos(m) is invariant under m -> -m and m -> i*m, m^4=q.
    m=q.^.25; y=1+cosh(m).*cos(m);
end
function y=beam_n(q)
% [cosh(m)sin(m)-sinh(m)cos(m)]/m^3; regular at m=0.
    small=abs(q)<1e-3; y=zeros(size(q));
    v=q(small); y(small)=2/3-v/315+v.^2/623700-v.^3/5108103000;
    m=q(~small).^.25;
    y(~small)=(cosh(m).*sin(m)-sinh(m).*cos(m))./m.^3;
end
