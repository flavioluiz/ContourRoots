function [Delta,G,meta] = coupled_beam_model(varargin)
%COUPLED_BEAM_MODEL Exact Euler-Bernoulli beam coupled to a finite oscillator.
%   [DELTA,G,META]=COUPLED_BEAM_MODEL('Mass',.1,'Stiffness',1.2362,...)
%   DELTA(s) is an entire characteristic function, without modal truncation.
%   Pass it to characteristic_roots(...,'AssumeAnalytic',true).
%
%   Options (consistent SI units, or all dimensionless):
%     Length             beam length L (1)
%     FlexuralRigidity   EI (1)
%     MassPerLength      rho*A (1)
%     Mass               attachment mass m (.1)
%     Stiffness          connection stiffness k (1.2362363368)
%     Damping            connection damping c (.05)
%     Topology           'absorber' (default) or 'grounded-tip'
%
%   absorber: independent mass z(t) connected to tip y=w(L,t) by k,c
%   in parallel, with no connection to ground. The beam is clamped at x=0.
%   grounded-tip: mass rigidly follows the tip; k,c connect it to ground.
%   The beam itself is elastic and undamped; no rotary tip inertia.
%
%   G.TipForceToTip, G.MassForceToTip, G.MassForceToMass are analytic N/D
%   structs (last two are for absorber only). META.Beam is the bare tip
%   receptance. Characteristic roots include hidden modes in decoupled limits;
%   transfer_poles additionally removes input/output cancellations.
    ip=inputParser;
    addParameter(ip,'Length',1,@positive);
    addParameter(ip,'FlexuralRigidity',1,@positive);
    addParameter(ip,'MassPerLength',1,@positive);
    addParameter(ip,'Mass',.1,@nonnegative);
    addParameter(ip,'Stiffness',1.2362363368,@nonnegative);
    addParameter(ip,'Damping',.05,@nonnegative);
    addParameter(ip,'Topology','absorber',@(v) any(strcmpi(v,{'absorber','grounded-tip'})));
    parse(ip,varargin{:}); p=ip.Results;
    if strcmpi(p.Topology,'absorber') && p.Mass==0
        error('coupled_beam_model:Mass','The independent oscillator requires Mass>0.');
    end
    tau=sqrt(p.MassPerLength*p.Length^4/p.FlexuralRigidity);
    compliance=p.Length^3/p.FlexuralRigidity;
    mu=p.Mass/(p.MassPerLength*p.Length);
    kappa=p.Stiffness*compliance;
    gamma=p.Damping*compliance/tau;
    A=@(s) beam_A(-(tau*s).^2);
    B=@(s) beam_B(-(tau*s).^2);
    M=@(s) mu*(tau*s).^2;
    Z=@(s) kappa+gamma*tau*s;
    if strcmpi(p.Topology,'absorber')
        Delta=@(s) B(s).*(M(s)+Z(s))+A(s).*M(s).*Z(s);
        G.TipForceToTip=struct('Numerator',@(s) compliance*A(s).*(M(s)+Z(s)), ...
            'Denominator',Delta);
        G.MassForceToTip=struct('Numerator',@(s) compliance*A(s).*Z(s), ...
            'Denominator',Delta);
        G.MassForceToMass=struct('Numerator',@(s) compliance*(B(s)+A(s).*Z(s)), ...
            'Denominator',Delta);
        scaledResidual=@(s) abs(Delta(s))./(abs(B(s).*(M(s)+Z(s)))+ ...
            abs(A(s).*M(s).*Z(s))+realmin);
    else
        Delta=@(s) B(s)+A(s).*(M(s)+Z(s));
        G.TipForceToTip=struct('Numerator',@(s) compliance*A(s),'Denominator',Delta);
        scaledResidual=@(s) abs(Delta(s))./(abs(B(s))+abs(A(s).*(M(s)+Z(s)))+realmin);
    end
    meta=struct('parameters',p,'timeScale',tau,'massRatio',mu, ...
        'dimensionlessStiffness',kappa,'dimensionlessDamping',gamma, ...
        'Beam',struct('Numerator',@(s) compliance*A(s),'Denominator',B), ...
        'normalizedResidual',scaledResidual);
end

function v=beam_A(q)
% [cosh(b)sin(b)-sinh(b)cos(b)]/b^3, b^4=q, entire in q.
    small=abs(q)<1e-3; v=zeros(size(q));
    z=q(small); v(small)=2/3-z/315+z.^2/623700-z.^3/5108103000;
    b=q(~small).^.25;
    v(~small)=(cosh(b).*sin(b)-sinh(b).*cos(b))./b.^3;
end
function v=beam_B(q)
    b=q.^.25; v=1+cosh(b).*cos(b);
end
function yes=positive(x)
    yes=isnumeric(x) && isscalar(x) && isreal(x) && isfinite(x) && x>0;
end
function yes=nonnegative(x)
    yes=isnumeric(x) && isscalar(x) && isreal(x) && isfinite(x) && x>=0;
end
