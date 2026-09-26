function predicted = pade_critical_delays(N,D,Tmax,orders,varargin)
%PADE_CRITICAL_DELAYS Imaginary-axis crossings predicted by diagonal Pade.
%   PREDICTED = PADE_CRITICAL_DELAYS(N,D,TMAX,ORDERS) finds every T in
%   [0,TMAX] for which
%       D(i*w) + R_n(i*w,T) N(i*w) = 0,
%   where R_n is the diagonal [n/n] Pade approximation of exp(-s*T).
%
%   A diagonal Pade approximation is all-pass on the imaginary axis:
%   |R_n(i*w,T)|=1. Therefore it preserves the exact magnitude equation
%   |D(i*w)|=|N(i*w)|. It can still shift or miss critical delays because
%   its phase approximates -w*T only over a finite range.
%
%   Options:
%     'Tmin'          lower delay bound (default 0)
%     'GridSize'      phase-search samples (default 12000)
%     'Tolerance'     residual and clustering tolerance (default 1e-9)

    validateattributes(Tmax,{'numeric'},{'scalar','real','finite','nonnegative'});
    validateattributes(orders,{'numeric'},{'vector','integer','positive'});
    p=inputParser;
    addParameter(p,'Tmin',0,@(x) isscalar(x) && x>=0);
    addParameter(p,'GridSize',12000,@(x) isscalar(x) && x>=500);
    addParameter(p,'Tolerance',1e-9,@(x) isscalar(x) && x>0);
    parse(p,varargin{:}); opt=p.Results;
    assert(opt.Tmin<=Tmax,'Tmin must not exceed Tmax.');

    % The exact routine exposes every candidate frequency independently of
    % whether its exact phase branch falls inside the requested delay range.
    [~,info]=critical_delays(N,D,Tmax,'Tmin',opt.Tmin);
    frequencies=info.candidateFrequencies(:);
    rows=zeros(0,5);

    for order=unique(round(orders(:).'))
        [pn,qn]=pade_delay(1,order); % polynomials in z=s*T
        for jf=1:numel(frequencies)
            omega=frequencies(jf);
            target=-polyval(D,1i*omega)/polyval(N,1i*omega);
            target=target/abs(target);
            grid=linspace(opt.Tmin,Tmax,round(opt.GridSize));
            z=1i*omega*grid;
            ratio=polyval(pn,z)./polyval(qn,z);
            h=imag(ratio.*conj(target));
            brackets=find(h(1:end-1).*h(2:end)<=0);
            candidates=zeros(0,1);
            for ib=brackets(:).'
                a=grid(ib); b=grid(ib+1);
                fun=@(T) imag(pade_ratio(pn,qn,omega,T)*conj(target));
                if abs(h(ib))<10*eps
                    Tc=a;
                else
                    try
                        Tc=fzero(fun,[a b]);
                    catch
                        continue
                    end
                end
                R=pade_ratio(pn,qn,omega,Tc);
                residual=abs(polyval(D,1i*omega)+R*polyval(N,1i*omega));
                if real(R*conj(target))>0 && residual<opt.Tolerance*(1+abs(polyval(D,1i*omega)))
                    if isempty(candidates) || all(abs(Tc-candidates)>100*opt.Tolerance*(1+Tc))
                        candidates(end+1,1)=Tc; %#ok<AGROW>
                    end
                end
            end
            for k=1:numel(candidates)
                Tc=candidates(k);
                R=pade_ratio(pn,qn,omega,Tc);
                residual=abs(polyval(D,1i*omega)+R*polyval(N,1i*omega));
                phaseExact=-omega*Tc;
                phasePade=continuous_pade_phase(pn,qn,omega,Tc,600);
                rows(end+1,:)=[order omega Tc phasePade-phaseExact residual]; %#ok<AGROW>
            end
        end
    end

    if isempty(rows)
        predicted=table(zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
            'VariableNames',{'Order','Frequency','Delay','PhaseDefect','Residual'});
    else
        [~,idx]=sortrows(rows(:,[1 3 2]),[1 2 3]); rows=rows(idx,:);
        predicted=array2table(rows,'VariableNames', ...
            {'Order','Frequency','Delay','PhaseDefect','Residual'});
        predicted.Order=round(predicted.Order);
    end
end

function R=pade_ratio(pn,qn,omega,T)
    z=1i*omega*T;
    R=polyval(pn,z)/polyval(qn,z);
end

function phase=continuous_pade_phase(pn,qn,omega,T,npts)
    tau=linspace(0,T,max(20,npts));
    z=1i*omega*tau;
    values=polyval(pn,z)./polyval(qn,z);
    p=unwrap(angle(values));
    phase=p(end);
end
