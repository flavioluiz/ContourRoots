function predicted = pade_critical_delays(N,D,Tmax,orders,varargin)
%PADE_CRITICAL_DELAYS Imaginary-axis crossings predicted by diagonal Pade.
%   PREDICTED = PADE_CRITICAL_DELAYS(N,D,TMAX,ORDERS) finds every T in
%   [TMIN,TMAX] for which
%       D(i*w) + R_n(i*w*T) N(i*w) = 0,
%   where R_n = P_n/Q_n is the diagonal [n/n] Pade approximation of
%   exp(-s*T), for each order n in ORDERS.
%
%   A diagonal Pade approximation is all-pass on the imaginary axis,
%   |R_n(i*y)| = 1, so the crossing frequencies are those of the exact
%   equation (from CRITICAL_DELAYS). Its phase phi_n(y) = -2 arg Q_n(i*y)
%   decreases monotonically from 0 to -n*pi. For each frequency w with exact
%   phase offset theta (see CRITICAL_DELAYS), the branch k exists if and
%   only if theta + 2*pi*k < n*pi, and then the unique solution of
%   phi_n(w*T) = -(theta + 2*pi*k) is computed by bracketing and FZERO. No
%   sampling grid is used, so enlarging [TMIN,TMAX] never removes a
%   crossing found in a smaller interval.
%
%   Options:
%     'Tmin'          lower delay bound (default 0)
%     'Tolerance'     residual tolerance, relative to |D(i*w)| (default 1e-9)
%     'GridSize'      ignored; kept for compatibility with version 0.1
%
%   Output table columns: Order, Frequency, Delay, PhaseDefect (Pade phase
%   minus exact phase -w*T at the predicted delay), Residual.
%
%   See also CRITICAL_DELAYS, PADE_DELAY, PADE_CHARACTERISTIC.

    validateattributes(Tmax,{'numeric'},{'scalar','real','finite','nonnegative'});
    validateattributes(orders,{'numeric'},{'vector','integer','positive'});
    p=inputParser;
    addParameter(p,'Tmin',0,@(x) isscalar(x) && x>=0);
    addParameter(p,'GridSize',12000,@(x) isscalar(x));
    addParameter(p,'Tolerance',1e-9,@(x) isscalar(x) && x>0);
    parse(p,varargin{:}); opt=p.Results;
    assert(opt.Tmin<=Tmax,'Tmin must not exceed Tmax.');

    % Frequencies are independent of the delay interval.
    w = warning('off','critical_delays:PersistentImaginaryRoots');
    restore = onCleanup(@() warning(w));
    [~,info]=critical_delays(N,D,0);
    frequencies=info.candidateFrequencies(:);
    rows=zeros(0,5);

    for order=unique(round(orders(:).'))
        [~,qn]=pade_delay(1,order);                 % Q_n(z), z = s*T
        qRoots=roots(qn);                           % all in Re z < 0
        phi=@(y) pade_phase(qRoots,y);
        for jf=1:numel(frequencies)
            omega=frequencies(jf);
            Dw=polyval(D,1i*omega); Nw=polyval(N,1i*omega);
            theta=mod(-angle(-Dw/Nw),2*pi);
            k=0;
            while theta+2*pi*k < order*pi
                target=-(theta+2*pi*k);
                y=solve_phase(phi,target);
                Tc=y/omega;
                if Tc>=opt.Tmin && Tc<=Tmax
                    z=1i*omega*Tc;
                    [pn,qn1]=pade_delay(1,order);
                    R=polyval(pn,z)/polyval(qn1,z);
                    residual=abs(Dw+R*Nw);
                    if residual <= opt.Tolerance*max(abs(Dw),realmin)
                        rows(end+1,:)=[order omega Tc phi(y)+y residual]; %#ok<AGROW>
                    end
                end
                if Tc>Tmax, break; end
                k=k+1;
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

function ph=pade_phase(qRoots,y)
% Continuous phase of P_n(iy)/Q_n(iy) = conj(Q_n(iy))/Q_n(iy): -2 arg Q_n(iy),
% written as a sum of continuous terms (every root has Re < 0).
    ph=zeros(size(y));
    for r=qRoots(:).'
        ph=ph-2*(atan2(y-imag(r),-real(r))-atan2(-imag(r),-real(r)));
    end
end

function y=solve_phase(phi,target)
% phi decreases monotonically from 0 to -n*pi; find phi(y) = target.
    hi=1;
    while phi(hi)>target
        hi=2*hi;
    end
    y=fzero(@(t) phi(t)-target,[0 hi]);
end
