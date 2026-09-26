function [v,i]=response_kernels(m,t,power,o)
% Invert the regular part, with exact transport shift when metadata exists.
    tau=t-m.delay;
    % A sample at the delay time must use the exact right limit at tau=0,
    % not a rounding-error positive tau (e.g. 0.35-0.35 = 5.6e-17).
    tau(abs(tau)<=8*eps*max(1,abs(t)))=0;
    active=tau>=0; v=zeros(size(t)); i=response_info(t,o);
    i.errorEstimate(:)=0; i.periodError(:)=0; i.bandwidthError(:)=0; i.shiftError(:)=0; i.resolvedMask(:)=true;
    if any(active) && ~m.zero
        f=m.regular; initial=0; add=zeros(sum(active),1);
        if power==0
            initial=m.initial;
            if isfinite(initial) && initial~=0
                % Exact subtraction/addition, not a rational approximation.
                h0=initial; beta=1;
                if ~isempty(m.bound)&&isfinite(m.bound), beta=max(1,abs(m.bound)); end
                % The auxiliary pole must remain LEFT of an explicit line.
                if ~isempty(o.Abscissa), beta=max(beta,1-o.Abscissa); end
                f=@(s) response_values(m.regular,s)-h0./(s+beta);
                add=h0*exp(-beta*tau(active)); initial=0;
            end
        else
            f=@(s) response_values(m.regular,s)./s.^power;
            o.AssumeStable=false; % Integrated kernels need a line right of s=0.
        end
        b=m.bound;
        if power>0 && ~isempty(b), b=max(0,b); end
        if power>0 && ~isempty(o.Abscissa) && o.Abscissa<=0
            error('ContourRoots:ResponseDomain','Step/ramp inversion requires a positive Abscissa (the input adds a singularity at zero).');
        end
        [w,j]=response_invert(f,tau(active),o,b,initial); v(active)=w+add;
        i=j; i.time=t;
        for field={'errorEstimate','resolvedMask','periodError','bandwidthError','shiftError'}
            key=field{1}; if strcmp(key,'resolvedMask'), q=true(size(t)); else, q=zeros(size(t)); end
            q(active)=j.(key); i.(key)=q;
        end
    else
        i.stopReason='Exact zero regular response on this interval.';
    end
    i.converged=all(i.resolvedMask); if i.converged, i.status='converged'; end
end
