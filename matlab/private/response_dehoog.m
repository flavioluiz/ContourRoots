function [y,i]=response_dehoog(f,t,o,bound)
% Double-precision de Hoog inversion with degree, period and shift checks.
% Q-D recurrence adapted from mpmath (BSD-3-Clause); see THIRD_PARTY_NOTICES.
    i=response_info(t,o); y=NaN(size(t)); history=zeros(0,5);
    % Each level DOUBLES the degree, hence the sampled bandwidth 2M*2*pi/P.
    % A smaller ratio (e.g. 1.4) lets consecutive degrees share a false
    % plateau that omits lightly damped modes just above the band.
    degrees=[10 20 40 80 160 320];
    for k=1:numel(t)
        P=4*t(k); sigma=response_line(o,bound,P);
        i.abscissa=maxfinite(i.abscissa,sigma); amp=exp(max(0,sigma)*t(k));
        i.amplification=max(i.amplification,amp);
        if (sigma+.5/P)*t(k)>log(1e12), continue; end
        previous=NaN; success=0;
        for level=1:min(o.MaxRefinements,numel(degrees))
            M=degrees(level);
            if 4*M+1>o.MaxPoints || (4*M+1)^2*64>o.MaxMemoryMB*2^20, break; end
            [a,ea,ca]=one(f,t(k),P,sigma,M);
            % A modest period perturbation avoids oversampling the Q-D
            % table into severe near-dependence in double precision.
            L=ceil(1.25*M);
            [b,eb,cb]=one(f,t(k),1.25*P,sigma,L);
            [c,ec,cc]=one(f,t(k),1.25*P,sigma+.5/P,L);
            ep=abs(b-a); es=abs(c-b); ed=abs(b-previous);
            err=max([ep es ed])+eps*exp(max(0,sigma+.5/P)*t(k))*max([ca cb cc])*100;
            % MATLAB max omits NaNs: the first degree has NO previous value.
            % It must not count as a successful degree-convergence check.
            if any(~isfinite([a b c previous err])), err=Inf; end
            y(k)=b; i.errorEstimate(k)=err; i.periodError(k)=ep; i.bandwidthError(k)=ed; i.shiftError(k)=es;
            i.evaluations=i.evaluations+4*M+8*L+6; i.points=max(i.points,2*L+1);
            i.symmetryError=max([i.symmetryError ea eb ec]);
            history(end+1,:)=[t(k) M ep ed es]; %#ok<AGROW>
            if err<=o.AbsTol+o.RelTol*abs(b), success=success+1; else, success=0; end
            % Low degrees can share a false plateau before an oscillatory
            % mode enters the sampled band. Require a degree-40 guard probe
            % and two successive doublings in agreement. This is still
            % finite resolution, not a global spectral bound.
            if success>=2 && M>=40, i.resolvedMask(k)=true; break; end
            previous=b;
        end
    end
    i.history=history; i.converged=all(i.resolvedMask);
    i.stopReason='Degree/period/shift refinement incomplete or conditioning limit.';
    if i.converged, i.status='converged'; i.stopReason='Two successive independent refinement checks passed at every positive time.'; end
end
function [v,symerr,condition]=one(f,t,P,sigma,M)
    s=sigma+1i*(0:2*M).'*2*pi/P; fp=response_values(f,s);
    neg=response_values(f,conj(s)); symerr=max(abs(fp-conj(neg)))/max(realmin,max(abs(fp)));
    condition=2*sum(abs(fp))/P;
    if symerr>1e-9, error('ContourRoots:ResponseSymmetry','Transfer is not conjugate symmetric.'); end
    v=response_dehoog_value(fp,t,P,sigma,M);
end

function c=maxfinite(a,b)
    if isnan(a), c=b; else, c=max(a,b); end
end
