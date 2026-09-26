function [y,i]=response_dehoog(f,t,o,bound)
% Double-precision de Hoog inversion with degree, period and shift checks.
% Q-D recurrence adapted from mpmath (BSD-3-Clause); see THIRD_PARTY_NOTICES.
    i=response_info(t,o); y=NaN(size(t)); history=zeros(0,5);
    degrees=[12 20 28 40 56 80 112 160];
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
            % mode enters the sampled band. Require a degree-40 guard probe.
            % This is still finite resolution, not a global spectral bound.
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
    if all(fp==0), v=0; return; end
    e=complex(zeros(2*M+1,M+1)); q=complex(zeros(2*M,M));
    q(1,1)=fp(2)/(fp(1)/2); q(2:2*M,1)=fp(3:end)./fp(2:end-1);
    for r=1:M
        n=2*(M-r)+1;
        e(1:n,r+1)=q(2:n+1,r)-q(1:n,r)+e(2:n+1,r);
        if r<M, q(1:n-1,r+1)=q(2:n,r).*e(2:n,r+1)./e(1:n-1,r+1); end
    end
    d=complex(zeros(2*M+1,1)); d(1)=fp(1)/2;
    for r=1:M, d(2*r)=-q(1,r); d(2*r+1)=-e(1,r+1); end
    z=exp(2i*pi*t/P); ap=0; ac=d(1); bp=1; bc=1;
    for j=2:2*M
        an=ac+d(j)*ap*z; bn=bc+d(j)*bp*z;
        ap=ac; ac=an; bp=bc; bc=bn;
        scale=max([abs(ac),abs(ap),abs(bc),abs(bp)]);
        if scale>1e100, ac=ac/scale; ap=ap/scale; bc=bc/scale; bp=bp/scale; end
    end
    brem=(1+(d(2*M)-d(2*M+1))*z)/2;
    x=d(2*M+1)*z/brem;
    remainder=brem*x/(sqrt(1+x)+1); % stable sqrt(1+x)-1
    v=exp(sigma*t)*2/P*real((ac+remainder*ap)/(bc+remainder*bp));
end
function c=maxfinite(a,b)
    if isnan(a), c=b; else, c=max(a,b); end
end
