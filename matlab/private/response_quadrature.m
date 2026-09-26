function [y,i]=response_quadrature(f,t,o,bound)
% Direct adaptive Bromwich quadrature, independent of Fourier acceleration.
% Finite-band tail/line checks are estimates, not certified tail bounds.
    i=response_info(t,o); y=NaN(size(t)); history=zeros(0,5);
    for k=1:numel(t)
        tk=t(k);
        if ~isempty(o.Abscissa)
            sigma=o.Abscissa;
        elseif o.AssumeStable
            sigma=0;
        else
            sigma=max(0,bound)+1/tk;
        end
        if isnan(i.abscissa), i.abscissa=sigma; else, i.abscissa=max(i.abscissa,sigma); end
        amp=exp(max(0,sigma+.25/tk)*tk); i.amplification=max(i.amplification,amp);
        if ~isfinite(amp) || amp>1e12, continue; end
        W=32*pi/tk; previous=NaN; success=0;
        for level=1:o.MaxRefinements
            panels=ceil(2*W*tk/pi);
            if panels*32>o.MaxPoints || panels*32*128>o.MaxMemoryMB*2^20, break; end
            budget=min(o.MaxPoints,floor(o.MaxMemoryMB*2^20/128));
            [a,ea,na,sa]=one(W,sigma,1);
            [b,eb,nb,sb]=one(2*W,sigma,.25);
            [c,ec,nc,sc]=one(2*W,sigma+.25/tk,.25);
            tail=abs(b-a); shift=abs(c-b); degree=abs(b-previous);
            err=max([tail shift degree ea+eb ec+eb]);
            if any(~isfinite([a b c previous err])), err=Inf; end
            y(k)=b; i.errorEstimate(k)=err; i.bandwidthError(k)=tail;
            i.periodError(k)=0; i.shiftError(k)=shift;
            i.evaluations=i.evaluations+na+nb+nc; i.points=max([i.points na nb nc]);
            i.bandwidth=max(i.bandwidth,2*W); i.symmetryError=max([i.symmetryError sa sb sc]);
            history(end+1,:)=[tk 2*W tail degree shift]; %#ok<AGROW>
            if err<=o.AbsTol+o.RelTol*abs(b), success=success+1; else, success=0; end
            if success>=2, i.resolvedMask(k)=true; break; end
            previous=b; W=2*W;
        end
    end
    i.history=history; i.converged=all(i.resolvedMask);
    i.stopReason='Bromwich quadrature bandwidth/tolerance/shift checks incomplete or budget exhausted.';
    if i.converged
        i.status='converged'; i.stopReason='Two successive bandwidth/quadrature/shift checks passed.';
    end
    function [v,err,count,symerr]=one(cut,line,tighten)
        count=0; symerr=0; v=NaN; err=Inf;
        scale=exp(line*tk)/pi;
        knots=linspace(0,cut,ceil(cut*tk/pi)+1);
        try
            [q,e]=quadgk(@integrand,0,cut,'Waypoints',knots(2:end-1), ...
                'AbsTol',o.AbsTol*tighten/(16*scale), ...
                'RelTol',max(100*eps,o.RelTol*tighten/16), ...
                'MaxIntervalCount',max(650,floor(budget/32)));
            v=scale*q; err=scale*e;
        catch ex
            if ~strcmp(ex.identifier,'ContourRoots:QuadratureBudget'), rethrow(ex); end
        end
        function z=integrand(w)
            count=count+2*numel(w);
            if count>budget, error('ContourRoots:QuadratureBudget','Quadrature evaluation budget exhausted.'); end
            p=response_values(f,line+1i*w); n=response_values(f,line-1i*w);
            symmetry=max(abs(p(:)-conj(n(:))))/max(realmin,max(abs(p(:))));
            symerr=max(symerr,symmetry);
            if symmetry>1e-9, error('ContourRoots:ResponseSymmetry','Transfer is not conjugate symmetric.'); end
            z=real(p.*exp(1i*w*tk));
        end
    end
end
