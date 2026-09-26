function [y,i]=response_fft(f,t,o,bound)
% Separate period, bandwidth and line-shift checks; finite-band heuristic.
    i=response_info(t,o); y=NaN(size(t)); horizon=max(t);
    if numel(t)>1, dt=min(diff(t)); else, dt=horizon/32; end
    P=4*max(horizon,dt); N=2^nextpow2(max(1024,4*P/dt));
    sigma=response_line(o,bound,P); success=zeros(size(t)); history=zeros(0,5);
    i.abscissa=sigma; i.amplification=exp(max(0,sigma)*horizon);
    if max(0,sigma+.5/P)*horizon>log(1e12)
        i.stopReason='Exponential amplification exceeds 1e12; use shorter time blocks.'; return
    end
    for level=1:o.MaxRefinements
        sigma=response_line(o,bound,P);
        i.abscissa=sigma; i.amplification=exp(max(0,sigma)*horizon);
        if 4*N>o.MaxPoints || 4*N*128>o.MaxMemoryMB*2^20
            i.stopReason='FFT point/memory budget exhausted.'; break
        end
        [base,e1,c1]=sample(f,t,P,N,sigma);
        [per,e2,c2]=sample(f,t,2*P,2*N,sigma);
        [band,e3,c3]=sample(f,t,2*P,4*N,sigma);
        [shift,e4,c4]=sample(f,t,2*P,4*N,sigma+.5/P);
        ep=abs(per-base); eb=abs(band-per); es=abs(shift-band);
        floorError=eps*exp(max(0,sigma+.5/P)*horizon)*max([c1 c2 c3 c4])*20;
        err=max([ep eb es],[],2)+floorError;
        err(any(~isfinite([base per band shift]),2))=Inf;
        tol=o.AbsTol+o.RelTol*abs(band); good=err<=tol;
        y=band; i.errorEstimate=err; i.periodError=ep; i.bandwidthError=eb; i.shiftError=es;
        i.resolvedMask=good; i.points=4*N; i.evaluations=i.evaluations+11*N+4;
        i.fftPeriod=2*P; i.internalStep=P/(2*N); i.bandwidth=2*pi*N/P;
        i.symmetryError=max([i.symmetryError e1 e2 e3 e4]);
        history(end+1,:)=[2*P 4*N max(ep) max(eb) max(es)]; %#ok<AGROW>
        success(good)=success(good)+1; success(~good)=0;
        if all(success>=2), i.stopReason='Two successive period/bandwidth/shift checks passed.'; break; end
        if any(ep>tol), P=2*P; N=2*N; end
        N=2*N;
        i.stopReason='Refinement limit reached.';
    end
    i.history=history;
    i.resolvedMask=i.resolvedMask & success>=2;
    i.converged=all(i.resolvedMask); if i.converged, i.status='converged'; end
end
function [y,symerr,condition]=sample(f,t,P,N,sigma)
    dt=P/N; w=(2*pi/P)*[0:N/2 -N/2+1:-1].';
    v=response_values(f,sigma+1i*w);
    symerr=max(abs(v(2:N/2)-conj(v(end:-1:N/2+2))))/max(realmin,max(abs(v)));
    symerr=max(symerr,abs(imag(v(1)))/max(realmin,abs(v(1))));
    % Even-N Nyquist represents both endpoint half-weights of the integral.
    vn=response_values(f,sigma-1i*w(N/2+1));
    symerr=max(symerr,abs(vn-conj(v(N/2+1)))/max(realmin,abs(vn)));
    if symerr>1e-9, error('ContourRoots:ResponseSymmetry','Transfer is not conjugate symmetric; real-response API cannot discard this imaginary part.'); end
    v(N/2+1)=(v(N/2+1)+vn)/2;
    condition=sum(abs(v))/P;
    z=ifft(v)/dt;
    y=interp1((0:N-1).'*dt,real(z),t,'pchip').*exp(sigma*t);
end
