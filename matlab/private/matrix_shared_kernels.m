function bank=matrix_shared_kernels(models,block,t,o,extra)
% One refinement decision for every channel AND both integrated kernels.
% Evaluations return undelayed full pages; delays are exact time shifts and
% direct terms are removed before inversion. Snapshot schema stays scalar.
    nc=numel(models); orders=1; if strcmp(o.Interpolation,'foh'), orders=1:2; end
    ns=nc*numel(orders); nt=numel(t); tau=zeros(nt,nc); active=false(nt,nc);
    tol=repmat(extra.AbsTol,numel(models)/numel(extra.AbsTol),1).';
    tol=repmat(tol,1,numel(orders)); bounds=0; D=zeros(1,nc);
    infos=cell(1,ns); values=zeros(nt,ns); enabled=false(1,nc);
    for channel=1:nc
        m=models{channel}; tau(:,channel)=t-m.delay;
        tau(abs(tau(:,channel))<=8*eps*max(1,abs(t)),channel)=0;
        active(:,channel)=tau(:,channel)>0 & ~m.zero; enabled(channel)=any(active(:,channel)); D(channel)=m.D;
        if ~isempty(m.bound), bounds=max(bounds,m.bound); end
        for power=orders
            signal=channel+(power-1)*nc; infos{signal}=response_info(t,o);
            for key={'errorEstimate','periodError','bandwidthError','shiftError'}
                infos{signal}.(key{1})(:)=0;
            end
            infos{signal}.resolvedMask(:)=true;
            infos{signal}.stopReason='Exact zero regular response on this interval.';
        end
    end
    lineOptions=o; lineOptions.AssumeStable=false;
    if any(enabled)&&~isempty(o.Abscissa)&&o.Abscissa<=0
        error('ContourRoots:ResponseDomain','Step/ramp inversion requires a positive Abscissa.');
    end
    if any(enabled)
        if strcmp(o.Method,'fft'), fft_bank(); else, dehoog_bank(); end
    end
    bank=cell(size(models));
    for channel=1:nc
        co=o; co.AbsTol=extra.AbsTol(mod(channel-1,size(models,1))+1); co.MaxMemoryMB=o.MaxMemoryMB/4;
        si=finish(infos{channel}); R=[]; ri=[];
        if numel(orders)==2, R=values(:,channel+nc); ri=finish(infos{channel+nc}); end
        bank{channel}=response_snapshot(models{channel},t,co,values(:,channel),si,R,ri);
    end

    function fft_bank()
        horizon=max(tau(active)); dt=min(diff(t));
        P=4*max(horizon,dt); N=2^nextpow2(max(1024,4*P/dt));
        mask=repmat(active,1,numel(orders)); success=zeros(nt,ns);
        for j=1:ns, infos{j}.resolvedMask=~mask(:,j); end
        for level=1:o.MaxRefinements
            sigma=response_line(lineOptions,bounds,P);
            if max(0,sigma+.5/P)*horizon>log(1e12)
                stop('Exponential amplification exceeds 1e12; use shorter time blocks.'); break
            end
            % Full matrix samples plus scalar FFT temporaries, in addition to
            % the separately reserved cache and stored channel diagnostics.
            if 4*N>o.MaxPoints || 4*N*(32*nc+128)>o.MaxMemoryMB*2^20/4
                stop('Shared FFT point/memory budget exhausted.'); break
            end
            [base,e1,c1]=fft_sample(P,N,sigma);
            [per,e2,c2]=fft_sample(2*P,2*N,sigma);
            [band,e3,c3]=fft_sample(2*P,4*N,sigma);
            [shift,e4,c4]=fft_sample(2*P,4*N,sigma+.5/P);
            ep=abs(per-base); eb=abs(band-per); es=abs(shift-band);
            floorError=eps*exp(max(0,sigma+.5/P)*horizon)*max([c1;c2;c3;c4],[],1)*20;
            err=max(max(ep,eb),es)+floorError;
            err(~isfinite(base)|~isfinite(per)|~isfinite(band)|~isfinite(shift))=Inf;
            err(~mask)=0; good=err<=tol+o.RelTol*abs(band);
            success(good)=success(good)+1; success(~good)=0; values=band;
            for j=1:ns
                q=infos{j}; q.abscissa=sigma; q.amplification=exp(max(0,sigma)*horizon);
                q.errorEstimate=err(:,j); q.periodError=ep(:,j); q.bandwidthError=eb(:,j); q.shiftError=es(:,j);
                q.resolvedMask=~mask(:,j)|(good(:,j)&success(:,j)>=2);
                q.points=4*N; q.evaluations=q.evaluations+(11*N+4)*any(mask(:,j));
                q.fftPeriod=2*P; q.internalStep=P/(2*N); q.bandwidth=2*pi*N/P;
                q.symmetryError=max([q.symmetryError e1(j) e2(j) e3(j) e4(j)]);
                q.history(end+1,:)=[2*P 4*N max(ep(:,j)) max(eb(:,j)) max(es(:,j))];
                infos{j}=q;
            end
            if all(success(mask)>=2)
                stop('Two successive common-grid period/bandwidth/shift checks passed.'); return
            end
            threshold=tol+o.RelTol*abs(band);
            if any(ep(mask)>threshold(mask))
                P=2*P; N=2*N;
            end
            N=2*N; stop('Shared refinement limit reached.');
        end
    end

    function [y,sym,condition]=fft_sample(P,N,sigma)
        dt=P/N; w=(2*pi/P)*[0:N/2 -N/2+1:-1].';
        s=sigma+1i*w; pages=block(s); negative=block(sigma-1i*w(N/2+1));
        y=zeros(nt,ns); sym=zeros(1,ns); condition=sym;
        for c=find(enabled)
            for p=orders
                j=c+(p-1)*nc;
                v=(pages(:,c)-D(c))./s.^p;
                vn=(negative(c)-D(c))/(sigma-1i*w(N/2+1))^p;
                sym(j)=max(abs(v(2:N/2)-conj(v(end:-1:N/2+2))))/max(realmin,max(abs(v)));
                sym(j)=max([sym(j),abs(imag(v(1)))/max(realmin,abs(v(1))), ...
                    abs(vn-conj(v(N/2+1)))/max(realmin,abs(vn))]);
                check_symmetry(sym(j)); v(N/2+1)=(v(N/2+1)+vn)/2;
                condition(j)=sum(abs(v))/P; z=ifft(v)/dt; at=active(:,c); tt=tau(at,c);
                y(at,j)=interp1((0:N-1).'*dt,real(z),tt,'pchip').*exp(sigma*tt);
            end
        end
    end

    function dehoog_bank()
        degrees=[10 20 40 80 160 320];
        for k=1:nt
            channels=find(active(k,:)); if isempty(channels), continue; end
            ix=reshape(channels(:)+(orders-1)*nc,1,[]);
            P=4*max(tau(k,channels)); sigma=response_line(lineOptions,bounds,P);
            for j=ix, infos{j}.resolvedMask(k)=false; end
            if (sigma+.5/P)*max(tau(k,channels))>log(1e12), continue; end
            previous=NaN(1,ns); success=zeros(1,ns);
            for level=1:min(o.MaxRefinements,numel(degrees))
                M=degrees(level); L=ceil(1.25*M);
                if 4*M+1>o.MaxPoints || (4*M+1)^2*64+(4*M+1)*32*nc>o.MaxMemoryMB*2^20/4, break; end
                [a,ea,ca]=dehoog_sample(k,channels,P,sigma,M);
                [b,eb,cb]=dehoog_sample(k,channels,1.25*P,sigma,L);
                [c,ec,cc]=dehoog_sample(k,channels,1.25*P,sigma+.5/P,L);
                ep=abs(b-a); es=abs(c-b); ed=abs(b-previous);
                amp=exp(max(0,sigma+.5/P)*max(tau(k,channels)));
                err=max([ep;es;ed],[],1)+eps*amp*max([ca;cb;cc],[],1)*100;
                err(~isfinite(a)|~isfinite(b)|~isfinite(c)|~isfinite(previous)|~isfinite(err))=Inf;
                good=err<=tol+o.RelTol*abs(b); success(good)=success(good)+1; success(~good)=0;
                for j=ix
                    values(k,j)=b(j); q=infos{j}; q.errorEstimate(k)=err(j);
                    q.periodError(k)=ep(j); q.bandwidthError(k)=ed(j); q.shiftError(k)=es(j);
                    q.evaluations=q.evaluations+4*M+8*L+6; q.points=max(q.points,2*L+1);
                    if isnan(q.abscissa), q.abscissa=sigma; else, q.abscissa=max(q.abscissa,sigma); end
                    q.amplification=max(q.amplification,amp); q.symmetryError=max([q.symmetryError ea(j) eb(j) ec(j)]);
                    q.history(end+1,:)=[tau(k,mod(j-1,nc)+1) M ep(j) ed(j) es(j)];
                    q.resolvedMask(k)=success(j)>=2&&M>=40; infos{j}=q;
                end
                if all(success(ix)>=2)&&M>=40, break; end
                previous=b;
            end
        end
        stop('Common-grid degree/period/shift checks; inspect per-channel resolvedMask.');
    end

    function [v,sym,condition]=dehoog_sample(k,channels,P,sigma,M)
        s=sigma+1i*(0:2*M).'*2*pi/P; pages=block(s); negative=block(conj(s));
        v=zeros(1,ns); sym=v; condition=v;
        for c=channels
            for p=orders
                j=c+(p-1)*nc; fp=(pages(:,c)-D(c))./s.^p; neg=(negative(:,c)-D(c))./conj(s).^p;
                sym(j)=max(abs(fp-conj(neg)))/max(realmin,max(abs(fp))); check_symmetry(sym(j));
                condition(j)=2*sum(abs(fp))/P;
                v(j)=response_dehoog_value(fp,tau(k,c),P,sigma,M);
            end
        end
    end
    function stop(reason)
        for j=1:ns, if enabled(mod(j-1,nc)+1), infos{j}.stopReason=reason; end, end
    end
end
function q=finish(q)
    q.converged=all(q.resolvedMask); if q.converged, q.status='converged'; end
end
function check_symmetry(e)
    if e>1e-9, error('ContourRoots:ResponseSymmetry','Transfer is not conjugate symmetric.'); end
end
