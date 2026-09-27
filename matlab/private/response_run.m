function [y,t,i]=response_run(kind,G,u,t,args,plotDefault,data)
    if nargin<7, data=[]; end
    [o,t]=response_options(kind,t,args,plotDefault);
    if isempty(data), m=response_model(G,o); else, m=data.model; end
    if strcmp(kind,'impulse') && ~m.known && ~o.RegularImpulse
        error('ContourRoots:ResponseImpulse','For an opaque model assert RegularImpulse=true after excluding hidden Dirac terms; supply Feedthrough when present.');
    end
    switch kind
        case {'impulse','inverse'}
            [y,i]=response_kernels(m,t,0,o);
            if m.D~=0, i.singularTerms=struct('time',m.delay,'order',0,'weight',m.D); end
        case 'step'
            [y,i]=response_kernels(m,t,1,o); y=y+m.D*(t>=m.delay);
        case 'lsim'
            u=input_values(u,t);
            convSigma=0;
            if ~isempty(m.bound)&&isfinite(m.bound), convSigma=min(20/t(end),max(0,m.bound)); end
            if strcmp(o.Interpolation,'zoh'), weights=[u(1);diff(u)]; scale=max(1,sum(abs(weights)));
            else
                slopes=diff(u)./diff(t); weights=[slopes(1);diff(slopes);0];
                scale=max(1,abs(u(1))+sum(abs(weights)));
            end
            ko=o; ko.AbsTol=o.AbsTol/(4*scale); ko.RelTol=o.RelTol/(4*scale);
            if strcmp(o.Interpolation,'zoh')
                [S,i]=response_get_kernel(m,t,1,ko,data);
                y=response_convolve(S,weights,t,convSigma,o);
                err=error_conv(i.errorEstimate,abs(weights),t,o);
                resolved=all(i.resolvedMask);
            else
                [R,i]=response_get_kernel(m,t,2,ko,data);
                y=response_convolve(R,weights,t,convSigma,o);
                err=error_conv(i.errorEstimate,abs(weights),t,o);
                resolved=all(i.resolvedMask);
                if u(1)~=0
                    [S,j]=response_get_kernel(m,t,1,ko,data); y=y+u(1)*S;
                    err=err+abs(u(1))*j.errorEstimate;
                    resolved=resolved&&all(j.resolvedMask); i.evaluations=i.evaluations+j.evaluations;
                    i.stepKernelInfo=j;
                end
            end
            y=y+m.D*held(u,t,t-m.delay,o.Interpolation);
            i.kernelErrorEstimate=i.errorEstimate; i.errorEstimate=err;
            i.resolvedMask=resolved & isfinite(y) & err<=o.AbsTol+o.RelTol*abs(y);
            i.converged=all(i.resolvedMask); i.status='unresolved'; if i.converged, i.status='converged'; end
            if ~i.converged, i.stopReason='Kernel inversion or propagated convolution error did not meet tolerance.'; end
            i.interpolation=o.Interpolation;
            i.inputErrorEstimate=NaN;
            i.assumptions{end+1}='Input is the chosen interpolant of supplied samples; unsampled input error is not estimated.';
            i.kernelConvention='Integrated step/ramp basis; linear convolution, no dt factor.';
            i.kernelReused=~isempty(data); i.kernelPreparationEvaluations=0;
            if ~isempty(data)
                i.kernelPreparationEvaluations=data.info.evaluations;
                if i.converged
                    i.stopReason='Stored kernel errors propagated through this input met output tolerances; no inversion performed.';
                end
                i.assumptions{end+1}='Fixed prepared kernels: no new inversion; kernel history describes preparation, output errors are rechecked for this input.';
            end
    end
    i.singularityBound=m.bound; i.domainSource=m.domainSource;
    i.assumptions=[m.assumptions i.assumptions];
    i.assumptions{end+1}='Convergence is conditional on the asserted inversion domain, not a global stability or error certificate.';
    if ~i.converged && o.Warn
        warning('ContourRoots:ResponseUnresolved','Response unresolved: %s Inspect info.resolvedMask and refinement history.',i.stopReason);
    end
    if o.Display
        fprintf('%s: %s; %d/%d samples resolved, %d evaluations.\n',kind,i.status,sum(i.resolvedMask),numel(t),i.evaluations);
    end
    if o.Plot, response_plot(kind,t,y,u,i,o); end
end
function u=input_values(u,t)
    if isa(u,'function_handle')
        try
            v=u(t);
        catch
            v=arrayfun(u,t);
        end
        if isscalar(v), v=arrayfun(u,t); end
        u=v;
    end
    validateattributes(u,{'numeric'},{'vector','real','finite','numel',numel(t)});
    u=double(u(:));
end
function y=held(u,t,q,method)
    y=zeros(size(q)); active=q>=0 & q<=t(end);
    if strcmp(method,'zoh'), method='previous'; else, method='linear'; end
    y(active)=interp1(t,u,q(active),method);
end
function e=error_conv(e,w,t,o)
    if any(~isfinite(e))
        % Avoid Inf*0 and do not disguise unresolved kernels.
        e=inf(size(t)); return
    end
    e=max(0,response_convolve(e,w,t,0,o));
end
