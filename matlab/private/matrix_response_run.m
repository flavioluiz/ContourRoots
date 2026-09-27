function [y,t,info]=matrix_response_run(kind,M,u,t,args,plotDefault,data)
    if nargin<7, data=[]; end
    reused=~isempty(data);
    if reused
        if data.schema~=1, error('ContourRoots:KernelModel','Unsupported matrix kernel schema.'); end
        allowed={'AbsTol','RelTol','MaxPoints','MaxMemoryMB','Interpolation','Plot','Parent','Warn','Display'};
        if mod(numel(args),2), error('ContourRoots:KernelOption','Use name/value pairs.'); end
        for k=1:2:numel(args)
            if ~any(strcmpi(args{k},allowed)), error('ContourRoots:KernelOption','Only output tolerances, convolution budgets and presentation may change on reuse.'); end
        end
        sz=data.size;
        args=[{'Method',data.options.Method,'Interpolation',data.interpolation, ...
            'MaxPoints',data.options.MaxPoints,'MaxMemoryMB',data.options.MaxMemoryMB} args];
    else
        if ~isscalar(M), error('ContourRoots:MatrixRepresentation','Use one matrix model.'); end
        sz=M.Size;
    end
    [o,t,~,extra]=matrix_response_options(kind,t,args,plotDefault,sz);
    if reused
        if ~isequal(t,data.time), error('ContourRoots:KernelGrid','Times must exactly match K.Time.'); end
        if ~strcmp(o.Interpolation,data.interpolation), error('ContourRoots:KernelInterpolation','Interpolation must match K.Interpolation.'); end
    else
        [models,stats]=matrix_response_context(M,o,extra);
    end
    if strcmp(kind,'lsim')
        if isa(u,'function_handle'), u=u(t); end
        validateattributes(u,{'numeric'},{'real','finite','size',[numel(t) sz(2)]});
        y=zeros(numel(t),sz(1)); err=y; mask=true(size(y)); magnitude=y;
    else
        y=zeros(numel(t),sz(1),sz(2)); err=y; mask=true(size(y));
    end
    channelInfo=cell(sz); terms=struct('output',{},'input',{},'time',{},'order',{},'weight',{});
    for j=1:sz(2), for i=1:sz(1)
        co=o; co.Plot=false; co.Display=false; co.Warn=false; co.Parent=[];
        co.AbsTol=extra.AbsTol(i); co.MaxMemoryMB=o.MaxMemoryMB/4;
        if strcmp(kind,'lsim'), co.AbsTol=co.AbsTol/sz(2); end
        ca=options_cell(co);
        if reused
            % Use the same scalar convolution with stored numerical metadata.
            [v,~,ci]=response_run(kind,[],u(:,j),t,ca,false,data.channels{i,j});
        else
            ui=[]; if strcmp(kind,'lsim'), ui=u(:,j); end
            [v,~,ci]=response_run(kind,[],ui,t,ca,false,[],models{i,j});
        end
        channelInfo{i,j}=ci;
        for k=1:numel(ci.singularTerms)
            q=ci.singularTerms(k); terms(end+1)=struct('output',i,'input',j,'time',q.time,'order',q.order,'weight',q.weight); %#ok<AGROW>
        end
        if strcmp(kind,'lsim')
            y(:,i)=y(:,i)+v; err(:,i)=err(:,i)+ci.errorEstimate;
            mask(:,i)=mask(:,i)&ci.resolvedMask; magnitude(:,i)=magnitude(:,i)+abs(v);
        else
            y(:,i,j)=v; err(:,i,j)=ci.errorEstimate; mask(:,i,j)=ci.resolvedMask;
        end
    end, end
    if strcmp(kind,'lsim')
        err=err+8*eps*sz(2)*magnitude; % include floating-point summation, especially cancellation
    end
    mask=mask & isfinite(y) & err<=reshape(extra.AbsTol,1,[])+o.RelTol*abs(y);
    if reused
        info=struct('evaluations',0,'factorEvaluations',0,'factorizations',0,'linearSolves',0,'rhsColumns',0,'cacheHits',0);
        info.kernelPreparationEvaluations=data.info.evaluations;
    else
        info=stats(); info.kernelPreparationEvaluations=0;
    end
    info.kind=['matrix_' kind]; info.size=sz; info.converged=all(mask(:)); info.certified=false;
    info.status='unresolved'; if info.converged, info.status='converged'; end
    info.errorEstimate=err; info.resolvedMask=mask; info.channelInfo=channelInfo;
    info.kernelReused=reused; info.singularTerms=terms; info.time=t;
    info.assumptions={'Causal real zero-state response; stated branches/domain and per-channel direct terms.', ...
        'Convergence is conditional, not a global spectrum or error certificate.'};
    if strcmp(kind,'lsim'), info.interpolation=o.Interpolation; info.inputErrorEstimate=NaN; end
    info.stopReason='Channel errors summed per output and compared with the total response (including cancellation).';
    if ~info.converged&&o.Warn, warning('ContourRoots:ResponseUnresolved','Matrix response unresolved; inspect channelInfo and resolvedMask.'); end
    if o.Display, fprintf('%s: %s, %d matrix-node evaluations.\n',info.kind,info.status,info.evaluations); end
    if o.Plot
        ax=o.Parent; if isempty(ax), figure; ax=axes; end
        plot(ax,t,reshape(y,numel(t),[])); grid(ax,'on'); xlabel(ax,'Time'); ylabel(ax,'Response');
        title(ax,strrep(info.kind,'_',' '));
    end
end
function args=options_cell(o)
    names=fieldnames(o); args=reshape([names struct2cell(o)].',1,[]);
end
