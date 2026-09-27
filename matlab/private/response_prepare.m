function data=response_prepare(G,t,args,m)
% Prepare a fixed, input-independent numerical snapshot. No persistent cache.
    [o,t]=response_options('lsim',t,[{'AbsTol',1e-8,'RelTol',1e-6} args],false);
    if o.Plot || ~isempty(o.Parent)
        error('ContourRoots:KernelOption','CKERNEL does not plot; pass Plot/Parent to CLSIM instead.');
    end
    storage=numel(t)*(256+80*o.MaxRefinements);
    if storage>o.MaxMemoryMB*2^20
        error('ContourRoots:ResponseBudget','Prepared kernels and diagnostics exceed MaxMemoryMB.');
    end
    if nargin<4, m=response_model(G,o); end
    [S,si]=response_kernels(m,t,1,o); R=[]; ri=[];
    if strcmp(o.Interpolation,'foh'), [R,ri]=response_kernels(m,t,2,o); end
    if ~si.converged || (~isempty(ri)&&~ri.converged)
        error('ContourRoots:KernelUnresolved', ...
            'Kernel preparation did not converge. Increase effort, change the method or check the domain; no reusable object was created.');
    end
    count=si.evaluations; if ~isempty(ri), count=count+ri.evaluations; end
    info=struct('converged',true,'certified',false,'method',o.Method, ...
        'evaluations',count,'interpolation',o.Interpolation, ...
        'AbsTol',o.AbsTol,'RelTol',o.RelTol,'step',si,'ramp',ri, ...
        'singularityBound',m.bound,'domainSource',m.domainSource, ...
        'assumptions',{m.assumptions},'feedthrough',m.D,'delay',m.delay);
    % Numeric metadata only: no mutable closure, LTI handle or symbolic engine.
    m=rmfield(m,{'full','regular'});
    data=struct('schema',1,'time',t,'options',o,'model',m, ...
        'step',S,'ramp',R,'stepInfo',si,'rampInfo',ri,'info',info);
    if o.Display
        fprintf('ckernel: %s prepared; %d samples, %d transfer evaluations.\n',o.Interpolation,numel(t),count);
    end
end
