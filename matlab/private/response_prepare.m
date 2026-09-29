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
    data=response_snapshot(m,t,o,S,si,R,ri);
end
