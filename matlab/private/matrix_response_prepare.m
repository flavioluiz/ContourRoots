function data=matrix_response_prepare(M,t,args)
    if ~isa(M,'ContourRootsModel')||~isscalar(M), error('ContourRoots:MatrixRepresentation','Use one CMIMO/CDYN model.'); end
    [shared,args]=shared_option(args);
    [o,t,~,extra]=matrix_response_options('lsim',t,[{'AbsTol',1e-8,'RelTol',1e-6} args],false,M.Size);
    if o.Plot||~isempty(o.Parent), error('ContourRoots:KernelOption','CKERNEL does not plot.'); end
    if shared && strcmp(o.Method,'quadrature')
        error('ContourRoots:KernelOption','SharedGrid supports fft and dehoog; use SharedGrid=false for quadrature.');
    end
    [models,stats,block]=matrix_response_context(M,o,extra); bank=cell(M.Size);
    if shared
        bank=matrix_shared_kernels(models,block,t,o,extra);
    else
        for j=1:M.Size(2), for i=1:M.Size(1)
            co=o; co.AbsTol=extra.AbsTol(i); co.MaxMemoryMB=o.MaxMemoryMB/4;
            bank{i,j}=response_prepare([],t,options_cell(co),models{i,j});
        end, end
    end
    info=stats(); info.converged=true; info.certified=false; info.size=M.Size;
    info.channelInfo=cellfun(@(d) d.info,bank,'UniformOutput',false);
    info.sharedGrid=shared;
    info.AbsTol=extra.AbsTol; info.RelTol=o.RelTol;
    info.interpolation=o.Interpolation; info.method=o.Method;
    info.channelMetadata=struct('InputNames',{M.Options.InputNames},'OutputNames',{M.Options.OutputNames}, ...
        'InputUnits',{M.Options.InputUnits},'OutputUnits',{M.Options.OutputUnits});
    info.assumptions={'One numeric snapshot for all channels; common user inversion-domain assertions when required.'};
    data=struct('schema',1,'time',t,'size',M.Size,'interpolation',o.Interpolation, ...
        'options',o,'channels',{bank},'info',info);
end
function args=options_cell(o)
    names=fieldnames(o); args=reshape([names struct2cell(o)].',1,[]);
end

function [shared,args]=shared_option(args)
    shared=false;
    if mod(numel(args),2), error('ContourRoots:ResponseOption','Use name/value pairs.'); end
    keep=true(size(args));
    for k=1:2:numel(args)
        name=args{k};
        if (ischar(name)||(isstring(name)&&isscalar(name)))&&strcmpi(name,'SharedGrid')
            value=args{k+1};
            if ~(isscalar(value)&&(islogical(value)||isnumeric(value))&&isreal(value)&&any(value==[0 1]))
                error('ContourRoots:KernelOption','SharedGrid must be logical.');
            end
            shared=logical(value); keep(k:k+1)=false;
        end
    end
    args=args(keep);
end
