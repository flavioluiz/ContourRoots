function data=matrix_response_prepare(M,t,args)
    if ~isa(M,'ContourRootsModel')||~isscalar(M), error('ContourRoots:MatrixRepresentation','Use one CMIMO/CDYN model.'); end
    [o,t,~,extra]=matrix_response_options('lsim',t,[{'AbsTol',1e-8,'RelTol',1e-6} args],false,M.Size);
    if o.Plot||~isempty(o.Parent), error('ContourRoots:KernelOption','CKERNEL does not plot.'); end
    [models,stats]=matrix_response_context(M,o,extra); bank=cell(M.Size);
    for j=1:M.Size(2), for i=1:M.Size(1)
        co=o; co.AbsTol=extra.AbsTol(i); co.MaxMemoryMB=o.MaxMemoryMB/4;
        bank{i,j}=response_prepare([],t,options_cell(co),models{i,j});
    end, end
    info=stats(); info.converged=true; info.certified=false; info.size=M.Size;
    info.channelInfo=cellfun(@(d) d.info,bank,'UniformOutput',false);
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
