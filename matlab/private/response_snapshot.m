function data=response_snapshot(m,t,o,S,si,R,ri)
% Preserve the schema-1 numerical snapshot for either preparation engine.
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
