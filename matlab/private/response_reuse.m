function [y,t,i]=response_reuse(data,u,t,args,plotDefault)
% Restrict reuse to exactly the numerical snapshot; never reinvert silently.
    if data.schema~=1
        error('ContourRoots:KernelModel','Unsupported kernel snapshot; prepare it again with CKERNEL.');
    end
    allowed={'AbsTol','RelTol','MaxPoints','MaxMemoryMB','Interpolation','Plot','Parent','Warn','Display'};
    if mod(numel(args),2), error('ContourRoots:KernelOption','Use name/value pairs.'); end
    for j=1:2:numel(args)
        name=args{j};
        if ~(ischar(name)||(isstring(name)&&isscalar(name))) || ~any(strcmpi(name,allowed))
            error('ContourRoots:KernelOption','Only output tolerances, convolution budgets and presentation can change; prepare a new kernel for model/inversion changes.');
        end
    end
    defaults={'Method',data.options.Method,'Interpolation',data.options.Interpolation, ...
        'MaxPoints',data.options.MaxPoints,'MaxMemoryMB',data.options.MaxMemoryMB};
    args=[defaults args]; [o,t]=response_options('lsim',t,args,plotDefault);
    if ~isequal(t,data.time)
        error('ContourRoots:KernelGrid','Times must exactly match K.Time; prepare a new kernel for another grid.');
    end
    if ~strcmp(o.Interpolation,data.options.Interpolation)
        error('ContourRoots:KernelInterpolation','Interpolation must match K.Interpolation; prepare a new kernel for another hold.');
    end
    [y,t,i]=response_run('lsim',[],u,t,args,plotDefault,data);
end
