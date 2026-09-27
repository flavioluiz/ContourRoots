function [o,t,a,extra]=matrix_response_options(kind,t,args,plotDefault,sz)
% Matrix extensions are stripped before using the unchanged scalar parser.
    a=args; extra=struct('MaxEvaluations',Inf,'AbsTol',[],'Feedthrough',[],'InitialValue',[]);
    if mod(numel(args),2), error('ContourRoots:ResponseOption','Use name/value pairs.'); end
    keep=true(size(args));
    for k=1:2:numel(args)
        name=args{k};
        if ~(ischar(name)||(isstring(name)&&isscalar(name))), error('ContourRoots:ResponseOption','Option names must be text.'); end
        switch lower(char(name))
            case 'maxevaluations'
                extra.MaxEvaluations=args{k+1}; keep(k:k+1)=false;
            case 'abstol'
                v=args{k+1}; validateattributes(v,{'numeric'},{'vector','real','positive','finite'});
                if ~(isscalar(v)||numel(v)==sz(1)), error('ContourRoots:MatrixShape','AbsTol must be scalar or one value per output.'); end
                extra.AbsTol=v(:); a{k+1}=min(v);
            case {'feedthrough','initialvalue'}
                key='Feedthrough'; if strcmpi(name,'InitialValue'), key='InitialValue'; end
                v=args{k+1}; validateattributes(v,{'numeric'},{'real','finite','size',sz});
                extra.(key)=v; a{k+1}=[];
        end
    end
    a=a(keep);
    try
        [o,t]=response_options(kind,t,a,plotDefault);
    catch err
        if strcmp(err.identifier,'ContourRoots:ResponseBudget')
            error('ContourRoots:MatrixBudget','%s',err.message);
        end
        rethrow(err)
    end
    if isempty(extra.AbsTol), extra.AbsTol=o.AbsTol; end
    if isscalar(extra.AbsTol), extra.AbsTol=repmat(extra.AbsTol,sz(1),1); end
    % Default Inf, like the scalar API (which bounds work only per inversion,
    % through MaxPoints). A finite value is an explicit shared budget.
    validateattributes(extra.MaxEvaluations,{'numeric'},{'scalar','positive','nonnan'});
    if isfinite(extra.MaxEvaluations) && extra.MaxEvaluations~=fix(extra.MaxEvaluations)
        error('ContourRoots:MatrixOption','MaxEvaluations must be a positive integer or Inf.');
    end
    % Reserve space for all channels/diagnostics, not just one scalar output.
    if numel(t)*prod(sz)*(256+80*o.MaxRefinements)>o.MaxMemoryMB*2^20/2
        error('ContourRoots:MatrixBudget','Matrix output/kernel diagnostics exceed the shared memory budget.');
    end
end
