function [o,t] = response_options(kind,t,args,plotDefault)
    validateattributes(t,{'numeric'},{'real','vector','nonempty','finite','nonnegative'});
    t=double(t(:));
    if any(diff(t)<=0), error('ContourRoots:ResponseTime','Times must be strictly increasing.'); end
    o=struct('Method','fft','SingularityBound',[],'Abscissa',[], ...
        'AssumeStable',false,'AbsTol',1e-6,'RelTol',1e-4, ...
        'MaxRefinements',8,'MaxPoints',2^20,'MaxMemoryMB',256, ...
        'Interpolation','foh','Feedthrough',[],'RegularImpulse',false, ...
        'InitialValue',[],'Plot',plotDefault,'Parent',[], 'Warn',true,'Display',false);
    if mod(numel(args),2), error('ContourRoots:ResponseOption','Use name/value pairs.'); end
    names=fieldnames(o);
    for j=1:2:numel(args)
        key=args{j};
        if ~(ischar(key) || (isstring(key)&&isscalar(key)))
            error('ContourRoots:ResponseOption','Option names must be text.');
        end
        ix=find(strcmpi(key,names),1);
        if isempty(ix), error('ContourRoots:ResponseOption','Unknown option %s.',key); end
        o.(names{ix})=args{j+1};
    end
    o.Method=validatestring(o.Method,{'fft','dehoog','quadrature'});
    o.Interpolation=validatestring(o.Interpolation,{'foh','zoh'});
    for key={'AssumeStable','RegularImpulse','Plot','Warn','Display'}
        v=o.(key{1});
        if ~(isscalar(v) && (islogical(v)||isnumeric(v)) && isreal(v) && any(v==[0 1]))
            error('ContourRoots:ResponseOption','%s must be logical.',key{1});
        end
        o.(key{1})=logical(v);
    end
    for key={'AbsTol','RelTol','MaxMemoryMB'}
        validateattributes(o.(key{1}),{'numeric'},{'scalar','real','finite','positive'});
    end
    for key={'MaxRefinements','MaxPoints'}
        validateattributes(o.(key{1}),{'numeric'},{'scalar','real','finite','integer','positive'});
    end
    for key={'SingularityBound','Abscissa','Feedthrough','InitialValue'}
        if ~isempty(o.(key{1}))
            validateattributes(o.(key{1}),{'numeric'},{'scalar','real','finite'});
        end
    end
    if ~isempty(o.Parent) && ~(isscalar(o.Parent)&&isgraphics(o.Parent,'axes'))
        error('ContourRoots:ResponseOption','Parent must be an axes.');
    end
    if ~isempty(o.Abscissa) && ~isempty(o.SingularityBound) && o.Abscissa<=o.SingularityBound
        error('ContourRoots:ResponseDomain','Abscissa must exceed SingularityBound.');
    end
    uniform=numel(t)<3 || max(abs(diff(t)-mean(diff(t))))<=1e-9*max(diff(t));
    if (strcmp(o.Method,'fft') || strcmp(kind,'lsim')) && ~uniform
        error('ContourRoots:ResponseTime','FFT/CLSIM needs a uniform time grid.');
    end
    if strcmp(kind,'lsim') && (numel(t)<2 || t(1)~=0)
        error('ContourRoots:ResponseTime','CLSIM requires at least two times, starting at zero.');
    end
    if numel(t)*128>o.MaxMemoryMB*2^20
        error('ContourRoots:ResponseBudget','Requested output/diagnostic arrays exceed MaxMemoryMB.');
    end
end
