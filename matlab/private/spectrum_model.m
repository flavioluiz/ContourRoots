function [target,dt,other,analytic,source] = spectrum_model(model,opt)
% Adapt inputs without replacing nonrational terms by rational approximants.
    poles=strcmpi(opt.Mode,'poles'); other=[]; dt=[];
    source='unverified black box'; analytic=false;
    if isa(model,'symfun'), model=formula(model); end
    if isa(model,'sym')
        if ~isscalar(model), error('complex_spectrum:SISO','Use a scalar symbolic expression.'); end
        vars=symvar(model);
        if numel(vars)>1
            error('complex_spectrum:Parameters','Substitute numeric values for all parameters except the complex variable.');
        end
        if isempty(vars), vars=sym('s'); end
        [N,D]=numden(model);
        if poles, a=D; b=N; else, a=N; b=D; end
        if isequal(a,sym(0)) || isequal(b,sym(0))
            error('complex_spectrum:Degenerate','The zero function or a zero denominator has no isolated zero/pole problem.');
        end
        analytic=entire(a) && entire(b);
        if ~analytic && ~opt.AssumeAnalytic
            error('complex_spectrum:AnalyticContract', ...
                ['Cannot establish analytic factors automatically (possible sqrt/log, branch point or inner pole). ', ...
                 'Provide regularized numerator/denominator handles and AssumeAnalytic=true after checking the region.']);
        end
        target=matlabFunction(a,'Vars',vars); other=matlabFunction(b,'Vars',vars);
        dt=matlabFunction(diff(a,vars),'Vars',vars);
        if analytic, source='symbolic entire-function grammar'; else, source='user analytic assertion'; end
        analytic=true;
    elseif isstruct(model) && isscalar(model) && all(isfield(model,{'Numerator','Denominator'}))
        [N,dN,an]=factor(model.Numerator); [D,dD,ad]=factor(model.Denominator);
        if poles, target=D; dt=dD; other=N; else, target=N; dt=dN; other=D; end
        analytic=an && ad;
        if analytic, source='analytic polynomial/symbolic factors'; end
        if opt.AssumeAnalytic, analytic=true; source='user analytic assertion for N and D'; end
        if ~analytic
            % Use the ratio for exploratory mode, so local poles cannot be
            % mistaken for zeros merely because a factor was provided.
            if poles, target=@(z) D(z)./N(z); else, target=@(z) N(z)./D(z); end
            other=[]; dt=[];
        end
    elseif isa(model,'function_handle')
        if poles, target=@(z) 1./model(z); else, target=model; dt=opt.Derivative; end
        if opt.AssumeAnalytic, analytic=true; source='user analytic assertion for target'; end
    elseif isa(model,'tf') || isa(model,'ss') || isa(model,'zpk')
        if ~isequal(size(model),[1 1]) || model.Ts~=0
            error('complex_spectrum:SISO','Only continuous-time scalar SISO models are supported.');
        end
        if isa(model,'tf') || ~hasdelay(model)
            sys=tf(model); [N,D]=tfdata(sys,'v');
            delay=sys.InputDelay+sys.OutputDelay+sys.IODelay;
            pair=struct('Numerator',@(z) polyval(N,z).*exp(-delay*z), ...
                        'Denominator',D);
            asserted=opt; asserted.AssumeAnalytic=true;
            [target,dt,other,analytic,~]=spectrum_model(pair,asserted);
            source='rational LTI factors with external transport delay';
        else
            f=@(z) evalfr(model,z);
            if poles, target=@(z) 1./f(z); else, target=f; end
            source='evalfr of internal-delay LTI model (exploratory)';
        end
    else
        error('complex_spectrum:Input','Use sym, symfun, function handle, N/D struct or continuous SISO tf/ss/zpk.');
    end
end

function [f,df,analytic]=factor(a)
    df=[]; analytic=false;
    if isnumeric(a)
        validateattributes(a,{'numeric'},{'vector','nonempty','finite'});
        if all(a==0), error('complex_spectrum:Degenerate','Factors must not be identically zero.'); end
        f=@(z) polyval(a,z); d=polyder(a); df=@(z) polyval(d,z); analytic=true;
    elseif isa(a,'function_handle')
        f=a;
    elseif isa(a,'sym') || isa(a,'symfun')
        if isa(a,'symfun'), a=formula(a); end
        vars=symvar(a);
        if numel(vars)>1 || ~isscalar(a), error('complex_spectrum:Parameters','Each factor must depend on only one variable.'); end
        if isempty(vars), vars=sym('s'); end
        if isequal(a,sym(0)), error('complex_spectrum:Degenerate','Factors must not be identically zero.'); end
        analytic=entire(a); f=matlabFunction(a,'Vars',vars); df=matlabFunction(diff(a,vars),'Vars',vars);
    else
        error('complex_spectrum:Factor','A factor must be a handle, symbolic scalar or polynomial coefficient vector.');
    end
end

function yes=entire(a)
% Deliberately conservative: sufficient, not necessary, conditions.
    if isempty(symvar(a))
        yes=isfinite(double(a)); return
    end
    if isSymType(a,'variable'), yes=true; return; end
    c=children(a);
    if isSymType(a,'power')
        p=c{2};
        yes=isempty(symvar(p)) && isAlways(p>=0) && isAlways(p==floor(p)) && entire(c{1});
    elseif isSymType(a,'plus | times | exp | sin | cos | sinh | cosh')
        yes=all(cellfun(@entire,c));
    else
        yes=false;
    end
end
