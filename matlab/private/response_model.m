function m=response_model(G,o)
% Dedicated transfer adapter; deliberately independent of spectrum_model.
    m=struct('full',[],'regular',[],'bound',[],'domainSource','user assertion', ...
        'D',0,'delay',0,'initial',NaN,'known',false,'zero',false, ...
        'assumptions',{{'Causal real SISO response, zero prehistory, consistent branches.'}});
    if isa(G,'tf') || isa(G,'ss') || isa(G,'zpk')
        if ~isequal(size(G),[1 1]) || G.Ts~=0
            error('ContourRoots:ResponseModel','Use a continuous-time SISO model.');
        end
        if isa(G,'tf') || ~hasdelay(G)
            q=tf(G); [n,d]=tfdata(q,'v'); m=rational(n,d,m);
            m.delay=q.InputDelay+q.OutputDelay+q.IODelay;
            if m.delay<0, error('ContourRoots:ResponseModel','Negative transport delay is noncausal.'); end
        else
            m.full=@(s) evalfr(G,s);
        end
    elseif isnumeric(G)
        validateattributes(G,{'numeric'},{'scalar','real','finite'}); m=rational(G,1,m);
    elseif isa(G,'function_handle')
        m.full=G;
    elseif isa(G,'sym') || isa(G,'symfun')
        if isa(G,'symfun'), G=formula(G); end
        checkvars(G); m.full=tohandle(G);
    elseif isstruct(G) && isscalar(G) && all(isfield(G,{'Numerator','Denominator'}))
        n=G.Numerator; d=G.Denominator;
        vars=union(symbols(n),symbols(d));
        if numel(vars)>1, error('ContourRoots:ResponseModel','N and D must use the same single symbolic variable.'); end
        if isnumeric(n)&&isnumeric(d), m=rational(n,d,m);
        else
            nf=tohandle(n); df=tohandle(d);
            m.full=@(s) response_values(nf,s)./response_values(df,s);
        end
    else
        error('ContourRoots:ResponseModel','Use a handle, ndpair, scalar gain, symbolic scalar or continuous SISO LTI model.');
    end
    if ~isempty(o.Feedthrough)
        if m.known && abs(o.Feedthrough-m.D)>1e-12*(1+abs(m.D))
            error('ContourRoots:ResponseModel','Feedthrough contradicts exact model metadata.');
        end
        m.D=o.Feedthrough;
    elseif ~m.known
        m.assumptions{end+1}='Feedthrough=0 assumed: supplied transfer has a locally integrable ordinary impulse response, with no hidden Dirac terms.';
    end
    if ~m.known
        f=m.full; D=m.D; m.regular=@(s) response_values(f,s)-D;
    end
    if ~isempty(o.InitialValue)
        if m.known && abs(o.InitialValue-m.initial)>1e-12*(1+abs(m.initial))
            error('ContourRoots:ResponseModel','InitialValue contradicts the analytical ordinary impulse right limit.');
        end
        m.initial=o.InitialValue;
    end
    if ~isempty(o.SingularityBound)
        if ~isempty(m.bound) && o.SingularityBound<m.bound
            error('ContourRoots:ResponseDomain','SingularityBound contradicts rational denominator poles.');
        end
        m.bound=o.SingularityBound; m.domainSource='user half-plane assertion';
    elseif o.AssumeStable
        if ~isempty(m.bound) && m.bound>=0
            error('ContourRoots:ResponseDomain','AssumeStable contradicts rational denominator poles.');
        end
        m.bound=0; m.domainSource='user integrable-impulse assertion';
    end
    if isempty(m.bound) && isempty(o.Abscissa)
        error('ContourRoots:ResponseDomain','Provide SingularityBound or Abscissa; a bounded pole search cannot establish the inversion half-plane.');
    end
    if ~isempty(o.Abscissa) && ~isempty(m.bound) && o.Abscissa<=m.bound
        error('ContourRoots:ResponseDomain','Abscissa must lie to the right of the model bound.');
    end
end
function m=rational(n,d,m)
    validateattributes(n,{'numeric'},{'vector','nonempty','real','finite'});
    validateattributes(d,{'numeric'},{'vector','nonempty','real','finite'});
    n=trim(n); d=trim(d);
    if all(d==0), error('ContourRoots:ResponseModel','Zero denominator.'); end
    if numel(n)>numel(d), error('ContourRoots:ResponseModel','Improper transfers with impulse derivatives are not supported.'); end
    if numel(n)==numel(d), m.D=n(1)/d(1); nr=n-m.D*d; else, nr=n; end
    nr=trim(nr); m.zero=all(nr==0);
    m.initial=0;
    if numel(nr)==numel(d)-1, m.initial=nr(1)/d(1); end
    m.full=@(s) polyval(n,s)./polyval(d,s);
    m.regular=@(s) polyval(nr,s)./polyval(d,s);
    p=roots(d); if isempty(p), m.bound=-Inf; else, m.bound=max(real(p)); end
    m.known=true; m.domainSource='rational denominator (conservative before cancellations)';
end
function a=trim(a)
    a=double(a(:).'); j=find(a~=0,1); if isempty(j), a=0; else, a=a(j:end); end
end
function f=tohandle(a)
    if isa(a,'function_handle')
        f=a;
    elseif isnumeric(a)
        validateattributes(a,{'numeric'},{'vector','nonempty','real','finite'}); f=@(s) polyval(a,s);
    elseif isa(a,'sym') || isa(a,'symfun')
        if isa(a,'symfun'), a=formula(a); end
        checkvars(a); v=symvar(a); if isempty(v), v=sym('s'); end
        f=matlabFunction(a,'Vars',v);
    else
        error('ContourRoots:ResponseModel','Invalid numerator/denominator factor.');
    end
end
function checkvars(a)
    if ~isscalar(a)||numel(symvar(a))>1
        error('ContourRoots:ResponseModel','Substitute parameters; use one scalar complex variable.');
    end
end
function v=symbols(a)
    v={};
    if isa(a,'sym')||isa(a,'symfun')
        if isa(a,'symfun'), a=formula(a); end
        checkvars(a); x=symvar(a); if ~isempty(x), v={char(x)}; end
    end
end
