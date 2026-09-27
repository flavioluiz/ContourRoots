function [kind,sz,dims,f,o]=matrix_model(kind,f,varargin)
% Constructor validation only: never evaluate opaque functions here.
    if ~any(strcmp(kind,{'transfer','structured'}))
        error('ContourRoots:MatrixRepresentation','Construct models with CMIMO or CDYN.');
    end
    o=struct('Size',[],'Dimensions',[],'Batched',false,'DomainCheck',[], ...
        'Feedthrough',[],'Delay',[],'InitialValue',[], ...
        'InputNames',{{}},'OutputNames',{{}},'InputUnits',{{}},'OutputUnits',{{}}, ...
        'DomainDescription','User-declared domain; guards do not certify analyticity.');
    if mod(numel(varargin),2), error('ContourRoots:MatrixOption','Use name/value pairs.'); end
    names=fieldnames(o);
    for k=1:2:numel(varargin)
        key=varargin{k}; ix=find(strcmpi(key,names),1);
        if isempty(ix), error('ContourRoots:MatrixOption','Unknown matrix option.'); end
        o.(names{ix})=varargin{k+1};
    end
    validateattributes(o.Batched,{'logical','numeric'},{'scalar','binary'}); o.Batched=logical(o.Batched);
    if ~(ischar(o.DomainDescription)&&isrow(o.DomainDescription) || isstring(o.DomainDescription)&&isscalar(o.DomainDescription))
        error('ContourRoots:MatrixDomain','DomainDescription must be scalar text.');
    end
    o.DomainDescription=char(o.DomainDescription);
    if ~isempty(o.DomainCheck)&&~isa(o.DomainCheck,'function_handle')
        error('ContourRoots:MatrixDomain','DomainCheck must be a scalar-node logical guard.');
    end
    dims=[];
    if strcmp(kind,'transfer')
        if ~isempty(o.Dimensions), error('ContourRoots:MatrixOption','Use Size for CMIMO.'); end
        if isa(f,'symfun'), f=formula(f); end
        if isa(f,'tf')||isa(f,'ss')||isa(f,'zpk')
            if f.Ts~=0, error('ContourRoots:MatrixRepresentation','Only continuous-time LTI models are supported.'); end
            descriptor_check(f);
            q=cell(size(f)); for j=1:size(f,2), for i=1:size(f,1), q{i,j}=f(i,j); end, end; f=q;
        end
        if iscell(f)
            if isempty(f)||~ismatrix(f), error('ContourRoots:MatrixShape','Channels must be a nonempty rectangular cell array.'); end
            kind='channels'; sz=size(f);
            vars={};
            for k=1:numel(f), checkchannel(f{k}); vars=union(vars,symbol_names(f{k})); end
            if numel(vars)>1, error('ContourRoots:MatrixRepresentation','All symbolic channels must use one common variable.'); end
        elseif isnumeric(f)||isa(f,'sym')
            checkmatrix(f); sz=size(f); f=convert(f);
        elseif isa(f,'function_handle')
            sz=o.Size; checkdims(sz,2);
        else
            error('ContourRoots:MatrixRepresentation','Use a matrix handle, numeric/symbolic matrix, LTI model or channel cell array.');
        end
        if ~isempty(o.Size)&&~isequal(double(o.Size(:).'),double(sz))
            error('ContourRoots:MatrixShape','Size contradicts the model dimensions.');
        end
    else
        if ~isempty(o.Size), error('ContourRoots:MatrixOption','Use Dimensions for CDYN.'); end
        dims=o.Dimensions;
        if isempty(dims)
            dims=nan(1,3);
            for k=1:4
                v=f{k}; if isa(v,'symfun'), v=formula(v); end
                if isnumeric(v)||isa(v,'sym')
                    checkmatrix(v); a=size(v);
                    switch k
                        case 1, if a(1)~=a(2), error('ContourRoots:MatrixShape','H must be square.'); end; dims=merge(dims,[a(1) NaN NaN]);
                        case 2, dims=merge(dims,[a(1) NaN a(2)]);
                        case 3, dims=merge(dims,[a(2) a(1) NaN]);
                        case 4, dims=merge(dims,[NaN a(1) a(2)]);
                    end
                end
            end
        end
        checkdims(dims,3); dims=double(dims(:).'); sz=dims(2:3);
        shapes={[dims(1) dims(1)],[dims(1) dims(3)],[dims(2) dims(1)],sz};
        vars={};
        for k=1:4
            v=f{k}; if isa(v,'symfun'), v=formula(v); end
            if isa(v,'sym'), vars=union(vars,arrayfun(@char,symvar(v),'UniformOutput',false)); end
            if ~isa(v,'function_handle')
                checkmatrix(v);
                if ~isequal(size(v),shapes{k}), error('ContourRoots:MatrixShape','Factor dimensions contradict Dimensions.'); end
            end
            f{k}=convert(v);
        end
        if numel(vars)>1, error('ContourRoots:MatrixRepresentation','All factors must use the same single symbolic variable.'); end
    end
    sz=double(sz(:).');
    for key={'Feedthrough','Delay','InitialValue'}
        v=o.(key{1});
        if ~isempty(v)
            validateattributes(v,{'numeric'},{'real','finite','size',sz});
            if strcmp(key{1},'Delay')&&any(v(:)<0), error('ContourRoots:MatrixRepresentation','Delays must be nonnegative.'); end
        end
    end
    for key={'InputNames','InputUnits','OutputNames','OutputUnits'}
        v=o.(key{1}); if isempty(v), continue; end
        n=sz(1); if startsWith(key{1},'Input'), n=sz(2); end
        if ~(iscellstr(v)||isstring(v))||numel(v)~=n, error('ContourRoots:MatrixShape','Channel labels/units must match the declared dimensions.'); end
        o.(key{1})=cellstr(v);
    end
    o.Size=sz; o.Dimensions=dims;
end
function d=merge(d,a)
    ix=isfinite(a); bad=ix & isfinite(d) & d~=a;
    if any(bad), error('ContourRoots:MatrixShape','Inconsistent factor dimensions.'); end
    d(ix)=a(ix);
end
function checkdims(a,n)
    if ~isnumeric(a)||~isreal(a)||~isvector(a)||numel(a)~=n||any(~isfinite(a)|a<1|a~=fix(a))
        error('ContourRoots:MatrixShape','Declare positive integer dimensions; handles are never probed.');
    end
end
function checkmatrix(v)
    if ~(isnumeric(v)||isa(v,'sym'))||isempty(v)||~ismatrix(v)
        error('ContourRoots:MatrixRepresentation','Factors must be nonempty numeric/symbolic matrices or handles.');
    end
    if isnumeric(v)&&any(~isfinite(v(:))), error('ContourRoots:MatrixRepresentation','Constants must be finite.'); end
end
function f=convert(f)
    if isa(f,'sym')
        v=symvar(f); if numel(v)>1, error('ContourRoots:MatrixRepresentation','Substitute parameters before construction.'); end
        if isempty(v), f=double(f); else, f=matlabFunction(f,'Vars',v); end
    elseif isnumeric(f)
        f=double(f);
    end
end
function checkchannel(v)
    if isnumeric(v)
        if ~isscalar(v)||~isfinite(v), error('ContourRoots:MatrixRepresentation','Numeric channels are scalar gains; use NDPAIR for coefficient vectors.'); end
    elseif isa(v,'sym')||isa(v,'symfun')
        if isa(v,'symfun'), v=formula(v); end
        if ~isscalar(v)||numel(symvar(v))>1, error('ContourRoots:MatrixRepresentation','Use scalar symbolic channels with parameters substituted.'); end
    elseif isa(v,'tf')||isa(v,'ss')||isa(v,'zpk')
        if ~isequal(size(v),[1 1])||v.Ts~=0, error('ContourRoots:MatrixRepresentation','Each LTI channel must be continuous SISO.'); end
        descriptor_check(v);
    elseif ~(isa(v,'function_handle')||(isstruct(v)&&isscalar(v)&&all(isfield(v,{'Numerator','Denominator'}))))
        error('ContourRoots:MatrixRepresentation','Unsupported scalar channel.');
    end
end
function descriptor_check(v)
    if isa(v,'ss')
        [~,~,~,~,E]=dssdata(v);
        if ~isempty(E)&&rank(full(E))<size(E,1)
            error('ContourRoots:MatrixRepresentation','Singular descriptor models require a separate support contract.');
        end
    end
end
function names=symbol_names(v)
    names={};
    if isa(v,'symfun'), v=formula(v); end
    if isa(v,'sym'), names=arrayfun(@char,symvar(v),'UniformOutput',false);
    elseif isstruct(v)&&all(isfield(v,{'Numerator','Denominator'}))
        names=union(symbol_names(v.Numerator),symbol_names(v.Denominator));
    end
end
