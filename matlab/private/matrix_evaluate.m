function [G,info]=matrix_evaluate(M,s,varargin)
    if ~isa(M,'ContourRootsModel')||~isscalar(M), error('ContourRoots:MatrixRepresentation','Use one CMIMO/CDYN model.'); end
    validateattributes(s,{'numeric'},{'vector','nonempty','finite'}); s=s(:); ns=numel(s);
    o=struct('MaxEvaluations',1e6,'MaxMemoryMB',256,'ApplyDelay',true);
    if mod(numel(varargin),2), error('ContourRoots:MatrixOption','Use name/value pairs.'); end
    for k=1:2:numel(varargin)
        names=fieldnames(o); ix=find(strcmpi(varargin{k},names),1);
        if isempty(ix), error('ContourRoots:MatrixOption','Unknown evaluation option.'); end
        o.(names{ix})=varargin{k+1};
    end
    validateattributes(o.MaxEvaluations,{'numeric'},{'scalar','positive','integer','finite'});
    validateattributes(o.MaxMemoryMB,{'numeric'},{'scalar','positive','finite'});
    if ns>o.MaxEvaluations, error('ContourRoots:MatrixBudget','Node budget exceeded.'); end
    sz=M.Size; work=prod(sz);
    if strcmp(M.Representation,'structured'), n=M.Dimensions(1); work=work+3*n*n+n*sum(sz); end
    bytes=16*prod(sz)*ns; budget=o.MaxMemoryMB*2^20;
    if bytes+32*work>budget, error('ContourRoots:MatrixBudget','Output and one matrix solve exceed MaxMemoryMB.'); end
    chunk=max(1,floor((budget-bytes)/(32*work)));
    G=complex(zeros(sz(1),sz(2),ns));
    info=struct('evaluations',ns,'factorEvaluations',0,'factorizations',0,'linearSolves',0,'rhsColumns',0,'certified',false);
    if ~isempty(M.Options.DomainCheck)
        for k=1:ns
            yes=M.Options.DomainCheck(s(k));
            if ~islogical(yes)||~isscalar(yes)||~yes, error('ContourRoots:MatrixDomain','A requested node violates DomainCheck.'); end
        end
    end
    for start=1:chunk:ns
        ix=start:min(ns,start+chunk-1); z=s(ix);
        switch M.Representation
            case 'structured'
                n=M.Dimensions(1); shapes={[n n],[n sz(2)],[sz(1) n],sz}; v=cell(1,4);
                for a=1:4
                    [v{a},calls]=factor(M.Factors{a},z,shapes{a},M.Options.Batched);
                    info.factorEvaluations=info.factorEvaluations+calls;
                end
                for a=1:numel(ix)
                    H=v{1}(:,:,a); X=H\v{2}(:,:,a);
                    % Reject unusable solves instead of converting them to zeros.
                    if any(~isfinite(X(:))), error('ContourRoots:MatrixEvaluation','Nonfinite structured solve at a required node.'); end
                    G(:,:,ix(a))=v{3}(:,:,a)*X+v{4}(:,:,a);
                end
                info.factorizations=info.factorizations+numel(ix); info.linearSolves=info.linearSolves+numel(ix);
                info.rhsColumns=info.rhsColumns+numel(ix)*sz(2);
            case 'channels'
                for j=1:sz(2), for i=1:sz(1)
                    G(i,j,ix)=reshape(channel(M.Factors{i,j},z),1,1,[]);
                end, end
                info.factorEvaluations=info.factorEvaluations+numel(ix)*prod(sz);
            otherwise
                [G(:,:,ix),calls]=factor(M.Factors,z,sz,M.Options.Batched);
                info.factorEvaluations=info.factorEvaluations+calls;
        end
    end
    if o.ApplyDelay && ~isempty(M.Options.Delay)
        G=G.*exp(-M.Options.Delay.*reshape(s,1,1,[]));
    end
    if any(~isfinite(G(:))), error('ContourRoots:MatrixEvaluation','Transfer values must be finite at required regular nodes.'); end
end
function [v,calls]=factor(f,s,sz,batched)
    calls=0;
    if isnumeric(f), v=repmat(full(f),1,1,numel(s)); return; end
    if batched
        v=f(s); calls=1; check(v,sz,numel(s));
    else
        v=complex(zeros(sz(1),sz(2),numel(s)));
        for k=1:numel(s), q=f(s(k)); calls=calls+1; check(q,sz,1); v(:,:,k)=q; end
    end
end
function check(v,sz,n)
    if ~isnumeric(v)||size(v,1)~=sz(1)||size(v,2)~=sz(2)||size(v,3)~=n||ndims(v)>3||any(~isfinite(v(:)))
        error('ContourRoots:MatrixShape','Evaluator must return finite rows-by-columns-by-nodes values in the declared layout.');
    end
end
function v=channel(f,s)
    if isnumeric(f), v=f+zeros(size(s));
    elseif isa(f,'function_handle'), v=response_values(f,s);
    elseif isa(f,'sym')||isa(f,'symfun')
        if isa(f,'symfun'), f=formula(f); end
        vars=symvar(f); if isempty(vars), v=double(f)+zeros(size(s)); else, v=response_values(matlabFunction(f,'Vars',vars),s); end
    elseif isstruct(f), v=channel_factor(f.Numerator,s)./channel_factor(f.Denominator,s);
    else, v=arrayfun(@(z) evalfr(f,z),s); end
end
function v=channel_factor(f,s)
    if isnumeric(f), v=polyval(f,s); else, v=channel(f,s); end
end
