function [models,stats]=matrix_response_context(M,o,extra)
% One bounded cache per operation, shared across channels and kernel orders.
% Never persistent, never copied into a prepared kernel snapshot.
    sz=M.Size; originals=cell(sz); models=cell(sz);
    for j=1:sz(2), for i=1:sz(1)
        co=o;
        for metadataKey={'Feedthrough','InitialValue'}
            name=metadataKey{1}; v=M.Options.(name); w=extra.(name);
            if ~isempty(v)&&~isempty(w)&&any(v(:)~=w(:))
                error('ContourRoots:ResponseModel','Response metadata contradicts the matrix model.');
            end
            if ~isempty(w), v=w; end
            if ~isempty(v), co.(name)=v(i,j); else, co.(name)=[]; end
        end
        if strcmp(M.Representation,'channels'), source=M.Factors{i,j};
        elseif strcmp(M.Representation,'transfer')&&isnumeric(M.Factors), source=M.Factors(i,j);
        else, source=@(s) NaN(size(s)); end % metadata only, never evaluated
        m=response_model(source,co);
        m.assumptions{end+1}=M.Options.DomainDescription;
        if ~isempty(M.Options.Delay), m.delay=m.delay+M.Options.Delay(i,j); end
        originals{i,j}=m;
        m.full=@(s) values(s,i,j); D=m.D;
        full=m.full; m.regular=@(s) response_values(full,s)-D;
        models{i,j}=m;
    end, end
    cache=containers.Map('KeyType','char','ValueType','any');
    capacity=max(1,floor(o.MaxMemoryMB*2^20/4/(16*prod(sz)+400)));
    count=struct('evaluations',0,'factorEvaluations',0,'factorizations',0,'linearSolves',0,'rhsColumns',0,'cacheHits',0);
    stats=@snapshot;
    function result=snapshot
        result=count;
    end
    function v=values(s,i,j)
        if strcmp(M.Representation,'channels')
            % Independent scalar channels: evaluate only the requested one.
            % Their adaptive grids rarely coincide, so filling whole pages
            % would multiply the work by ny*nu for little sharing.
            if count.evaluations+numel(s)>extra.MaxEvaluations
                error('ContourRoots:MatrixBudget','Shared matrix-node evaluation budget exhausted.');
            end
            if ~isempty(M.Options.DomainCheck)
                for k=1:numel(s)
                    ok=M.Options.DomainCheck(s(k));
                    if ~islogical(ok)||~isscalar(ok)||~ok, error('ContourRoots:MatrixDomain','Node violates DomainCheck.'); end
                end
            end
            v=response_values(originals{i,j}.full,s);
            count.evaluations=count.evaluations+numel(s);
            count.factorEvaluations=count.factorEvaluations+numel(s);
            return
        end
        v=zeros(size(s)); z=s(:);
        % Bounded chunks also apply when a scalar inversion requests a huge grid.
        for first=1:min(1024,capacity):numel(z)
            ix=first:min(numel(z),first+min(1024,capacity)-1); zz=z(ix);
            keys=arrayfun(@key,zz,'UniformOutput',false);
            hit=isKey(cache,keys); count.cacheHits=count.cacheHits+sum(hit);
            for k=find(hit(:)).', q=cache(keys{k}); v(ix(k))=q(i,j); end
            missing=find(~hit); if isempty(missing), continue; end
            [nodes,~,map]=unique(zz(missing));
            if count.evaluations+numel(nodes)>extra.MaxEvaluations
                error('ContourRoots:MatrixBudget','Shared matrix-node evaluation budget exhausted.');
            end
            if cache.Count+numel(nodes)>capacity, cache=containers.Map('KeyType','char','ValueType','any'); end
            if strcmp(M.Representation,'channels')
                if ~isempty(M.Options.DomainCheck)
                    for k=1:numel(nodes)
                        ok=M.Options.DomainCheck(nodes(k));
                        if ~islogical(ok)||~isscalar(ok)||~ok, error('ContourRoots:MatrixDomain','Node violates DomainCheck.'); end
                    end
                end
                pages=complex(zeros(sz(1),sz(2),numel(nodes)));
                for jj=1:sz(2), for ii=1:sz(1)
                    pages(ii,jj,:)=reshape(response_values(originals{ii,jj}.full,nodes),1,1,[]);
                end, end
                work=struct('evaluations',numel(nodes),'factorEvaluations',numel(nodes)*prod(sz), ...
                    'factorizations',0,'linearSolves',0,'rhsColumns',0);
            else
                [pages,work]=matrix_evaluate(M,nodes,'MaxMemoryMB',o.MaxMemoryMB/4,'ApplyDelay',false);
            end
            for name={'evaluations','factorEvaluations','factorizations','linearSolves','rhsColumns'}
                field=name{1}; count.(field)=count.(field)+work.(field);
            end
            for k=1:numel(nodes), cache(key(nodes(k)))=pages(:,:,k); end
            for k=1:numel(missing), v(ix(missing(k)))=pages(i,j,map(k)); end
        end
    end
end
function k=key(z)
    k=[num2hex(real(z)) num2hex(imag(z))];
end
