function [roots,info]=matrix_modes_search(H,region,varargin)
% Separate, dense, bounded matrix backend. No scalar determinant callbacks.
    o=struct('AssumeAnalytic',false,'Derivative',[],'DomainCheck',[], ...
        'Scaling','fixed','RootTolerance',1e-8,'RankTolerance',1e-8, ...
        'SeedPoints',[],'ContourPoints',12,'ContourRefinements',7, ...
        'MaxDepth',16,'MaxCells',2048,'MaxIterations',60, ...
        'MaxEvaluations',100000,'MaxMemoryMB',256, ...
        'Display',false,'Plot',false,'Warn',true);
    if mod(numel(varargin),2), error('ContourRoots:MatrixOption','Use name/value pairs.'); end
    names=fieldnames(o);
    for k=1:2:numel(varargin)
        ix=find(strcmpi(varargin{k},names),1);
        if isempty(ix), error('ContourRoots:MatrixOption','Unknown CMODES option.'); end
        o.(names{ix})=varargin{k+1};
    end
    validateattributes(region,{'numeric'},{'real','finite','vector','numel',4});
    region=double(region(:).');
    if region(1)>=region(2)||region(3)>=region(4)
        error('ContourRoots:MatrixDomain','Region must have positive width and height.');
    end
    for key={'AssumeAnalytic','Display','Plot','Warn'}
        validateattributes(o.(key{1}),{'logical','numeric'},{'scalar','binary'});
    end
    if ~o.AssumeAnalytic
        error('ContourRoots:MatrixAnalytic','Declare AssumeAnalytic=true for a regular analytic H on the closed region.');
    end
    for key={'RootTolerance','RankTolerance','MaxMemoryMB'}
        validateattributes(o.(key{1}),{'numeric'},{'scalar','real','finite','positive'});
    end
    for key={'ContourPoints','ContourRefinements','MaxCells','MaxIterations','MaxEvaluations'}
        validateattributes(o.(key{1}),{'numeric'},{'scalar','real','finite','integer','positive'});
    end
    validateattributes(o.MaxDepth,{'numeric'},{'scalar','real','finite','integer','nonnegative'});
    if ~isempty(o.SeedPoints), validateattributes(o.SeedPoints,{'numeric'},{'finite','vector'}); end
    if ~(ischar(o.Scaling)&&isrow(o.Scaling) || isstring(o.Scaling)&&isscalar(o.Scaling)) || ...
            ~any(strcmp(o.Scaling,{'fixed','none'}))
        error('ContourRoots:MatrixOption','Scaling must be fixed or none.');
    end
    for key={'Derivative','DomainCheck'}
        if ~isempty(o.(key{1}))&&~isa(o.(key{1}),'function_handle')
            error('ContourRoots:MatrixOption','Derivative and DomainCheck must be handles.');
        end
    end
    modelGuard=[]; n=[];
    if isa(H,'ContourRootsModel')
        if ~strcmp(H.Representation,'structured')
            error('ContourRoots:MatrixRepresentation','CMODES needs H or a CDYN model, not CMIMO.');
        end
        n=H.Dimensions(1); modelGuard=H.Options.DomainCheck; H=H.Factors{1};
    end
    if ~(isa(H,'function_handle')||isnumeric(H))
        error('ContourRoots:MatrixRepresentation','Use a matrix handle, numeric matrix or CDYN model.');
    end
    roots=complex(zeros(0,1)); localCounts=zeros(0,1); mult=zeros(0,1);
    nullity=zeros(0,1); radii=zeros(0,1); residuals=zeros(0,1); leftResiduals=zeros(0,1);
    rawResiduals=zeros(0,1); rawLeft=zeros(0,1); slopes=zeros(0,1);
    right={}; left={}; singularValues={}; unresolved=zeros(0,4); clusters=zeros(0,4);
    evaluations=0; factorizations=0; linearSolves=0; cells=0; svds=0;
    L=[]; R=[]; referenceNorm=NaN; referenceNode=NaN; stopReason='';
    history=struct('box',{},'count',{},'ok',{},'samples',{},'minPivotRatio',{});
    outer=emptycount();
    traceInfo=struct('available',~isempty(o.Derivative),'value',NaN,'ok',false, ...
        'reason','No analytic derivative supplied; trace cross-check omitted.');
    try
        % Pick an admissible, nonsingular reference. Transformations are then
        % frozen, preserving analyticity and the zeros throughout the search.
        refs=[region(1)+1i*region(3), region(2)+1i*region(4), ...
            mean(region(1:2))+1i*mean(region(3:4))];
        for referenceCandidate=refs
            A=raw(referenceCandidate,H); l=ones(n,1); r=ones(n,1);
            if strcmp(o.Scaling,'fixed')
                row=max(abs(A),[],2); row(row==0)=1;
                l=1./max(row,realmin); B=l.*A;
                col=max(abs(B),[],1).'; col(col==0)=1;
                r=1./max(col,realmin);
            end
            B=(l.*A).*r.';
            if all(isfinite(B(:))) && rcond(B)>eps
                L=l; R=r; referenceNorm=norm(B,'fro'); referenceNode=referenceCandidate; break
            end
        end
        if isempty(L)
            stopReason='No numerically regular scaling reference; regularity is unresolved.';
            unresolved=region;
        else
            outer=countbox(region,true);
            if outer.ok && outer.count>=0
                visit(region,outer,0);
            else
                unresolved=region; stopReason='Outer contour unresolved, singular, or inconsistent with analytic H.';
            end
        end
    catch e
        if ~strcmp(e.identifier,'ContourRoots:MatrixBudget'), rethrow(e); end
        unresolved=region; stopReason=e.message;
    end
    [~,ix]=sortrows([-real(roots) imag(roots)],[1 2]); roots=roots(ix);
    countComplete=outer.ok&&outer.count>=0;
    locationComplete=countComplete&&isempty(unresolved)&&isempty(clusters)&&sum(localCounts)==outer.count;
    multiplicityComplete=locationComplete&&all(isfinite(mult));
    complete=locationComplete&&multiplicityComplete;
    status='unresolved'; if complete, status='numerically_complete'; stopReason='All contour counts reconciled.'; end
    if isempty(stopReason), stopReason='Unresolved cells or local clusters remain.'; end
    info=struct('kind','matrix_characteristic','status',status,'complete',complete, ...
        'certified',false,'locationComplete',locationComplete,'multiplicityComplete',multiplicityComplete, ...
        'countComplete',countComplete,'count',outer.count,'region',region, ...
        'localCounts',localCounts(ix),'multiplicity',mult(ix),'nullity',nullity(ix), ...
        'locationRadius',radii(ix),'unresolvedBoxes',unresolved,'clusters',clusters, ...
        'residuals',residuals(ix),'leftResiduals',leftResiduals(ix), ...
        'unscaledResiduals',rawResiduals(ix),'unscaledLeftResiduals',rawLeft(ix), ...
        'simpleModeSlope',slopes(ix),'rightVectors',{right(ix)},'leftVectors',{left(ix)}, ...
        'singularValues',{singularValues(ix)},'evaluations',evaluations, ...
        'factorizations',factorizations,'linearSolves',linearSolves,'svds',svds, ...
        'cells',cells,'history',history,'traceCheck',traceInfo,'stopReason',stopReason, ...
        'scaling',struct('method',o.Scaling,'left',L,'right',R, ...
          'referenceNode',referenceNode,'referenceNorm',referenceNorm), ...
        'analyticSource','User assertion for H only','options',o);
    info.assumptions={'H is regular and analytic on a neighborhood of the closed rectangle.', ...
        'Numerical contour resolution is conditional, not a proof against undersampled oscillations.'};
    info.residualNormalization='Singular residuals of L*H*R divided by its fixed regular-reference Frobenius norm; original-coordinate absolute residuals also reported.';
    info.warnings={};
    if ~complete, info.warnings={stopReason}; end
    if ~isempty(clusters)
        info.warnings{end+1}='Local count, nullity, coincidence or projected derivative do not resolve clustered/defective structure; multiplicity is NaN.';
    end

    function A=raw(z,f)
        if ~inside(z,region), error('ContourRoots:MatrixDomain','Evaluation would leave the search rectangle.'); end
        guards={modelGuard,o.DomainCheck};
        for ig=1:2
            if isempty(guards{ig}), continue; end
            yes=guards{ig}(z);
            if ~(islogical(yes)&&isscalar(yes)&&yes)
                error('ContourRoots:MatrixDomain','DomainCheck rejected a required node.');
            end
        end
        if evaluations>=o.MaxEvaluations, error('ContourRoots:MatrixBudget','MaxEvaluations exhausted.'); end
        evaluations=evaluations+1;
        if isa(f,'function_handle'), A=f(z); else, A=f; end
        if ~isnumeric(A)||isempty(A)||~ismatrix(A)||size(A,1)~=size(A,2)||any(~isfinite(A(:)))
            error('ContourRoots:MatrixShape','H and Derivative must return finite square numeric matrices.');
        end
        if isempty(n), n=size(A,1); end
        if ~isequal(size(A),[n n]), error('ContourRoots:MatrixShape','Matrix size changed during the search.'); end
        % Conservative bound for dense LU/SVD/bordered workspaces. Callback
        % allocations and vendor BLAS workspaces are not controlled here.
        if 40*16*n*n>o.MaxMemoryMB*1024^2
            error('ContourRoots:MatrixBudget','Dense matrix work exceeds MaxMemoryMB.');
        end
        A=full(double(A));
    end
    function A=value(z), A=(L.*raw(z,H)).*R.'; end
    function D=derivative(z,box)
        if ~isempty(o.Derivative), D=(L.*raw(z,o.Derivative)).*R.'; return; end
        margin=min([real(z)-box(1),box(2)-real(z),imag(z)-box(3),box(4)-imag(z)]);
        h=min(1e-4*max(1,abs(z)),margin/3);
        if h<=16*eps*max(1,abs(z)), D=nan(n); return; end
        D=(value(z+h)-value(z-h)-1i*(value(z+1i*h)-value(z-1i*h)))/(4*h);
    end
    function c=countbox(box,checkTrace)
        c=emptycount(); previous=NaN; stable=0;
        for level=0:o.ContourRefinements
            np=o.ContourPoints*2^level;
            if np*4*128>o.MaxMemoryMB*1024^2, error('ContourRoots:MatrixBudget','Contour workspace exceeds MaxMemoryMB.'); end
            t=(0:np-1)/np; x0=box(1); x1=box(2); y0=box(3); y1=box(4);
            z=[x0+(x1-x0)*t+1i*y0, x1+1i*(y0+(y1-y0)*t), ...
                x1-(x1-x0)*t+1i*y1, x0+1i*(y1-(y1-y0)*t)];
            phase=zeros(size(z)); logs=phase; tr=phase; minRatio=Inf;
            for j=1:numel(z)
                A=value(z(j)); [a,b,p]=lu(A); factorizations=factorizations+1;
                d=diag(b);
                if any(d==0)||any(~isfinite(d)), record(); return; end
                % Permutation determinant is exactly +/-1; count inversions
                % rather than taking the determinant of the permutation.
                [~,perm]=max(p,[],2); parity=0;
                for ip=1:n, parity=parity+sum(perm(ip+1:end)<perm(ip)); end
                phase(j)=sum(angle(d))+pi*mod(parity,2);
                logs(j)=sum(log(abs(d)));
                minRatio=min(minRatio,min(abs(d))/max(abs(d)));
                if checkTrace && ~isempty(o.Derivative)
                    D=derivative(z(j),region);
                    tr(j)=trace(b\(a\(p*D))); linearSolves=linearSolves+1;
                end
            end
            nxt=[2:numel(z) 1]; dp=angle(exp(1i*(phase(nxt)-phase)));
            dm=logs(nxt)-logs; winding=sum(dp)/(2*pi); total=round(winding);
            resolved=max(abs(dp))<pi/3&&max(abs(dm))<2&&abs(winding-total)<1e-7;
            if checkTrace&&~isempty(o.Derivative)
                tv=sum((tr+tr(nxt)).*(z(nxt)-z)/2)/(2i*pi);
                traceInfo.value=tv; traceInfo.ok=isfinite(tv)&&abs(tv-total)<1e-4*max(1,abs(total));
                traceInfo.reason='Trapezoidal trace integral compared with LU winding on the outer contour.';
                resolved=resolved&&traceInfo.ok;
            end
            if resolved&&total==previous, stable=stable+1; else, stable=0; end
            previous=total; c.samples=numel(z); c.minPivotRatio=minRatio;
            if stable>=2
                c.ok=true; c.count=total;
                if total~=0, c.center=sum((z+z(nxt))/2.*(dm+1i*dp))/(2i*pi*total); end
                record(); return
            end
        end
        record();
        function record()
            history(end+1)=struct('box',box,'count',c.count,'ok',c.ok, ...
                'samples',c.samples,'minPivotRatio',c.minPivotRatio);
        end
    end
    function visit(box,c,depth)
        if cells>=o.MaxCells, error('ContourRoots:MatrixBudget','MaxCells exhausted.'); end
        cells=cells+1;
        if c.count==0, return; end
        seeds=[c.center,mean(box(1:2))+1i*mean(box(3:4)),o.SeedPoints(:).'];
        for z0=seeds
            if ~isfinite(z0)||~strictinside(z0,box), continue; end
            [z,ok]=refine(z0,box,c.count);
            if ok && accept(z,c.count,box), return; end
        end
        if depth>=o.MaxDepth, unresolved(end+1,:)=box; return; end
        axes=[1 2]; if diff(box(3:4))>diff(box(1:2)), axes=[2 1]; end
        for ax=axes
            for fraction=[.5 .4384471872 .5732050808]
                q=2*ax-1; mid=box(q)+fraction*(box(q+1)-box(q)); a=box; b=box;
                a(q+1)=mid; b(q)=mid; ca=countbox(a,false); cb=countbox(b,false);
                if ca.ok&&cb.ok&&ca.count>=0&&cb.count>=0&&ca.count+cb.count==c.count
                    visit(a,ca,depth+1); visit(b,cb,depth+1); return
                end
            end
        end
        unresolved(end+1,:)=box;
    end
    function [z,ok]=refine(z,box,m)
        ok=false;
        for it=1:o.MaxIterations
            A=value(z); [U,S,V]=svd(A); svds=svds+1; sigma=S(end,end);
            D=derivative(z,box); if any(~isfinite(D(:))), return; end
            v=V(:,end); w=U(:,end); slope=w'*D*v;
            if sigma==0, ok=true; return; end
            if m==1
                % Bordered Newton, with a fixed normalization for this step.
                K=[A D*v;v' 0]; rhs=[-A*v;0];
                if rcond(K)<eps, return; end
                delta=K\rhs; linearSolves=linearSolves+1; factorizations=factorizations+1;
                step=-delta(end);
            else
                if rcond(A)>eps
                    tr=trace(A\D); linearSolves=linearSolves+1; factorizations=factorizations+1;
                    step=m/tr;
                elseif abs(slope)>eps*norm(D,'fro')
                    step=sigma/slope;
                else
                    % A nearly singular defective candidate still needs a
                    % local count; it is NOT classified as a simple mode.
                    ok=sigma/referenceNorm<1e-12; return
                end
            end
            if ~isfinite(step), return; end
            cap=max([diff(box(1:2)),diff(box(3:4))])/2;
            if abs(step)>cap, step=step*cap/abs(step); end
            accepted=false;
            for bt=0:16
                trial=z-step/2^bt;
                if ~strictinside(trial,box), continue; end
                At=value(trial); st=svd(At); svds=svds+1;
                if st(end)<=sigma*(1+1e-6), accepted=true; break; end
            end
            if ~accepted, return; end
            z=trial;
            if abs(step/2^bt)<min(1e-12,o.RootTolerance/100)*(1+abs(z)) && st(end)/referenceNorm<1e-10
                ok=true; return
            end
        end
    end
    function yes=accept(z,expected,parent)
        yes=false; radius=o.RootTolerance*(1+abs(z));
        if any(abs(roots-z)<=2*radius), return; end
        box=[real(z)-radius real(z)+radius imag(z)-radius imag(z)+radius];
        % Containment is essential to conservation: a local count cannot
        % borrow roots from a neighboring cell to satisfy this parent's count.
        if box(1)<=parent(1)||box(2)>=parent(2)||box(3)<=parent(3)||box(4)>=parent(4), return; end
        c=countbox(box,false);
        if ~c.ok||c.count~=expected||c.count<1, return; end
        A=value(z); [U,S,V]=svd(A); svds=svds+1; sig=diag(S);
        g=sum(sig<=o.RankTolerance*referenceNorm);
        if g==0||sig(end)/referenceNorm>1e-10, return; end
        W=U(:,end-g+1:end); V=V(:,end-g+1:end);
        vr=R.*V; wl=conj(L).*W;
        vr=vr./vecnorm(vr); wl=wl./vecnorm(wl);
        original=raw(z,H); D=derivative(z,region);
        smin=NaN;
        if expected==1 && g==1, smin=abs(W'*D*V)/referenceNorm; end
        multiplicity=NaN;
        % A nonsingular derivative on the nullspaces supports semisimple
        % multiplicity. Higher chains and unresolved clusters stay explicit.
        coupling=svd(W'*D*V); svds=svds+1;
        coincident=expected==1 || norm(A*V,'fro')<=100*eps*referenceNorm;
        if expected==g && coincident && all(isfinite(coupling)) && min(coupling)>o.RankTolerance*max(norm(D,'fro'),realmin)
            multiplicity=expected;
        else
            clusters(end+1,:)=box;
        end
        roots(end+1,1)=z; localCounts(end+1,1)=expected; mult(end+1,1)=multiplicity;
        nullity(end+1,1)=g; radii(end+1,1)=sqrt(2)*radius;
        residuals(end+1,1)=norm(A*V,'fro')/referenceNorm;
        leftResiduals(end+1,1)=norm(W'*A,'fro')/referenceNorm;
        rawResiduals(end+1,1)=norm(original*vr,'fro'); rawLeft(end+1,1)=norm(wl'*original,'fro');
        slopes(end+1,1)=smin; right{end+1,1}=vr; left{end+1,1}=wl; singularValues{end+1,1}=sig;
        yes=true;
    end
end
function c=emptycount()
    c=struct('ok',false,'count',NaN,'center',NaN,'samples',0,'minPivotRatio',NaN);
end
function yes=inside(z,b)
    yes=real(z)>=b(1)&&real(z)<=b(2)&&imag(z)>=b(3)&&imag(z)<=b(4);
end
function yes=strictinside(z,b)
    yes=real(z)>b(1)&&real(z)<b(2)&&imag(z)>b(3)&&imag(z)<b(4);
end
