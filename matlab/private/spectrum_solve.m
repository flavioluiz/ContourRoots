function [r,m,info] = spectrum_solve(f,df,box,opt,exploratory)
% Isolate analytic zeros using conserved contour counts and subdivision.
    r=complex(zeros(0,1)); m=zeros(0,1); radii=zeros(0,1);
    unresolved=zeros(0,4); cells=0;
    outer=spectrum_count(f,box,opt);
    if exploratory
        % A meromorphic black box may have Z-P=0 despite containing poles.
        % Newton candidates are locally checked, but completeness is unknown.
        [x,y]=ndgrid(linspace(box(1),box(2),opt.GridSize(1)), ...
                     linspace(box(3),box(4),opt.GridSize(2)));
        starts=[opt.SeedPoints(:); x(:)+1i*y(:)];
        for j=1:numel(starts)
            [z,ok]=spectrum_newton(f,df,starts(j),box,1,opt);
            if ok, accept_local(z,Inf); end
        end
    elseif outer.ok && outer.count>=0
        visit(box,outer,0);
    else
        unresolved=box;
    end
    [~,idx]=sortrows([-real(r) imag(r)],[1 2]);
    r=r(idx); m=m(idx); radii=radii(idx);
    complete=~exploratory && outer.ok && isempty(unresolved) && sum(m)==outer.count;
    info=struct('complete',complete,'count',outer.count,'contourResolved',outer.ok, ...
        'unresolvedBoxes',unresolved,'cellsVisited',cells,'locationRadius',radii, ...
        'residuals',abs(spectrum_values(f,r)),'exploratory',exploratory);

    function visit(b,c,depth)
        cells=cells+1;
        if cells>opt.MaxCells
            unresolved(end+1,:)=b; return
        end
        if c.count==0, return; end
        center=mean(b(1:2))+1i*mean(b(3:4));
        starts=[c.center; center; opt.SeedPoints(:)];
        starts=starts(isfinite(starts) & real(starts)>b(1) & real(starts)<b(2) ...
            & imag(starts)>b(3) & imag(starts)<b(4));
        for jj=1:numel(starts)
            [z,ok]=spectrum_newton(f,df,starts(jj),b,c.count,opt);
            if ok && accept_local(z,c.count), return; end
        end
        if depth>=opt.MaxDepth
            unresolved(end+1,:)=b; return
        end
        % Shift split lines when a zero lies on a proposed internal boundary.
        fractions=[.5 .4384471872 .5732050808];
        for axis=1:2
            if diff(b(1:2))>=diff(b(3:4)), dim=axis; else, dim=3-axis; end
            for frac=fractions
                if dim==1
                    mid=b(1)+frac*diff(b(1:2));
                    left=[b(1) mid b(3:4)]; right=[mid b(2:4)];
                else
                    mid=b(3)+frac*diff(b(3:4));
                    left=[b(1:3) mid]; right=[b(1:2) mid b(4)];
                end
                a=spectrum_count(f,left,opt); d=spectrum_count(f,right,opt);
                if a.ok && d.ok && a.count>=0 && d.count>=0 && a.count+d.count==c.count
                    visit(left,a,depth+1); visit(right,d,depth+1); return
                end
            end
        end
        unresolved(end+1,:)=b;
    end

    function accepted=accept_local(z,expected)
        accepted=false;
        radius=opt.RootTolerance*(1+abs(z));
        if any(abs(r-z)<4*radius), return; end
        local=[real(z)-radius real(z)+radius imag(z)-radius imag(z)+radius];
        if local(1)<=box(1) || local(2)>=box(2) || local(3)<=box(3) || local(4)>=box(4)
            return
        end
        count=spectrum_count(f,local,opt);
        if count.ok && count.count>0 && (isinf(expected) || count.count==expected)
            r(end+1,1)=z; m(end+1,1)=count.count; radii(end+1,1)=sqrt(2)*radius;
            accepted=true;
        end
    end
end
