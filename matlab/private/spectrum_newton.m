function [z,ok] = spectrum_newton(f,df,z,box,m,opt)
% Damped Newton, optionally with multiplicity correction.
%   The iteration never evaluates f outside the closed cell BOX: a trial
%   step that would leave the cell is shortened, and the finite-difference
%   stencil is kept inside it. Functions that are analytic only on a
%   neighborhood of the search region (for example with a branch cut just
%   outside it) are therefore never evaluated where they are undefined.
    ok = false;
    width = max([diff(box(1:2)) diff(box(3:4))]);
    for k=1:opt.MaxIterations
        v = f(z);
        if ~isfinite(v), return; end
        if v==0, ok = true; return; end
        if isempty(df)
            % Four complex directions: O(h^4), valid for analytic functions.
            h = 1e-4*max(1,abs(z));
            dist = min([real(z)-box(1) box(2)-real(z) imag(z)-box(3) box(4)-imag(z)]);
            h = max(min(h,0.9*dist),1e-12*max(1,abs(z)));
            d = (f(z+h)-f(z-h)-1i*(f(z+1i*h)-f(z-1i*h)))/(4*h);
        else
            d = df(z);
        end
        if ~isfinite(d) || d==0, return; end
        step = m*v/d;
        if abs(step)>width/2, step=step*(width/2)/abs(step); end
        alpha = 1; vt = Inf;
        for backtrack=1:12
            trial=z-alpha*step;
            if inside(trial,box)
                vt=f(trial);
                if isfinite(vt) && abs(vt)<=abs(v), break; end
            end
            alpha=alpha/2;
        end
        if ~inside(trial,box) || ~isfinite(vt) || abs(vt)>abs(v), return; end
        z=trial;
        if abs(alpha*step)<=opt.RootTolerance*(1+abs(z))
            ok=true; return
        end
    end
end

function yes = inside(z,box)
    yes = real(z)>=box(1) && real(z)<=box(2) && imag(z)>=box(3) && imag(z)<=box(4);
end
