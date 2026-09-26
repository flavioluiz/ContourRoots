function [z,ok] = spectrum_newton(f,df,z,box,m,opt)
% Damped Newton, optionally with multiplicity correction.
    ok = false;
    width = max([diff(box(1:2)) diff(box(3:4))]);
    for k=1:opt.MaxIterations
        v = f(z);
        if ~isfinite(v), return; end
        if v==0, ok = true; return; end
        if isempty(df)
            h = 1e-4*max(1,abs(z));
            % Four complex directions: O(h^4), valid for analytic functions.
            d = (f(z+h)-f(z-h)-1i*(f(z+1i*h)-f(z-1i*h)))/(4*h);
        else
            d = df(z);
        end
        if ~isfinite(d) || d==0, return; end
        step = m*v/d;
        if abs(step)>width/2, step=step*(width/2)/abs(step); end
        alpha = 1;
        for backtrack=1:12
            trial=z-alpha*step; vt=f(trial);
            if isfinite(vt) && abs(vt)<=abs(v), break; end
            alpha=alpha/2;
        end
        if ~isfinite(vt) || abs(vt)>abs(v), return; end
        z=trial;
        if real(z)<box(1) || real(z)>box(2) || imag(z)<box(3) || imag(z)>box(4)
            return
        end
        if abs(alpha*step)<=opt.RootTolerance*(1+abs(z))
            ok=true; return
        end
    end
end
