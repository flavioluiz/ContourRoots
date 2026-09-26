function [r,residuals] = aeroelastic_pk_roots(U,model,seeds)
%AEROELASTIC_PK_ROOTS Real 2-variable Newton solve of the p-k equations.
%   The entire aerodynamic matrix is evaluated at i*imag(s)*b/U, as in
%   Kaiser & Quero (2022), Eq. (4). No contour completeness claim is made.
    r=complex(NaN(size(seeds))); residuals=inf(size(seeds));
    for k=1:numel(seeds)
        x=[real(seeds(k));imag(seeds(k))];
        for iter=1:60
            f=equations(x);
            if norm(f)<1e-11, break; end
            h=1e-5*max(1,norm(x)); J=zeros(2);
            for j=1:2
                d=zeros(2,1); d(j)=h;
                J(:,j)=(equations(x+d)-equations(x-d))/(2*h);
            end
            if rcond(J)<1e-12, break; end
            step=J\f; improved=false;
            for bt=0:12
                y=x-step/2^bt;
                if y(2)>0 && norm(equations(y))<norm(f)
                    x=y; improved=true; break
                end
            end
            if ~improved, break; end
        end
        residuals(k)=norm(equations(x));
        if residuals(k)<1e-8, r(k)=x(1)+1i*x(2); end
    end
    if any(~isfinite(r)), error('aeroelastic:PK','A p-k iteration failed.'); end
    if numel(r)>1 && abs(r(1)-r(2))<1e-5
        error('aeroelastic:PK','Two starts converged to the same p-k root.');
    end
    function f=equations(x)
        v=det(aeroelastic_matrix(x(1)+1i*x(2),U,model,'pk'))/det(model.K);
        f=[real(v);imag(v)];
    end
end
