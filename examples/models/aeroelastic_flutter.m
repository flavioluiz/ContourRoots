function result = aeroelastic_flutter(model,method)
%AEROELASTIC_FLUTTER Local flutter candidates from Re/Im det D(i*w,U)=0.
%   Base MATLAB only. Several independent starts, no published flutter
%   value used as a seed. This is NOT a proof that every crossing was found.
%   Confirm onset separately with a contour sweep and static divergence.
    if nargin<2, method='exact'; end
    wh=sqrt(model.K(1,1)/model.M(1,1));
    wa=sqrt(model.K(2,2)/model.M(2,2)); ub=model.b*wa;
    found=zeros(0,3);
    for wf=[.65 1 1.35]
        for uf=[1 2 4]
            x=log([wf*sqrt(wh*wa);uf*ub]);
            for iter=1:60
                f=equations(x);
                if norm(f)<1e-12, break; end
                h=1e-5; J=zeros(2);
                for j=1:2
                    d=zeros(2,1); d(j)=h;
                    J(:,j)=(equations(x+d)-equations(x-d))/(2*h);
                end
                if rcond(J)<1e-12, break; end
                step=J\f; step=step/max(1,norm(step));
                improved=false;
                for bt=0:14
                    y=x-step/2^bt;
                    if exp(y(2))<.05*ub || exp(y(2))>20*ub || ...
                            exp(y(1))<.05*wh || exp(y(1))>5*wa, continue; end
                    if norm(equations(y))<norm(f), improved=true; x=y; break; end
                end
                if ~improved, break; end
            end
            if norm(equations(x))<1e-9 && exp(x(2))>.05*ub
                w=exp(x(1)); U=exp(x(2));
                if isempty(found) || all(abs(found(:,2)-U)>1e-5*(1+U))
                    found(end+1,:)=[w U norm(equations(x))]; %#ok<AGROW>
                end
            end
        end
    end
    if isempty(found), error('aeroelastic:Flutter','No converged flutter candidate.'); end
    found=sortrows(found,2); w=found(1,1); U=found(1,2); s=1i*w;
    h=1e-4; hu=1e-4*U;
    if strcmpi(method,'pk')
        speed=NaN; % Not a holomorphic characteristic function.
    else
        ds=(aeroelastic_delta(s+h,U,model,method)-aeroelastic_delta(s-h,U,model,method))/(2*h);
        du=(aeroelastic_delta(s,U+hu,model,method)-aeroelastic_delta(s,U-hu,model,method))/(2*hu);
        speed=real(-du/ds);
    end
    result=struct('U',U,'omega',w,'k',w*model.b/U,'residual',found(1,3), ...
        'crossingSpeed',speed,'candidates',found,'method',method);
    function f=equations(x)
        dynamic=aeroelastic_matrix(1i*exp(x(1)),exp(x(2)),model,method);
        v=det(dynamic)/det(model.K); f=[real(v);imag(v)];
    end
end
