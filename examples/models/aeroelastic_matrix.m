function D = aeroelastic_matrix(s,U,model,method)
%AEROELASTIC_MATRIX Dynamic stiffness of a 2-DOF typical section.
%   D=s^2*M+s*Cs+K-A(s*b/U); s must be a scalar.
%   Methods: exact (default), jones, quasisteady, pk. The pk function is
%   NONANALYTIC in s and must never be passed to contour root counting.
%   A follows Kaiser & Quero (2022), Eq. (33), with e=a.
%   At U=0 the still-fluid added mass remains; rho=0 is the dry structure.
    if nargin<4, method='exact'; end
    validateattributes(s,{'numeric'},{'scalar','finite'});
    validateattributes(U,{'numeric'},{'scalar','real','finite','nonnegative'});
    b=model.b; a=model.a; r=model.rho;
    A2=[-1 a*b;a*b -(1/8+a^2)*b^2];
    if U==0
        A=pi*r*b^2*s^2*A2;
    else
        p=s*b/U;
        switch lower(method)
            case 'exact', C=theodorsen_laplace(p);
            case 'jones', C=1-.165*p/(p+.0455)-.335*p/(p+.3);
            case 'quasisteady', C=1;
            case 'pk'
                p=1i*imag(p); C=theodorsen_laplace(p);
            otherwise, error('aeroelastic:Method','Unknown aerodynamic method.');
        end
        A1=[-2*C (-1-2*C*(.5-a))*b; ...
            2*C*(.5+a)*b (.5-a)*(2*C*(.5+a)-1)*b^2];
        A0=[0 -2*C*b;0 2*C*(.5+a)*b^2];
        A=pi*r*U^2*(p^2*A2+p*A1+A0);
    end
    D=s^2*model.M+s*model.Cs+model.K-A;
end
