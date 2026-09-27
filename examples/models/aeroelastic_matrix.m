function D = aeroelastic_matrix(s,U,model,method)
%AEROELASTIC_MATRIX Dynamic stiffness of a 2-DOF typical section.
%   D=s^2*M+s*Cs+K-A(s*b/U); s must be a scalar.
%   Methods: exact (default), jones, quasisteady, pk. The pk function is
%   NONANALYTIC in s and must never be passed to contour root counting.
%   A is THEODORSEN_LOADS (Kaiser & Quero 2022, Eq. (33), with e=a).
%   At U=0 the still-fluid added mass remains; rho=0 is the dry structure.
    if nargin<4, method='exact'; end
    validateattributes(s,{'numeric'},{'scalar','finite'});
    validateattributes(U,{'numeric'},{'scalar','real','finite','nonnegative'});
    A=theodorsen_loads(s,U,model.b,model.a,model.rho,method);
    D=s^2*model.M+s*model.Cs+model.K-A;
end
