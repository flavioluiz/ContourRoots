function A = theodorsen_loads(s,U,b,a,rho,method)
%THEODORSEN_LOADS Unsteady air loads on a thin airfoil strip, Laplace domain.
%   A = THEODORSEN_LOADS(S,U,B,A,RHO) returns the 2x2 matrix such that
%       [downward force per span; nose-up moment per span] = A*[h; alpha],
%   h = downward plunge at the elastic axis, alpha = nose-up pitch, for a
%   strip of semichord B with the elastic axis A*B aft of midchord, in
%   incompressible flow of density RHO at speed U. S is a scalar.
%   Noncirculatory (added mass) and circulatory terms follow Kaiser & Quero
%   (2022), Eq. (33), with Theodorsen's function C(p), p = S*B/U, evaluated
%   EXACTLY from Bessel functions (THEODORSEN_LAPLACE). Analytic in S except
%   on the negative real axis (branch cut).
%   U = 0 keeps the still-fluid added mass; RHO = 0 gives A = 0 (vacuum).
%
%   A = THEODORSEN_LOADS(...,METHOD) with METHOD = 'exact' (default),
%   'jones' (rational approximation), 'quasisteady' (C = 1) or 'pk'
%   (C evaluated at i*Im(p): NOT analytic, never use it in CROOTS).
    if nargin<6, method='exact'; end
    validateattributes(s,{'numeric'},{'scalar','finite'});
    validateattributes(U,{'numeric'},{'scalar','real','finite','nonnegative'});
    A2=[-1 a*b;a*b -(1/8+a^2)*b^2];
    if rho==0, A=zeros(2); return; end
    if U==0, A=pi*rho*b^2*s^2*A2; return; end
    p=s*b/U;
    switch lower(method)
        case 'exact', C=theodorsen_laplace(p);
        case 'jones', C=1-.165*p/(p+.0455)-.335*p/(p+.3);
        case 'quasisteady', C=1;
        case 'pk', p=1i*imag(p); C=theodorsen_laplace(p);
        otherwise, error('aeroelastic:Method','Unknown aerodynamic method.');
    end
    A1=[-2*C (-1-2*C*(.5-a))*b; ...
        2*C*(.5+a)*b (.5-a)*(2*C*(.5+a)-1)*b^2];
    A0=[0 -2*C*b;0 2*C*(.5+a)*b^2];
    A=pi*rho*U^2*(p^2*A2+p*A1+A0);
end
