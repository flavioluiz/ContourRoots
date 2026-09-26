function value = aeroelastic_delta(s,U,model,method)
%AEROELASTIC_DELTA Characteristic function of the 2-DOF aeroelastic section.
%   V = AEROELASTIC_DELTA(S,U,MODEL) returns det(D(s,U)) for every element
%   of S, where D(s,U) = s^2*M + s*Cs + K - A(s,U) is the dynamic stiffness
%   matrix (see AEROELASTIC_MATRIX). Its roots are the aeroelastic modes, as
%   the roots of det(sI-A) are the poles of a state-space model. Use it with
%   CROOTS and 'AssumeAnalytic',true, in a region that does not touch the
%   negative real axis (branch cut of Theodorsen's function), e.g.
%       model = aeroelastic_section('nasa');
%       croots(@(s) aeroelastic_delta(s,170,model),[1e-3 60 -250 250], ...
%           'AssumeAnalytic',true)          % unstable roots at 170 ft/s
%
%   V = AEROELASTIC_DELTA(S,U,MODEL,METHOD) uses METHOD = 'exact' (default),
%   'jones' or 'quasisteady' aerodynamics. 'pk' is refused: it is not an
%   analytic function of s (use AEROELASTIC_PK_ROOTS).
%
%   The determinant is computed for S^-1*D*S^-1 with S = diag(sqrt(diag(K))):
%   a constant scaling that improves conditioning and does not move the
%   roots. Do not replace it by abs(det(D)), norm(D) or singular values,
%   which are not analytic.
    if nargin<4, method='exact'; end
    if strcmpi(method,'pk')
        error('aeroelastic:Nonanalytic','Use a real two-variable solver for p-k.');
    end
    scale=sqrt(diag(model.K));
    value=arrayfun(@evaluate,s);
    function v=evaluate(z)
        d=aeroelastic_matrix(z,U,model,method)./(scale*scale.');
        v=det(d);
    end
end
