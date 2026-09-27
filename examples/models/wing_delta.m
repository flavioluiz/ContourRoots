function d = wing_delta(s,U,wing)
%WING_DELTA Characteristic function of the continuous wing, det K(s,U).
%   D = WING_DELTA(S,U,WING) evaluates det(WING_MATRIX(s,U,WING)) for every
%   element of S. Its zeros are the aeroelastic modes (poles) of the wing
%   at airspeed U, the analogue of det(sI - A) for a state-space model,
%   but for the CONTINUOUS structure with EXACT Theodorsen aerodynamics:
%   no finite elements, no modal truncation, no rational approximation.
%   Use it with CROOTS and 'AssumeAnalytic',true in a rectangle that does not
%   touch the negative real axis (the branch cut of Theodorsen's function):
%       wing = wing_model('goland');
%       croots(@(s) wing_delta(s,150,wing),[1e-3 60 -400 400],'AssumeAnalytic',true)
%   Do not replace it by abs(det), norm or singular values (not analytic).
    d=zeros(size(s),'like',s+1i);
    for j=1:numel(s), d(j)=det(full(wing_matrix(s(j),U,wing))); end
end
