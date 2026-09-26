function C = theodorsen_laplace(p)
%THEODORSEN_LAPLACE Generalized Theodorsen function C(p), p = s*b/U.
%   On the imaginary axis, p = i*k, it is the classical Theodorsen function
%   of the reduced frequency k. It is analytic except on the negative real
%   axis (branch cut).
%   C = K1(p)/(K0(p)+K1(p)), p=s*b/U. No rational approximation.
%   Negative real arguments are rejected: select a side of the branch cut.
%   C(0)=1 is the static limit, NOT an analytic extension through zero.
    validateattributes(p,{'numeric'},{'finite'});
    if any(imag(p(:))==0 & real(p(:))<0)
        error('theodorsen:BranchCut','Negative real arguments lie on the branch cut.');
    end
    C=ones(size(p)); nz=p~=0;
    % Both Bessel functions are multiplied by exp(p), which cancels.
    q=p(nz); k0=besselk(0,q,1); k1=besselk(1,q,1);
    C(nz)=1./(1+k0./k1);
end
