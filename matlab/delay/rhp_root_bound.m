function R = rhp_root_bound(N,D)
%RHP_ROOT_BOUND Radius containing every root of D(s)+N(s)exp(-sT) with Re s>=0.
%   For Re(s) >= 0 and T >= 0, |exp(-sT)| <= 1, so any root satisfies
%   |D(s)| <= |N(s)|. If deg N < deg D (retarded case), this set is bounded:
%   with D normalized to be monic, |D(s)| >= |s|^m - sum|d_k||s|^k and
%   |N(s)| <= sum|n_k||s|^k. R is the positive root of the resulting
%   Cauchy-type polynomial, and the bound holds for every delay T >= 0.

    N = N(:).'; D = D(:).';
    N = N(find(N~=0,1):end); D = D(find(D~=0,1):end);
    m = numel(D)-1;
    assert(numel(N)-1 < m,'The bound requires deg N < deg D (retarded case).');
    N = N/D(1); D = D/D(1);
    c = -abs(D); c(1) = 1;                       % |s|^m - sum |d_k| |s|^k
    c(end-numel(N)+1:end) = c(end-numel(N)+1:end) - abs(N);
    r = roots(c);
    r = real(r(abs(imag(r)) < 1e-10 & real(r) > 0));
    R = max([r; 0]);
end
