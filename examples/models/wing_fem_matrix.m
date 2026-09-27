function D = wing_fem_matrix(s,U,fem,loads)
%WING_FEM_MATRIX Reduced dynamic stiffness s^2 M + K - A(s,U) of WING_FEM.
%   LOADS = 'laplace' (default) uses THEODORSEN_LOADS; 'hankel' uses the
%   independent harmonic formulas THEODORSEN_LOADS_HANKEL (s must then lie
%   on the imaginary axis). Validation oracle only.
    if nargin<4, loads='laplace'; end
    D=s^2*fem.M+fem.K;
    for j=1:numel(fem.strips)
        e=fem.strips(j);
        if strcmpi(loads,'hankel'), Q=theodorsen_loads_hankel(imag(s),U,e.b,e.a,e.rho);
        else, Q=theodorsen_loads(s,U,e.b,e.a,e.rho); end
        coeff=[Q(1,1) Q(1,2) Q(2,1) Q(2,2)];
        for k=1:4, D=D-coeff(k)*fem.W{j,k}; end
    end
end
