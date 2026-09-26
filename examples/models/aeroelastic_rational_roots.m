function r=aeroelastic_rational_roots(U,m,method)
%AEROELASTIC_RATIONAL_ROOTS Independent finite-state comparison models.
%   Jones uses two lag states for C(p)=.5+.165*.0455/(p+.0455)
%   +.335*.3/(p+.3). Quasi-steady uses C=1 and no aerodynamic states.
%   The exact-Theodorsen calculation never calls this function for seeds.
    validateattributes(U,{'numeric'},{'scalar','real','finite','nonnegative'});
    b=m.b; a=m.a;
    Ma=pi*m.rho*b^2*[1 -a*b;-a*b (1/8+a^2)*b^2];
    Me=m.M+Ma;
    Bnc=pi*m.rho*U*b*[0 -b;0 -(.5-a)*b^2];
    liftVector=2*pi*m.rho*U*b*[-1;b*(.5+a)];
    velocity=[1 b*(.5-a)]; position=[0 U];
    switch lower(method)
        case 'jones', c0=.5; beta=[.0455 .3]*U/b; weights=[.165 .335].*beta;
        case 'quasisteady', c0=1; beta=zeros(1,0); weights=zeros(1,0);
        otherwise, error('aeroelastic:Method','Use jones or quasisteady.');
    end
    nl=numel(beta);
    A=[zeros(2) eye(2) zeros(2,nl); ...
        Me\(-m.K+liftVector*c0*position) Me\(Bnc-m.Cs+liftVector*c0*velocity) Me\(liftVector*weights); ...
        ones(nl,1)*position ones(nl,1)*velocity -diag(beta)];
    r=eig(A);
end
