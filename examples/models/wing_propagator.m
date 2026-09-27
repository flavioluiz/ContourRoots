function E = wing_propagator(s,U,strip,scale)
%WING_PROPAGATOR Exact transfer matrix of one wing strip, z(right)=E*z(left).
%   Inside a strip of constant properties the Laplace-transformed equations
%   of bending, torsion and Theodorsen air loads,
%       EI w'''' + (s^2 mu) w + (s^2 S) alpha - Q11 w - Q12 alpha = 0
%      -GJ alpha'' + (s^2 S) w + (s^2 Ialpha) alpha - Q21 w - Q22 alpha = 0,
%   (Q = THEODORSEN_LOADS, ' = d/dy along the span) are six first-order
%   ODEs with constant coefficients for the state
%       z = [w; theta = w'; alpha; V = -M'; M = EI w''; T = GJ alpha'],
%   so z(y) = expm(A*y) z(0) EXACTLY: no spatial discretization. The state
%   is scaled by the fixed vector SCALE (z = diag(SCALE)*zbar), which only
%   improves conditioning. E is entire in s for a dry strip, and analytic
%   except on the negative real axis with Theodorsen aerodynamics.
    validateattributes(s,{'numeric'},{'scalar','finite'});
    Q=theodorsen_loads(s,U,strip.b,strip.a,strip.rho);
    D=s^2*[strip.mu strip.S;strip.S strip.Ialpha]-Q;
    A=zeros(6,'like',s+1i);
    A(1,2)=1;                  % w'     = theta
    A(2,5)=1/strip.EI;         % theta' = M/EI
    A(3,6)=1/strip.GJ;         % alpha' = T/GJ
    A(4,[1 3])=D(1,:);         % V'     = D11 w + D12 alpha
    A(5,4)=-1;                 % M'     = -V
    A(6,[1 3])=D(2,:);         % T'     = D21 w + D22 alpha
    A=(A.*scale.')./scale;     % scaled state
    E=expm(strip.length*A);
    if any(~isfinite(E(:)))
        error('wing:Overflow','Nonfinite strip propagator at s = %g%+gi.',real(s),imag(s));
    end
end
