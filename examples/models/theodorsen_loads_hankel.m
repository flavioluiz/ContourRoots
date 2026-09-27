function A = theodorsen_loads_hankel(omega,U,b,a,rho)
%THEODORSEN_LOADS_HANKEL Harmonic strip loads from the original lift/moment formulas.
%   Independent check of THEODORSEN_LOADS on the imaginary axis s = i*OMEGA:
%   Theodorsen's function is written with Hankel functions of the reduced
%   frequency k = OMEGA*B/U, and the lift and moment are assembled from
%   Theodorsen's original expressions (not from THEODORSEN_LAPLACE).
%   Same convention: [downward force; nose-up moment] = A*[h; alpha].
    if U<=0, error('wing:Harmonic','A positive airspeed is required.'); end
    s=1i*omega; k=omega*b/U;
    C=besselh(1,2,k)/(besselh(1,2,k)+1i*besselh(0,2,k));
    A=zeros(2);
    for j=1:2
        q=zeros(2,1); q(j)=1; h=q(1); alpha=q(2);
        downwash=s*h+(U+b*(.5-a)*s)*alpha;       % 3/4-chord normal velocity
        lift=pi*rho*b^2*(s^2*h+U*s*alpha-a*b*s^2*alpha) ...
            +2*pi*rho*U*b*C*downwash;
        moment=pi*rho*b^2*(a*b*s^2*h-U*b*(.5-a)*s*alpha ...
            -b^2*(1/8+a^2)*s^2*alpha)+2*pi*rho*U*b^2*(.5+a)*C*downwash;
        A(:,j)=[-lift;moment];
    end
end
