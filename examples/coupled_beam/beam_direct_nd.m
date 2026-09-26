%% Poles of a beam coupled to a mass-spring-damper, written directly as N/D
% Transfer function from a force F at the beam tip to the tip displacement
% y, for a clamped Euler-Bernoulli beam with a tuned mass absorber:
%
%            N(s)      Nb(s) (m s^2 + c s + k)
%   G(s) = ------ = -------------------------------------------
%            D(s)    Db(s) (m s^2 + c s + k) + Nb(s) m s^2 (c s + k)
%
% Nb/Db is the receptance of the bare beam. Normalized beam: L = EI = rho*A = 1.
% See docs/tutorials/08_coupled_beam.md for the derivation.
% No Symbolic Math Toolbox and no tf object are needed.
% Run setup_contourroots once per session before this script.

% Oscillator parameters (normalized)
m = 0.1;
k = 1.2362363368;
c = 0.05;

% Receptance of the bare beam, Hb(s) = Nb(s)/Db(s)
beta = @(s) (-s.^2).^(1/4);
Nb = @beam_numerator;                          % local function, end of file
Db = @(s) 1 + cosh(beta(s)).*cos(beta(s));

% Spring-damper element
Zc = @(s) k + c*s;

% Coupled transfer function: tip displacement / tip force
N = @(s) Nb(s).*(m*s.^2 + Zc(s));
D = @(s) Db(s).*(m*s.^2 + Zc(s)) + Nb(s).*m.*s.^2.*Zc(s);
G = ndpair(N, D);

% Search rectangle: [Re_min Re_max Im_min Im_max]
region = [-60 1 -150 150];

% Poles of the transfer function, after cancellations with the numerator
[p,info] = cpoles(G, region, 'AssumeAnalytic', true);
assert(info.complete, 'The search left unresolved regions.');

% Show one pole of each conjugate pair, by increasing frequency
upper = p(imag(p) > 0);
[~,order] = sort(imag(upper));
disp(upper(order))

% Alternative: roots of the characteristic equation D(s) = 0.
% Here both give the same ten poles; in general they can differ (see text).
[r,infoRoots] = croots(D, region, 'AssumeAnalytic', true);
assert(infoRoots.complete, 'The root search is incomplete.');

figure
cpzmap(G, [-1 0.2 -130 130], 'AssumeAnalytic', true)

% Local function: numerator of the beam receptance,
% A(q) = [cosh(b) sin(b) - sinh(b) cos(b)] / b^3 with b^4 = q = -s^2.
% It is an entire function of s; the series avoids 0/0 near the origin.
function A = beam_numerator(s)
    q = -s.^2;
    A = zeros(size(q));
    small = abs(q) < 1e-3;

    z = q(small);
    A(small) = 2/3 - z/315 + z.^2/623700 - z.^3/5108103000;

    b = q(~small).^(1/4);
    A(~small) = (cosh(b).*sin(b) - sinh(b).*cos(b))./b.^3;
end
