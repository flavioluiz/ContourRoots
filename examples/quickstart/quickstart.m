%% ContourRoots quick start
% Run setup_contourroots once per MATLAB session before this script.
% Every block below can also be pasted into the Command Window.

%% 1. Roots of a characteristic equation with a delay
% F(s) = s^2 + s + 1 + (2s + 3) e^{-s} has infinitely many roots. We look
% for the ones in the rectangle -8 < Re(s) < 2, -20 < Im(s) < 20.
F = @(s) s.^2 + s + 1 + (2*s + 3).*exp(-s);
region = [-8 2 -20 20];                 % [xmin xmax ymin ymax]
[r,info] = croots(F, region, 'AssumeAnalytic', true);
disp(r)                                 % 7 roots, rightmost first
disp(info.status)                       % 'numerically_complete'

% 'AssumeAnalytic',true states that F has no poles, branch cuts or other
% singularities in (a neighborhood of) the rectangle. That is true here:
% F is built from polynomials and exp. ContourRoots cannot verify it for
% an arbitrary function handle, so you declare it.

%% 2. The same interface as ROOTS for polynomials
croots([1 0 -1], [-2 2 -1 1])           % roots of s^2 - 1: 1 and -1

%% 3. Poles and zeros of a transfer function
% Give the numerator and denominator separately with NDPAIR. Here
% G(s) = sinh(s/2)/sinh(s): sinh(s) vanishes at i*pi*k, but at the even k
% the numerator vanishes too and the pole cancels.
G = ndpair(@(s) sinh(s/2), @(s) sinh(s));
p = cpoles(G, [-1 1 -10 10], 'AssumeAnalytic', true)   % i*pi*(+-1, +-3)
z = czeros(G, [-1 1 -10 10], 'AssumeAnalytic', true)   % empty

%% 4. Pole-zero map
figure
cpzmap(G, [-1 1 -10 10], 'AssumeAnalytic', true)
