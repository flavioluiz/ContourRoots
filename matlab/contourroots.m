function contourroots()
%CONTOURROOTS Overview of the ContourRoots toolbox.
%   CONTOURROOTS prints the version, the main functions and where to find
%   the documentation.
%
%   ContourRoots finds roots of scalar analytic functions and poles and
%   zeros of scalar nonrational transfer functions (with exponentials,
%   hyperbolic functions, square roots after regularization, ...) inside a
%   rectangle of the complex plane. It evaluates the function as given,
%   without Padé or modal truncation, and checks completeness with the
%   argument principle.
%
%   Quick start
%       F = @(s) s.^2 + s + 1 + (2*s + 3).*exp(-s);
%       r = croots(F,[-8 2 -20 20],'AssumeAnalytic',true)
%
%       G = ndpair(@(s) sinh(s/2), @(s) sinh(s));
%       cpzmap(G,[-1 1 -10 10],'AssumeAnalytic',true)
%
%   See also CROOTS, CPOLES, CZEROS, CPZMAP, NDPAIR, CRITICAL_DELAYS.

    root = fileparts(fileparts(mfilename('fullpath')));
    fprintf('\nContourRoots %s - poles, zeros and roots of nonrational functions\n\n', ...
        contourroots_version());
    fprintf('  croots(F,region)    roots of an analytic function F\n');
    fprintf('  cpoles(G,region)    poles of a transfer function G\n');
    fprintf('  czeros(G,region)    zeros of a transfer function G\n');
    fprintf('  cpzmap(G,region)    pole-zero map\n');
    fprintf('  ndpair(N,D)         G = N/D from two analytic functions\n');
    fprintf('  critical_delays     stability crossings of D(s)+N(s)exp(-sT)\n\n');
    fprintf('  region = [xmin xmax ymin ymax] in the complex plane\n\n');
    fprintf('  help croots                 function reference\n');
    fprintf('  Tutorials and manual:       %s\n', fullfile(root,'docs'));
    fprintf('  Examples:                   %s\n\n', fullfile(root,'examples'));
end
