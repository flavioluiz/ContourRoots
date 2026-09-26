function [r,info] = croots(F,region,varargin)
%CROOTS Roots of a scalar analytic function inside a rectangle.
%   R = CROOTS(F,REGION) returns the roots (zeros) of F inside the
%   rectangle REGION = [xmin xmax ymin ymax] of the complex plane, that is
%   xmin < Re(s) < xmax and ymin < Im(s) < ymax. R is a column vector,
%   ordered by decreasing real part. A root of multiplicity m is listed
%   once; see INFO.multiplicity.
%
%   F can be
%     - a function handle, for example @(s) s.^2 + s + 1 + exp(-s)
%       (vectorized with .* ./ .^ is fastest, but not required);
%     - a coefficient vector in descending powers, as in ROOTS;
%     - a symbolic expression or symfun (Symbolic Math Toolbox);
%     - a numerator/denominator pair from NDPAIR (roots of the numerator
%       that are not cancelled by the denominator);
%     - a continuous-time SISO tf, ss or zpk model (Control System Toolbox).
%
%   [R,INFO] = CROOTS(...) also returns diagnostics. The most important
%   fields are
%     INFO.status        'numerically_complete', 'unresolved' or 'exploratory'
%     INFO.complete      true when the contour counts account for every root
%     INFO.multiplicity  multiplicity of each root in R
%     INFO.residuals     |F| at each root
%   See the documentation (docs/api/croots.md) for every field.
%
%   CROOTS(...,Name,Value) sets options. The most useful are
%     'AssumeAnalytic'  true declares that a function handle is analytic
%                       (holomorphic) on a neighborhood of the closed
%                       rectangle. This enables the completeness check.
%                       Default false for handles; symbolic, polynomial
%                       and NDPAIR inputs are checked automatically.
%     'Plot'            true draws the roots in the rectangle (false).
%     'Display'         true prints a table of results (false).
%     'Warn'            false silences the warning issued when the search
%                       is not numerically complete and INFO is not
%                       requested (true).
%   Advanced options (Derivative, SeedPoints, RootTolerance, ContourPoints,
%   ContourRefinements, MaxDepth, MaxCells, MaxIterations, GridSize,
%   Singularities) are described in COMPLEX_SPECTRUM.
%
%   The function is evaluated exactly as given: no Padé or other rational
%   approximation is used. Roots are floating-point numbers, and
%   "numerically complete" is a numerical check, not a formal proof.
%
%   Examples
%       % Characteristic equation of a system with a delay
%       F = @(s) s.^2 + s + 1 + (2*s + 3).*exp(-s);
%       [r,info] = croots(F,[-8 2 -20 20],'AssumeAnalytic',true);
%       info.status
%
%       % Same interface as ROOTS for polynomials
%       croots([1 0 -1],[-2 2 -1 1])        % returns 1 and -1
%
%   See also CPOLES, CZEROS, CPZMAP, NDPAIR, COMPLEX_SPECTRUM, ROOTS.

    if nargin < 2, region = []; end
    [r,info] = cr_search('croots','zeros',F,region,varargin,nargout<2);
end
