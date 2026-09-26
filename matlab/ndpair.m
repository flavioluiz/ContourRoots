function G = ndpair(N,D)
%NDPAIR Transfer function given by a numerator and a denominator.
%   G = NDPAIR(N,D) represents G(s) = N(s)/D(s), where N and D are each
%   one of
%     - a function handle of one complex variable, e.g. @(s) sinh(s);
%     - a coefficient vector in descending powers, e.g. [1 3 2];
%     - a symbolic expression or symfun (Symbolic Math Toolbox).
%   Give N and D as separate ANALYTIC (entire, or at least holomorphic in
%   the search region) functions. Do not divide them yourself: keeping the
%   two factors separate is what lets CPOLES and CZEROS count poles and
%   zeros reliably and detect cancellations.
%
%   G is a struct with fields Numerator and Denominator. It is accepted by
%   CPOLES, CZEROS, CROOTS and CPZMAP. G can be evaluated with
%   G.Numerator(s)./G.Denominator(s) when both factors are handles.
%
%   Examples
%       G = ndpair(@(s) exp(-s), [1 1]);            % e^{-s}/(s+1)
%       G = ndpair(@(s) sinh(s/2), @(s) sinh(s));   % distributed system
%
%   See also CPOLES, CZEROS, CPZMAP.

    if nargin < 2
        error('ContourRoots:ndpair','NDPAIR needs a numerator and a denominator.');
    end
    check(N,'numerator'); check(D,'denominator');
    G = struct('Numerator',{N},'Denominator',{D});
end

function check(a,name)
    ok = isa(a,'function_handle') || isa(a,'sym') || isa(a,'symfun') || ...
        (isnumeric(a) && isvector(a) && ~isempty(a) && all(isfinite(a(:))));
    if ~ok
        error('ContourRoots:ndpair', ['The %s must be a function handle, a ' ...
            'coefficient vector or a symbolic expression.'], name);
    end
    if isnumeric(a) && all(a(:)==0)
        error('ContourRoots:ndpair','The %s must not be identically zero.',name);
    end
end
