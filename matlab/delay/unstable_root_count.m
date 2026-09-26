function Z = unstable_root_count(N,D,T,nPerEdge)
%UNSTABLE_ROOT_COUNT Bound-based count of roots with Re(s) >= 0.
%   Z = UNSTABLE_ROOT_COUNT(N,D,T) counts the roots of
%       D(s) + N(s) exp(-s T) = 0
%   in the closed right half-plane, for polynomial coefficient vectors N, D
%   (descending powers) with deg N < deg D (retarded case).
%
%   RHP_ROOT_BOUND gives a radius R such that every root with Re(s) >= 0
%   satisfies |s| <= R, for every delay T >= 0. This bound is a mathematical
%   fact. The argument principle is then evaluated numerically on the
%   rectangle [0,R+1] x [-(R+1),R+1], so no search window has to be guessed
%   and no rational approximation is used. The count itself is a
%   floating-point contour computation, not an interval-arithmetic proof;
%   it is unreliable when a root lies on or very near Re(s) = 0, i.e. at or
%   near a critical delay.
%
%   Z = UNSTABLE_ROOT_COUNT(N,D,T,NPEREDGE) sets the number of samples per
%   rectangle edge (default 4000).
%
%   See also RHP_ROOT_BOUND, CRITICAL_DELAYS, DELAY_ROOT_COUNT.

    if nargin < 4, nPerEdge = 4000; end
    R = rhp_root_bound(N,D) + 1;
    Z = delay_root_count(N,D,T,[0 R],[-R R],nPerEdge);
end
