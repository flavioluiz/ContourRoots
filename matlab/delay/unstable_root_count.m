function [Z,info] = unstable_root_count(N,D,T,nPerEdge,maxPerEdge)
%UNSTABLE_ROOT_COUNT Number of roots with Re(s) > 0, by a bound-based count.
%   Z = UNSTABLE_ROOT_COUNT(N,D,T) counts the roots of
%       D(s) + N(s) exp(-s T) = 0
%   with Re(s) > 0, for polynomial coefficient vectors N, D (descending
%   powers) with deg N < deg D (retarded case).
%
%   RHP_ROOT_BOUND gives a radius R such that every root with Re(s) >= 0
%   satisfies |s| <= R, for every delay T >= 0 (a mathematical fact). The
%   argument principle is then evaluated on the rectangle
%   [0,R+1] x [-(R+1),R+1] with adaptive sampling (DELAY_ROOT_COUNT), so no
%   search window has to be guessed and no rational approximation is used.
%
%   The left edge of the rectangle is the imaginary axis. The count is
%   therefore INCONCLUSIVE (Z = NaN, with a warning) when a root lies on or
%   extremely close to the imaginary axis - at or near a critical delay, for
%   persistent imaginary roots, or for large delays, where many roots
%   approach the axis - or when the sampling budget is not enough. A count
%   that could not be resolved is never reported as zero. Z == 0 means
%   stability only when there is no root on the imaginary axis.
%
%   The count is a floating-point computation, not an interval-arithmetic
%   proof. For many delays, the crossing formula based on CRITICAL_DELAYS
%   is exact and much cheaper; this function is an independent check.
%
%   [Z,INFO] = UNSTABLE_ROOT_COUNT(...) also returns the diagnostics of
%   DELAY_ROOT_COUNT (resolved, samplesPerEdge, reason) and does not warn.
%   NPEREDGE (default 256) and MAXPEREDGE (default 2^22) set the initial and
%   maximum samples per edge.
%
%   See also RHP_ROOT_BOUND, CRITICAL_DELAYS, DELAY_ROOT_COUNT.

    if nargin < 4, nPerEdge = []; end
    if nargin < 5, maxPerEdge = []; end
    R = rhp_root_bound(N,D) + 1;
    [Z,~,info] = delay_root_count(N,D,T,[0 R],[-R R],nPerEdge,maxPerEdge);
    info.radius = R - 1;
    if isnan(Z) && nargout < 2
        warning('unstable_root_count:Inconclusive', ...
            'Count inconclusive for T = %g (%s). A root may lie on or near the imaginary axis.', ...
            T, info.reason);
    end
end
