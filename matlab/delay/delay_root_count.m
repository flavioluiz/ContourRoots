function [count,winding] = delay_root_count(N,D,T,xlimits,ylimits,nPerEdge)
%DELAY_ROOT_COUNT Count zeros in a rectangle using the argument principle.
%   The quasi-polynomial is entire, so its boundary winding number is the
%   number of enclosed zeros. A root very near the boundary makes any
%   finite-sampling count ill-conditioned; Winding reports the raw value.

    if nargin < 6, nPerEdge = 500; end
    nPerEdge = max(30,round(nPerEdge));
    xb = linspace(xlimits(1),xlimits(2),nPerEdge);
    yr = linspace(ylimits(1),ylimits(2),nPerEdge);
    z = [xb + 1i*ylimits(1), ...
         xlimits(2) + 1i*yr(2:end), ...
         fliplr(xb(1:end-1)) + 1i*ylimits(2), ...
         xlimits(1) + 1i*fliplr(yr(2:end-1))];
    f = polyval(D,z) + exp(-z*T).*polyval(N,z);
    winding = sum(angle(f(2:end)./f(1:end-1))) + angle(f(1)/f(end));
    winding = winding/(2*pi);
    count = max(0,round(winding));
end
