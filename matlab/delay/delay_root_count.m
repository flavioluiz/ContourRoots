function [count,winding,info] = delay_root_count(N,D,T,xlimits,ylimits,nPerEdge,maxPerEdge)
%DELAY_ROOT_COUNT Argument-principle count of the roots of D(s)+N(s)exp(-sT).
%   COUNT = DELAY_ROOT_COUNT(N,D,T,XLIM,YLIM) counts the roots of the
%   quasi-polynomial D(s) + N(s) exp(-s T) inside the rectangle
%   XLIM(1) < Re(s) < XLIM(2), YLIM(1) < Im(s) < YLIM(2).
%
%   The boundary is sampled adaptively: the number of samples per edge
%   starts at NPEREDGE (default 256) and is doubled, up to MAXPEREDGE
%   (default 2^22), until the rounded winding number is the same at three
%   consecutive levels and, at the last two, every phase increment is below
%   pi/3, every log-modulus increment below 2, and the winding number is
%   within 1e-7 of an integer.
%
%   If this does not happen (a root on or very close to the boundary, or a
%   function too oscillatory for the sampling budget), COUNT is NaN: the
%   count is INCONCLUSIVE. It is never replaced by zero.
%
%   [COUNT,WINDING,INFO] = DELAY_ROOT_COUNT(...) also returns the last
%   winding number and a structure with fields resolved, samplesPerEdge
%   and reason.
%
%   See also UNSTABLE_ROOT_COUNT, DELAY_ROOTS.

    if nargin < 6 || isempty(nPerEdge), nPerEdge = 256; end
    if nargin < 7 || isempty(maxPerEdge), maxPerEdge = 2^22; end
    nPerEdge = max(32,round(nPerEdge));
    box = [xlimits(:).' ylimits(:).'];
    count = NaN; winding = NaN; previous = NaN; stable = 0;
    info = struct('resolved',false,'samplesPerEdge',0,'reason','');
    n = nPerEdge;
    while n <= maxPerEdge
        [w,maxPhase,maxMag,ok] = boundary_winding(N,D,T,box,n);
        info.samplesPerEdge = n;
        if ~ok
            info.reason = 'zero, overflow or NaN on the boundary';
            return
        end
        winding = w;
        c = round(w);
        resolved = maxPhase < pi/3 && maxMag < 2 && abs(w-c) < 1e-7;
        if c == previous && resolved
            stable = stable + 1;
        else
            stable = 0;
        end
        previous = c;
        if stable >= 2
            count = c;
            info.resolved = true;
            return
        end
        n = 2*n;
    end
    info.reason = sprintf('not resolved with %d samples per edge',info.samplesPerEdge);
end

function [w,maxPhase,maxMag,ok] = boundary_winding(N,D,T,box,n)
% Winding number of F along the counterclockwise boundary, evaluated in
% chunks so that millions of samples are never held in memory at once.
    total = 4*n; chunk = 2^20;
    w = 0; maxPhase = 0; maxMag = 0; ok = true;
    first = []; last = [];
    for start = 0:chunk:total-1
        j = start:min(start+chunk,total)-1;
        f = evaluate(N,D,T,point(j,n,box));
        if any(~isfinite(f) | f == 0), ok = false; return; end
        if isempty(first), first = f(1); end
        if ~isempty(last), f = [last f]; end %#ok<AGROW>
        [w,maxPhase,maxMag] = accumulate(f,w,maxPhase,maxMag);
        last = f(end);
    end
    [w,maxPhase,maxMag] = accumulate([last first],w,maxPhase,maxMag);
    w = w/(2*pi);
end

function [w,maxPhase,maxMag] = accumulate(f,w,maxPhase,maxMag)
    ratio = f(2:end)./f(1:end-1);
    dphase = angle(ratio);
    w = w + sum(dphase);
    maxPhase = max([maxPhase abs(dphase)]);
    maxMag = max([maxMag abs(log(abs(ratio)))]);
end

function z = point(j,n,box)
% Samples j = 0..4n-1 of the boundary: bottom, right, top, left.
    edge = floor(j/n); t = mod(j,n)/n;
    x0 = box(1); x1 = box(2); y0 = box(3); y1 = box(4);
    z = complex(zeros(size(j)));
    e = edge == 0; z(e) = x0 + (x1-x0)*t(e) + 1i*y0;
    e = edge == 1; z(e) = x1 + 1i*(y0 + (y1-y0)*t(e));
    e = edge == 2; z(e) = x1 - (x1-x0)*t(e) + 1i*y1;
    e = edge == 3; z(e) = x0 + 1i*(y1 - (y1-y0)*t(e));
end

function f = evaluate(N,D,T,z)
    f = polyval(D,z) + exp(-z*T).*polyval(N,z);
end
