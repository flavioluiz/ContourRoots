function [z,info] = czeros(G,region,varargin)
%CZEROS Zeros of a scalar transfer function inside a rectangle.
%   Z = CZEROS(G,REGION) returns the zeros of the transfer function G in
%   REGION = [xmin xmax ymin ymax], after pole-zero cancellations: a zero
%   of the numerator that is cancelled by an equal (or higher) order zero
%   of the denominator is not a zero of G.
%
%   Inputs and options are the same as for CPOLES. For a numerator/
%   denominator pair G = NDPAIR(N,D), CZEROS searches the zeros of N and
%   removes those shared with D. For a single function handle, CZEROS is
%   the same search as CROOTS.
%
%   [Z,INFO] = CZEROS(...) also returns the diagnostics described in
%   CROOTS and CPOLES.
%
%   Example
%       G = ndpair(@(s) (s+1).*exp(-s), [1 3 2]);   % (s+1)e^{-s}/((s+1)(s+2))
%       z = czeros(G,[-3 1 -3 3],'AssumeAnalytic',true)  % empty: s=-1 cancels
%
%   See also CPOLES, CROOTS, CPZMAP, NDPAIR, COMPLEX_SPECTRUM, ZERO.

    if nargin < 2, region = []; end
    [z,info] = cr_search('czeros','zeros',G,region,varargin,nargout<2);
end
