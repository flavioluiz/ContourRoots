function [p,info] = cpoles(G,region,varargin)
%CPOLES Poles of a scalar transfer function inside a rectangle.
%   P = CPOLES(G,REGION) returns the poles of the transfer function G in
%   REGION = [xmin xmax ymin ymax]. Pole-zero cancellations are removed:
%   a zero of the denominator that is cancelled by an equal (or higher)
%   order zero of the numerator is not a pole of G. P is a column vector
%   ordered by decreasing real part.
%
%   The recommended input is a numerator/denominator pair of analytic
%   functions, G = NDPAIR(N,D). It also accepts a symbolic quotient, a
%   continuous-time SISO tf/ss/zpk model, or a single function handle.
%   With a single handle G, poles are searched as zeros of 1./G; this mode
%   is exploratory unless 'AssumeAnalytic',true, which then asserts that
%   1./G is analytic, i.e. that G has NO zeros in the rectangle.
%
%   [P,INFO] = CPOLES(...) also returns diagnostics. In addition to the
%   fields described in CROOTS:
%     INFO.cancelledLocations  denominator zeros removed by cancellation
%     INFO.cancellationOrders  how many orders were cancelled at each one
%
%   Options are the same as for CROOTS ('AssumeAnalytic', 'Plot',
%   'Display', 'Warn', and the advanced COMPLEX_SPECTRUM options).
%
%   Example
%       % Poles of sinh(s/2)/sinh(s): the even multiples of i*pi cancel
%       G = ndpair(@(s) sinh(s/2), @(s) sinh(s));
%       p = cpoles(G,[-1 1 -10 10],'AssumeAnalytic',true)
%       % p = i*pi*[3 1 -1 -3] (in some order)
%
%   See also CZEROS, CROOTS, CPZMAP, NDPAIR, COMPLEX_SPECTRUM, POLE.

    if nargin < 2, region = []; end
    [p,info] = cr_search('cpoles','poles',G,region,varargin,nargout<2);
end
