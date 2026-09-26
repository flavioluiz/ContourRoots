function [r,info] = characteristic_roots(Delta,region,varargin)
%CHARACTERISTIC_ROOTS Zeros of a scalar characteristic expression.
%   Same as CROOTS, without the incompleteness warning. The name is kept
%   for compatibility with earlier code. REGION = [xmin xmax ymin ymax].
%
%   See also CROOTS, COMPLEX_SPECTRUM.
    [r,info]=complex_spectrum(Delta,region,varargin{:},'Mode','zeros');
end
