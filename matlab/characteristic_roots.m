function [r,info] = characteristic_roots(Delta,region,varargin)
%CHARACTERISTIC_ROOTS Zeros of a scalar characteristic expression.
%   See COMPLEX_SPECTRUM. REGION = [xmin xmax ymin ymax].
    [r,info]=complex_spectrum(Delta,region,varargin{:},'Mode','zeros');
end
