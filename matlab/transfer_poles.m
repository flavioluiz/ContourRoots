function [p,info] = transfer_poles(G,region,varargin)
%TRANSFER_POLES Poles of a scalar transfer function after cancellations.
%   Same as CPOLES, without the incompleteness warning. The name is kept
%   for compatibility with earlier code. REGION = [xmin xmax ymin ymax].
%
%   See also CPOLES, COMPLEX_SPECTRUM.
    [p,info]=complex_spectrum(G,region,varargin{:},'Mode','poles');
end
