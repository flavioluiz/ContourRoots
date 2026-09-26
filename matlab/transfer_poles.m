function [p,info] = transfer_poles(G,region,varargin)
%TRANSFER_POLES Poles of a scalar transfer function after cancellations.
%   See COMPLEX_SPECTRUM. REGION = [xmin xmax ymin ymax].
    [p,info]=complex_spectrum(G,region,varargin{:},'Mode','poles');
end
