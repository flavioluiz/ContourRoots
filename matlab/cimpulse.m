function [g,tOut,info] = cimpulse(G,t,varargin)
%CIMPULSE Zero-state unit-impulse response of a nonrational SISO model.
%   [G,T,INFO] = CIMPULSE(MODEL,T,...) returns the regular response.
%   Dirac terms are in INFO.singularTerms, not finite-height samples.
%   Opaque models require RegularImpulse=true and a valid inversion-domain
%   assertion (SingularityBound or Abscissa). InitialValue may specify g(0+).
%   Method: 'fft' (uniform T), 'dehoog' or 'quadrature' (nonuniform T too).
%   No outputs plots. See CSTEP, CLSIM, CINVLAPLACE.
%   Explicit CMIMO/CDYN input returns Nt-by-ny-by-nu; singular terms are
%   tagged with output/input indices. Feedthrough/InitialValue are matrices.
    if isa(G,'ContourRootsModel')
        [g,tOut,info]=matrix_response_run('impulse',G,[],t,varargin,nargout==0);
    else
        [g,tOut,info]=response_run('impulse',G,[],t,varargin,nargout==0);
    end
end
