function [g,tOut,info] = cimpulse(G,t,varargin)
%CIMPULSE Zero-state unit-impulse response of a nonrational SISO model.
%   [G,T,INFO] = CIMPULSE(MODEL,T,...) returns the regular response.
%   Dirac terms are in INFO.singularTerms, not finite-height samples.
%   Opaque models require RegularImpulse=true and a valid inversion-domain
%   assertion (SingularityBound or Abscissa). InitialValue may specify g(0+).
%   Method: 'fft' (uniform T), 'dehoog' or 'quadrature' (nonuniform T too).
%   No outputs plots. See CSTEP, CLSIM, CINVLAPLACE.
    [g,tOut,info]=response_run('impulse',G,[],t,varargin,nargout==0);
end
