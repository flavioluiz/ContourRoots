function [y,tOut,info] = cstep(G,t,varargin)
%CSTEP Zero-state unit-step response without rationalizing the model.
%   [Y,T,INFO] = CSTEP(G,T,...) inverts G(s)/s directly.
%   G: handle, ndpair, scalar gain, symbolic scalar, or continuous SISO LTI.
%   For opaque models supply SingularityBound or Abscissa. Feedthrough=0
%   is the default assertion when no exact feedthrough metadata is known.
%   Options: Method ('fft'/'dehoog'/'quadrature'), AbsTol, RelTol, MaxPoints,
%   MaxRefinements, Feedthrough, Plot, Parent, Warn, Display.
%   No outputs plots. See CIMPULSE, CLSIM, CINVLAPLACE.
    [y,tOut,info]=response_run('step',G,[],t,varargin,nargout==0);
end
