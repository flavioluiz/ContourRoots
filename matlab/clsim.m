function [y,tOut,info] = clsim(G,u,t,varargin)
%CLSIM Zero-state response to a prescribed real input, no ODE time march.
%   [Y,T,INFO] = CLSIM(G,U,T,...) uses integrated response kernels.
%   T is uniform, starts at zero; U is samples or a handle evaluated on T.
%   Interpolation: 'foh' (default, piecewise linear) or 'zoh' (held samples).
%   The supplied interpolant, not an unknown continuous signal, is simulated.
%   The third output is diagnostics, NOT internal states. No outputs plots.
%   See CSTEP for shared options and CINVLAPLACE for a known input transform.
    [y,tOut,info]=response_run('lsim',G,u,t,varargin,nargout==0);
end
