function [y,tOut,info] = clsim(G,u,t,varargin)
%CLSIM Zero-state response to a prescribed real input, no ODE time march.
%   [Y,T,INFO] = CLSIM(G,U,T,...) uses integrated response kernels.
%   T is uniform, starts at zero; U is samples or a handle evaluated on T.
%   Interpolation: 'foh' (default, piecewise linear) or 'zoh' (held samples).
%   The supplied interpolant, not an unknown continuous signal, is simulated.
%   The third output is diagnostics, NOT internal states. No outputs plots.
%   CLSIM(K,U,T,...) reuses K=CKERNEL(G,T,...), with no new inversions.
%   The same grid/hold is required; output error is rechecked for every input.
%   See CSTEP for shared options and CINVLAPLACE for a known input transform.
%   Explicit CMIMO/CDYN models take Nt-by-nu inputs and return Nt-by-ny.
%   Errors are summed per output, including cancellation between channels.
    if isa(G,'ContourRootsModel')
        [y,tOut,info]=matrix_response_run('lsim',G,u,t,varargin,nargout==0);
    else
        [y,tOut,info]=response_run('lsim',G,u,t,varargin,nargout==0);
    end
end
