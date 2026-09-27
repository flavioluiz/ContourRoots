function [K,info] = ckernel(G,t,varargin)
%CKERNEL Prepare reusable step/ramp kernels for repeated CLSIM calls.
%   K = CKERNEL(G,T,...) evaluates the transfer function only during setup.
%   [Y,T,INFO] = CLSIM(K,U,T,...) then simulates another input without any
%   new transfer evaluations or inverse transforms. T and Interpolation
%   must match the preparation; K is a read-only value object, not a model.
%   FOH stores step AND ramp kernels, so later inputs may have U(1) ~= 0.
%   ZOH stores the step kernel only. Direct terms and delays are retained.
%
%   Preparation AbsTol/RelTol apply to kernels (defaults 1e-8 and 1e-6).
%   CLSIM still checks its own output tolerances for each input (defaults
%   1e-6 and 1e-4). Large inputs can fail that check: prepare tighter kernels
%   explicitly. Unresolved preparation raises ContourRoots:KernelUnresolved.
%   [K,INFO] also returns K.Info, including setup evaluation counts.
%
%   Example
%       t = (0:.02:4).';
%       K = ckernel(@(s) 1./(s+1+.5*exp(-s)),t,'SingularityBound',0);
%       [y,~,info] = clsim(K,sin(t),t);
%       assert(info.evaluations == 0 && info.kernelReused);
%
%   See also CLSIM, CSTEP, CINVLAPLACE.
%   Explicit CMIMO/CDYN input prepares every channel in a separate immutable
%   ContourRootsMatrixKernel bank. Shared matrix evaluations occur only here.
    if isa(G,'ContourRootsModel')
        K=ContourRootsMatrixKernel(G,t,varargin{:});
    else
        K=ContourRootsKernel(G,t,varargin{:});
    end
    info=K.Info;
end
