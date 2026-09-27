function M = cmimo(G,varargin)
%CMIMO Explicit nonrational transfer matrix; see ContourRootsModel.
%   M = CMIMO(G,'Size',[ny nu]) declares a matrix-valued handle.
%   M = CMIMO(CHANNELS) accepts a rectangular cell array of SISO models.
%   M = CMIMO(A) represents a constant numeric matrix, not polynomials.
%   See also CDYN, CEVAL, CCHANNEL, CLSIM, CKERNEL.
    M=ContourRootsModel('transfer',G,varargin{:});
end
