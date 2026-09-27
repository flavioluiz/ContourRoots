function M = cdyn(H,B,C,D,varargin)
%CDYN Structured nonrational model H(s)x=B(s)u, y=C(s)x+D(s)u.
%   M = CDYN(H,B,C,D,'Dimensions',[n ny nu]) retains the original factors.
%   Factors are numeric matrices, symbolic matrices or matrix handles.
%   CEVAL uses C*(H\B)+D, one block solve per node, never INV.
%   Dimensions may be inferred from constants; opaque handles are not probed.
%   Batched=true declares handle layouts rows-by-columns-by-numberOfNodes.
%   D is a model factor, NOT necessarily the high-frequency feedthrough;
%   declare Feedthrough separately if needed by the time-response adapter.
%   This interface does not implement matrix poles or transmission zeros.
%   See also CMIMO, CEVAL, CCHANNEL, CKERNEL.
    M=ContourRootsModel('structured',{H,B,C,D},varargin{:});
end
