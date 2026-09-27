function [H,info]=wing_hybrid(s,U,strip,scale)
%WING_HYBRID Physical ports [fL;qR]=H*[qL;fR] of one constant strip.
%   q=[w;w';alpha], p=[V;M;T], fL=-pL and fR=pR. Port forces are outward
%   boundary loads in these fixed coordinates. The same SCALE as WING_MATRIX
%   conditions the internal computation; H maps physical SI port variables.
%   A singular Epp is an unavailable hybrid chart, not an invalid implicit
%   element: INFO.implicit remains usable. No pivot regularization is done.
    E=wing_propagator(s,U,strip,scale);
    a=E(1:3,1:3); b=E(1:3,4:6); c=E(4:6,1:3); d=E(4:6,4:6);
    reciprocalCondition=rcond(d);
    info=struct('available',reciprocalCondition>1e-12,'rcond',reciprocalCondition, ...
        'propagator',E,'implicit',[-E eye(6)]);
    H=NaN(6); if ~info.available, return; end
    X=d\[c eye(3)];
    Hs=[X(:,1:3) -X(:,4:6);a-b*X(:,1:3) b*X(:,4:6)];
    H=(Hs.*scale([4:6 1:3]))./scale.';
end
