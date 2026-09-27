function [H,info]=wing_hybrid_join(A,B)
%WING_HYBRID_JOIN Eliminate the interface of two physical hybrid elements.
%   qR_A=qL_B and fR_A+fL_B=0, NOT H=A*B. With blocks ordered as
%   [fL;qR]=H*[qL;fR], the interface coordinate solves
%   (I+A22*B11)qI = A21*qL - A22*B12*fR.
%   An unavailable pivot is reported without regularization.
    validateattributes(A,{'numeric'},{'size',[6 6],'finite'});
    validateattributes(B,{'numeric'},{'size',[6 6],'finite'});
    K=eye(3)+A(4:6,4:6)*B(1:3,1:3);
    reciprocalCondition=rcond(K);
    info=struct('available',reciprocalCondition>1e-12,'rcond',reciprocalCondition);
    H=NaN(6); if ~info.available, return; end
    X=K\[A(4:6,1:3) -A(4:6,4:6)*B(1:3,4:6)];
    FL=[A(1:3,1:3) zeros(3)]-A(1:3,4:6)*(B(1:3,1:3)*X+[zeros(3) B(1:3,4:6)]);
    QR=B(4:6,1:3)*X+[zeros(3) B(4:6,4:6)];
    H=[FL;QR]; info.interfaceMap=X;
end
