function [w,alpha,z] = wing_mode_shape(lambda,U,wing,x,y)
%WING_MODE_SHAPE Continuous mode shape from a null vector of WING_MATRIX.
%   [W,ALPHA] = WING_MODE_SHAPE(LAMBDA,U,WING,X,Y) returns the deflection
%   w(y) and twist alpha(y) of the mode at LAMBDA, at the span positions Y
%   (0 <= y <= L). X is a null vector of WING_MATRIX(LAMBDA,U,WING,PIECES),
%   for example CMODES info.rightVectors{k}, for any number of PIECES.
%   Only the root state x(1:6) is used: it is propagated along the span with
%   the exact strip propagators, so the shape is exact between the nodes
%   too (no interpolation, no modal basis). Z returns all six states
%   [w; w'; alpha; V; M; T] in physical units. The shape is complex: its
%   phase describes how bending and twist move relative to each other.
    validateattributes(y,{'numeric'},{'vector','real','>=',0,'<=',wing.L*(1+1e-12)});
    z0=wing.scale.*x(1:6);                      % physical root state
    edges=[0 cumsum([wing.strips.length])];
    z=zeros(6,numel(y));
    for k=1:numel(y)
        zk=z0;
        for j=1:numel(wing.strips)
            len=min(y(k),edges(j+1))-edges(j);
            if len<=0, break, end
            piece=wing.strips(j); piece.length=len;
            zk=wing_propagator(lambda,U,piece,ones(6,1))*zk;
        end
        z(:,k)=zk;
    end
    w=z(1,:); alpha=z(3,:);
end
