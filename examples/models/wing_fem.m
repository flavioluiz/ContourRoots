function fem = wing_fem(wing,elementsPerStrip,nModes)
%WING_FEM Finite-element model of the wing, used ONLY as an independent check.
%   FEM = WING_FEM(WING,NE,NMODES) discretizes each strip of WING into NE
%   beam elements (cubic Hermite bending, linear Saint-Venant torsion),
%   clamps the root, and projects onto the first NMODES vacuum modes. The
%   air loads are integrated with 4-point Gauss quadrature. This is the
%   classical discretize-then-truncate route that the continuous model of
%   WING_DELTA avoids; it converges to the same answer as NE and NMODES grow.
%   Use WING_FEM_MATRIX to evaluate the reduced dynamic stiffness matrix.
    if nargin<2, elementsPerStrip=4; end
    if nargin<3, nModes=20; end
    ne=numel(wing.strips)*elementsPerStrip; nd=3*(ne+1);
    K=zeros(nd); M=K; W=cell(numel(wing.strips),4);
    for j=1:numel(W), W{j}=zeros(nd); end
    xg=[-.861136311594053 -.339981043584856 .339981043584856 .861136311594053];
    wg=[.347854845137454 .652145154862546 .652145154862546 .347854845137454];
    for j=1:ne
        strip=ceil(j/elementsPerStrip); e=wing.strips(strip);
        l=e.length/elementsPerStrip; ix=3*(j-1)+(1:6); bend=[1 2 4 5]; twist=[3 6];
        kl=zeros(6); kl(bend,bend)=e.EI/l^3*[12 6*l -12 6*l; ...
            6*l 4*l^2 -6*l 2*l^2;-12 -6*l 12 -6*l;6*l 2*l^2 -6*l 4*l^2];
        kl(twist,twist)=e.GJ/l*[1 -1;-1 1]; K(ix,ix)=K(ix,ix)+kl;
        for g=1:4
            r=(xg(g)+1)/2; weight=wg(g)*l/2;
            N=[1-3*r^2+2*r^3 l*(r-2*r^2+r^3) 0 3*r^2-2*r^3 l*(-r^2+r^3) 0; ...
               0 0 1-r 0 0 r];
            M(ix,ix)=M(ix,ix)+weight*N.'*[e.mu e.S;e.S e.Ialpha]*N;
            c=0;
            for p=1:2
                for q=1:2
                    c=c+1; W{strip,c}(ix,ix)=W{strip,c}(ix,ix)+weight*N(p,:).'*N(q,:);
                end
            end
        end
    end
    free=4:nd; K=K(free,free); M=M(free,free);
    [V,d]=eig(K,M,'vector'); [d,order]=sort(real(d));
    V=V(:,order(1:min(nModes,numel(order))));
    for j=1:size(V,2), V(:,j)=V(:,j)/sqrt(V(:,j)'*M*V(:,j)); end
    for j=1:numel(W), W{j}=V.'*W{j}(free,free)*V; end
    fem=struct('K',V.'*K*V,'M',V.'*M*V,'W',{W},'omega',sqrt(d), ...
        'strips',wing.strips,'nElements',ne,'nModes',size(V,2));
end
