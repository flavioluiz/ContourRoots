function [poles,M,C,K] = coupled_beam_fem(parameters,nElements)
%COUPLED_BEAM_FEM Independent cubic-Hermite finite-element validation.
%   PARAMETERS is META.parameters returned by COUPLED_BEAM_MODEL.
%   Consistent element mass, clamped root, additional mass coordinate for
%   absorber topology. No modes from the analytical solver are used.
    validateattributes(nElements,{'numeric'},{'scalar','integer','>=',2});
    p=parameters; h=p.Length/nElements;
    ke=p.FlexuralRigidity/h^3*[12 6*h -12 6*h;6*h 4*h^2 -6*h 2*h^2; ...
        -12 -6*h 12 -6*h;6*h 2*h^2 -6*h 4*h^2];
    me=p.MassPerLength*h/420*[156 22*h 54 -13*h;22*h 4*h^2 13*h -3*h^2; ...
        54 13*h 156 -22*h;-13*h -3*h^2 -22*h 4*h^2];
    nd=2*(nElements+1); M=zeros(nd); K=M; C=M;
    for e=1:nElements
        idx=(2*e-1):(2*e+2); K(idx,idx)=K(idx,idx)+ke; M(idx,idx)=M(idx,idx)+me;
    end
    tip=nd-1;
    if strcmpi(p.Topology,'absorber')
        M(nd+1,nd+1)=p.Mass; K(nd+1,nd+1)=0; C(nd+1,nd+1)=0;
        v=zeros(nd+1,1); v(tip)=1; v(end)=-1;
        K=K+p.Stiffness*(v*v.'); C=C+p.Damping*(v*v.');
    else
        M(tip,tip)=M(tip,tip)+p.Mass;
        K(tip,tip)=K(tip,tip)+p.Stiffness;
        C(tip,tip)=C(tip,tip)+p.Damping;
    end
    % Delete exactly the two prescribed clamped degrees of freedom.
    M=M(3:end,3:end); C=C(3:end,3:end); K=K(3:end,3:end);
    n=size(M,1);
    poles=eig([zeros(n) eye(n);-K -C],[eye(n) zeros(n);zeros(n) M]);
end
