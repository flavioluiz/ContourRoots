function g = wing_transfer(s,U,wing,output,input)
%WING_TRANSFER Tip transfer function of the continuous wing, G(s) = Y(s)/F(s).
%   G = WING_TRANSFER(S,U,WING,OUTPUT,INPUT) returns, for every element of S,
%   the tip response OUTPUT to a tip load INPUT:
%       INPUT  1 = shear force (N), 2 = bending moment (N m), 3 = torque (N m)
%       OUTPUT 1 = deflection w (m), 2 = slope w' (rad), 3 = twist alpha (rad).
%   The load enters the free-tip equations of WING_MATRIX and K(s)x = f is
%   solved directly: no modal expansion, no rational fit. Zero feedthrough.
%   At high |s| each strip is split into more pieces (multiple shooting),
%   which leaves the model unchanged and keeps the solve well conditioned.
%   Use it with CSTEP/CLSIM (with a justified 'SingularityBound').
    validateattributes(output,{'numeric'},{'scalar','integer','>=',1,'<=',3});
    validateattributes(input,{'numeric'},{'scalar','integer','>=',1,'<=',3});
    g=zeros(size(s),'like',s+1i);
    e=wing.strips; longest=max([e.length]);
    for j=1:numel(s)
        % Multiple shooting: keep the growth exp(rate*piece length) of every
        % solution below about exp(3) across a piece. Bending waves always
        % have growing/decaying components (rate = wave number); torsion
        % waves grow only with Re(s). Exact for constant-property strips.
        kb=max((abs(s(j))^2*[e.mu]./[e.EI]).^0.25);
        kt=max(abs(real(s(j)))*sqrt([e.Ialpha]./[e.GJ]));
        substeps=max(1,ceil(max(kb,kt)*longest/3));
        K=wing_matrix(s(j),U,wing,substeps); f=zeros(size(K,1),1);
        f(end-3+input)=1/wing.scale(3+input);
        x=K\f;
        g(j)=wing.scale(output)*x(end-6+output);
    end
end
