function G=wing_transfer_matrix(s,U,wing,assembly)
%WING_TRANSFER_MATRIX All 3-by-3 tip channels, implicit or hybrid assembly.
%   Scalar-node evaluator. Each physical strip is split into identical
%   numerical subpieces to moderate exponential growth (same heuristic as
%   WING_TRANSFER, not a universal conditioning bound). No modal truncation.
    if nargin<4, assembly='implicit'; end
    assembly=validatestring(assembly,{'implicit','hybrid'});
    validateattributes(s,{'numeric'},{'scalar','finite'});
    e=wing.strips; kb=max((abs(s)^2*[e.mu]./[e.EI]).^0.25);
    kt=max(abs(real(s))*sqrt([e.Ialpha]./[e.GJ]));
    pieces=max(1,ceil(max(kb,kt)*max([e.length])/3));
    if strcmp(assembly,'implicit')
        K=wing_matrix(s,U,wing,pieces); B=zeros(size(K,1),3);
        B(end-2:end,:)=diag(1./wing.scale(4:6)); X=K\B;
        G=diag(wing.scale(1:3))*X(end-5:end-3,:);
    else
        H=[];
        for j=1:numel(e)
            piece=e(j); piece.length=piece.length/pieces;
            [A,ia]=wing_hybrid(s,U,piece,wing.scale);
            if ~ia.available, error('wing:HybridPivot','Element hybrid pivot unavailable; use implicit assembly.'); end
            % Identical subpieces: associative connection by binary powering.
            % This is NOT a power of the numeric hybrid matrix. Every operation
            % still enforces the exact port compatibility/equilibrium equations.
            block=wing_hybrid_repeat(A,pieces);
            if isempty(H), H=block;
            else
                [H,ij]=wing_hybrid_join(H,block);
                if ~ij.available, error('wing:HybridPivot','Interface hybrid pivot unavailable; use implicit assembly.'); end
            end
        end
        G=H(4:6,4:6); % qL=0, fR=input; output qR
    end
end
