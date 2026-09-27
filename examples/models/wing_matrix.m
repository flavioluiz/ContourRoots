function K = wing_matrix(s,U,wing,substeps)
%WING_MATRIX Boundary-value matrix of the clamped-free wing, K(s)*x = 0.
%   K = WING_MATRIX(S,U,WING) is the square matrix whose unknowns x are the
%   scaled states z = [w; w'; alpha; V; M; T] at the N+1 strip boundaries
%   (6 each). The rows state, in order:
%     clamped root:     w = w' = alpha = 0                 (3 equations)
%     each strip j:     z_(j+1) - E_j z_j = 0              (6 equations each)
%     free tip:         V = M = T = 0                      (3 equations)
%   K(s) is singular exactly at the modes of the wing: its determinant is the
%   characteristic function (WING_DELTA). Nothing is condensed or inverted,
%   so no artificial poles are created.
%
%   K = WING_MATRIX(S,U,WING,SUBSTEPS) splits every strip into SUBSTEPS equal
%   pieces (multiple shooting). For a strip of constant properties this is
%   exact, E = E_piece^SUBSTEPS, so the model is unchanged; it only keeps
%   each exponential moderate at high frequency (see WING_TRANSFER).
%   K is returned as a sparse matrix.
    if nargin<4, substeps=1; end
    n=numel(wing.strips); m=n*substeps; N=6*(m+1);
    props=cell2mat(struct2cell(wing.strips(:)));   % one column per strip
    I=[1 2 3]; J=[1 2 3]; V=ones(1,3);                  % clamped root
    for j=1:n
        if j==1 || any(props(:,j)~=props(:,j-1))       % identical strips share E
            piece=wing.strips(j); piece.length=piece.length/substeps;
            E=wing_propagator(s,U,piece,wing.scale);
            [ie,je,ve]=find(-E);
        end
        for k=1:substeps
            q=(j-1)*substeps+k;                         % piece number
            r0=3+6*(q-1); c0=6*(q-1);
            I=[I r0+ie.' r0+(1:6)]; J=[J c0+je.' c0+6+(1:6)]; %#ok<AGROW>
            V=[V ve.' ones(1,6)];                       %#ok<AGROW>
        end
    end
    I=[I N-2:N]; J=[J N-2:N]; V=[V ones(1,3)];          % free tip
    K=sparse(I,J,V,N,N);
end
