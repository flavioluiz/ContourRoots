function H=wing_hybrid_repeat(A,n)
%WING_HYBRID_REPEAT Connect N identical elements with O(log N) joins.
%   Associative port connection, not A^N. All pivots are checked.
    validateattributes(n,{'numeric'},{'scalar','positive','integer','finite'});
    H=[];
    while n>0
        if mod(n,2)
            if isempty(H), H=A; else, H=join(H,A); end
        end
        n=floor(n/2);
        if n>0, A=join(A,A); end
    end
end
function C=join(A,B)
    [C,info]=wing_hybrid_join(A,B);
    if ~info.available, error('wing:HybridPivot','Repeated-element hybrid pivot unavailable; use implicit assembly.'); end
end
