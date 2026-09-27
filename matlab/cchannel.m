function G = cchannel(M,output,input)
%CCHANNEL Extract a scalar evaluator from an explicit matrix transfer.
%   G = CCHANNEL(M,I,J) is a vectorized scalar function handle. Supply
%   scalar inversion-domain/feedthrough metadata when using it separately.
%   Matrix response calls retain channel metadata and share evaluations;
%   independent CCHANNEL handles do not share a cache.
    if ~isa(M,'ContourRootsModel')||~isscalar(M)
        error('ContourRoots:MatrixRepresentation','Use one CMIMO/CDYN model.');
    end
    validateattributes(output,{'numeric'},{'scalar','integer','>=',1,'<=',M.Size(1)});
    validateattributes(input,{'numeric'},{'scalar','integer','>=',1,'<=',M.Size(2)});
    G=@value;
    function v=value(s)
        pages=ceval(M,s); v=reshape(pages(output,input,:),size(s));
    end
end
