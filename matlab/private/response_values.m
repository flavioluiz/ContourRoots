function v=response_values(f,s)
% Scalar-only handles supported; bad shapes are never silently discarded.
    try
        v=f(s);
        if isscalar(v) && numel(s)>1, v=arrayfun(f,s); end
        if ~isequal(size(v),size(s)), v=arrayfun(f,s); end
    catch firstError
        if startsWith(firstError.identifier,'ContourRoots:Matrix'), rethrow(firstError); end
        try
            v=arrayfun(f,s);
        catch err
            if startsWith(err.identifier,'ContourRoots:Matrix'), rethrow(err); end
            error('ContourRoots:ResponseEvaluation','Transfer evaluation failed: %s',err.message);
        end
    end
    if ~isnumeric(v) || ~isequal(size(v),size(s)) || any(~isfinite(v(:)))
        error('ContourRoots:ResponseEvaluation','Transfer evaluations must be finite numeric scalars (check poles, cancellations and overflow).');
    end
    v=double(v);
end
