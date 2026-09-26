function v=response_values(f,s)
% Scalar-only handles supported; bad shapes are never silently discarded.
    try
        v=f(s);
        if isscalar(v) && numel(s)>1, v=arrayfun(f,s); end
        if ~isequal(size(v),size(s)), v=arrayfun(f,s); end
    catch
        try
            v=arrayfun(f,s);
        catch err
            error('ContourRoots:ResponseEvaluation','Transfer evaluation failed: %s',err.message);
        end
    end
    if ~isnumeric(v) || ~isequal(size(v),size(s)) || any(~isfinite(v(:)))
        error('ContourRoots:ResponseEvaluation','Transfer evaluations must be finite numeric scalars (check poles, cancellations and overflow).');
    end
    v=double(v);
end
