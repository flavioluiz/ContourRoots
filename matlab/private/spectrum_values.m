function v = spectrum_values(f,z)
% Evaluate a user function on an array of points.
%   Vectorized handles are called once on the whole array. A handle that
%   errors, or that returns a result of a different size (for example a
%   scalar-only function such as @(s) prod([s-1, s+1]), or a constant), is
%   evaluated point by point. A scalar result is never replicated: that
%   would silently treat a non-vectorized function as a constant.
    vectorized = false;
    if numel(z) > 1
        try
            v = f(z);
            vectorized = isnumeric(v) && isequal(size(v),size(z));
        catch
            vectorized = false;
        end
    end
    if ~vectorized
        v = zeros(size(z));
        for k = 1:numel(z)
            value = f(z(k));
            if ~(isnumeric(value) && isscalar(value))
                error('complex_spectrum:ScalarFunction', ...
                    'The model must return a numeric scalar for each complex scalar input.');
            end
            v(k) = value;
        end
    end
end
