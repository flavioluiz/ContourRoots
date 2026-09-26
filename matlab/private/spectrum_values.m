function v = spectrum_values(f,z)
% Support vectorized and scalar-only user handles (including evalfr).
    try
        v = f(z);
        if isscalar(v), v = repmat(v,size(z)); end
        if ~isequal(size(v),size(z)), error('shape'); end
    catch
        v = arrayfun(f,z);
    end
    if ~isnumeric(v) || ~isequal(size(v),size(z))
        error('complex_spectrum:ScalarFunction','The model must return a numeric scalar for each complex scalar input.');
    end
end
