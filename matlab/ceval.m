function [G,info] = ceval(M,s,varargin)
%CEVAL Evaluate an explicit matrix model at complex spectral nodes.
%   G = CEVAL(M,S) has size ny-by-nu-by-numel(S), in S(:) order.
%   A scalar node returns a matrix; singleton dimensions follow M.Size.
%   Options: MaxEvaluations (default 1e6), MaxMemoryMB (default 256).
%   INFO counts nodes, factor calls, factorizations and block linear solves.
%   Domain guards are enforced, but do not prove analyticity.
    [G,info]=matrix_evaluate(M,s,varargin{:});
end
