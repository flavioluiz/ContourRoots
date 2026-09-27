function [lambda,info] = cmodes(H,region,varargin)
%CMODES Characteristic values of a square analytic matrix H(s).
%   [LAMBDA,INFO] = CMODES(H,[xmin xmax ymin ymax],'AssumeAnalytic',true)
%   searches for singularities of H without explicitly forming det(H).
%   H is a scalar-node matrix handle, a constant matrix, or a CDYN model.
%   For CDYN only its H factor is searched: these are internal modes, not
%   necessarily poles of the observed transfer. CMIMO is not accepted.
%
%   The matrix must be regular and analytic on a neighborhood of the closed
%   rectangle and nonsingular on its boundary. AssumeAnalytic is required;
%   a domain guard checks evaluations but cannot prove this assumption.
%
%   INFO reports contour counts, local algebraic counts, numerical nullity,
%   left/right vectors, residuals, unresolved boxes and bounded work. Close
%   clusters and defective repeated roots can remain unresolved even when
%   their total count is known. Numerical completeness is not certification.
%
%   Options include Derivative, DomainCheck, Scaling ('fixed' or 'none'),
%   RootTolerance, RankTolerance, SeedPoints, ContourPoints,
%   ContourRefinements, MaxDepth, MaxCells, MaxIterations, MaxEvaluations,
%   MaxMemoryMB, Display, Plot and Warn. See docs/api/cmodes.md for defaults.
%
%   Example:
%     H = @(s) [s+1 .1*exp(-s);0 s+2];
%     [r,i] = cmodes(H,[-3 0 -1 1],'AssumeAnalytic',true);
%
%   See also CDYN, CMIMO, CROOTS, CEVAL.
    [lambda,info]=matrix_modes_search(H,region,varargin{:});
    if info.options.Display
        fprintf('CMODES: %s, outer count %g, %d isolated locations.\n', ...
            info.status,info.count,numel(lambda));
        disp(table(lambda,info.localCounts,info.multiplicity,info.nullity,info.residuals, ...
            'VariableNames',{'Location','LocalCount','Multiplicity','Nullity','Residual'}));
    end
    if info.options.Plot
        figure; plot(real(lambda),imag(lambda),'x','LineWidth',1.5);
        grid on; xlim(region(1:2)); ylim(region(3:4));
        xlabel('Re(s)'); ylabel('Im(s)'); title('Matrix characteristic values (not transfer poles)');
    end
    if nargout<2 && info.options.Warn && ~info.complete
        warning('ContourRoots:MatrixUnresolved','Matrix search is unresolved; request INFO for counts and unresolved boxes.');
    end
end
