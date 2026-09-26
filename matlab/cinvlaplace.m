function [y,tOut,info] = cinvlaplace(F,t,varargin)
%CINVLAPLACE Numerical inverse Laplace transform of a scalar real response.
%   [Y,T,INFO] = CINVLAPLACE(F,T,...) inverts the supplied COMPLETE F(s).
%   For a known input transform U(s), supply @(s) G(s).*U(s).
%   SingularityBound applies to F including the input; Abscissa explicitly
%   sets a valid line. Methods 'dehoog' and 'quadrature' allow nonuniform T.
%   An unknown value at zero is NaN/unresolved; use InitialValue only when
%   the analytical right limit is known. No automatic stability certificate.
    [y,tOut,info]=response_run('inverse',F,[],t,varargin,nargout==0);
end
