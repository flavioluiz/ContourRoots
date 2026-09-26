function [numDelay,denDelay] = pade_delay(T,n)
%PADE_DELAY Diagonal [n/n] Pade approximation of exp(-s*T).
%   Coefficients are returned in descending powers of s and this function
%   does not require Control System Toolbox.

    validateattributes(T,{'numeric'},{'scalar','real','finite','nonnegative'});
    validateattributes(n,{'numeric'},{'scalar','integer','nonnegative'});
    if n == 0 || T == 0
        numDelay = 1;
        denDelay = 1;
        return
    end
    k = 0:n;
    % log-gamma avoids intermediate factorial overflow for useful orders.
    logc = gammaln(2*n-k+1) + gammaln(n+1) ...
         - gammaln(2*n+1) - gammaln(k+1) - gammaln(n-k+1);
    asc = exp(logc).*T.^k;
    denDelay = fliplr(asc);
    numDelay = fliplr(asc.*(-1).^k);
    scale = denDelay(1);
    numDelay = numDelay/scale;
    denDelay = denDelay/scale;
end
