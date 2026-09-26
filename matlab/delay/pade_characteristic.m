function characteristic = pade_characteristic(N,D,T,n)
%PADE_CHARACTERISTIC Polynomial characteristic equation using [n/n] Pade.
%   C = PADE_CHARACTERISTIC(N,D,T,N) returns the coefficients (descending
%   powers) of D(s) Q_n(sT) + N(s) P_n(sT), where P_n/Q_n is the diagonal
%   Pade approximation of exp(-sT) (see PADE_DELAY). For deg N < deg D the
%   degree is deg D + n whenever T > 0.
%
%   Only exactly zero leading coefficients are removed: small coefficients
%   are legitimate (for small T and large n they can be many orders of
%   magnitude below the largest one). For large n the expanded polynomial
%   is badly conditioned; ROOTS of it may be inaccurate, and a balanced
%   state-space realization is preferable.

    [p,q] = pade_delay(T,n);
    a = conv(D,q);
    b = conv(N,p);
    L = max(numel(a),numel(b));
    characteristic = [zeros(1,L-numel(a)) a] + [zeros(1,L-numel(b)) b];
    first = find(characteristic ~= 0,1);
    if isempty(first)
        characteristic = 0;
    else
        characteristic = characteristic(first:end);
    end
end
