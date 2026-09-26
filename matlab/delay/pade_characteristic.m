function characteristic = pade_characteristic(N,D,T,n)
%PADE_CHARACTERISTIC Polynomial characteristic equation using [n/n] Pade.
%   D(s) Q_n(s) + N(s) P_n(s) = 0, where P_n/Q_n approximates exp(-sT).

    [p,q] = pade_delay(T,n);
    a = conv(D,q);
    b = conv(N,p);
    L = max(numel(a),numel(b));
    characteristic = [zeros(1,L-numel(a)) a] + [zeros(1,L-numel(b)) b];
    first = find(abs(characteristic) > eps(max(1,norm(characteristic,inf))),1);
    if isempty(first)
        characteristic = 0;
    else
        characteristic = characteristic(first:end);
    end
end
