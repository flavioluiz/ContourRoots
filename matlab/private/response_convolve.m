function y=response_convolve(kernel,weights,t,sigma,o)
% Weighted linear convolution of integrated kernels and basis coefficients.
% No dt factor: weights are jumps/slopes, not sampled input times g(t).
    n=numel(t); L=2^nextpow2(2*n-1);
    if L>o.MaxPoints || 128*L>o.MaxMemoryMB*2^20
        error('ContourRoots:ResponseBudget','Linear convolution exceeds MaxPoints/MaxMemoryMB.');
    end
    if ~isfinite(sigma), sigma=0; end
    if abs(sigma)*max(t)>log(1e12)
        error('ContourRoots:ResponseBudget','Convolution weighting exceeds 1e12; shorten the time interval.');
    end
    w=exp(-sigma*t);
    z=ifft(fft(kernel.*w,L).*fft(weights.*w,L));
    y=real(z(1:n))./w;
end
