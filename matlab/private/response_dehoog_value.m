function v=response_dehoog_value(fp,t,P,sigma,M)
% Q-D recurrence shared by scalar and matrix preparation (see notices).
    if all(fp==0), v=0; return; end
    e=complex(zeros(2*M+1,M+1)); q=complex(zeros(2*M,M));
    q(1,1)=fp(2)/(fp(1)/2); q(2:2*M,1)=fp(3:end)./fp(2:end-1);
    for r=1:M
        n=2*(M-r)+1;
        e(1:n,r+1)=q(2:n+1,r)-q(1:n,r)+e(2:n+1,r);
        if r<M, q(1:n-1,r+1)=q(2:n,r).*e(2:n,r+1)./e(1:n-1,r+1); end
    end
    d=complex(zeros(2*M+1,1)); d(1)=fp(1)/2;
    for r=1:M, d(2*r)=-q(1,r); d(2*r+1)=-e(1,r+1); end
    z=exp(2i*pi*t/P); ap=0; ac=d(1); bp=1; bc=1;
    for j=2:2*M
        an=ac+d(j)*ap*z; bn=bc+d(j)*bp*z;
        ap=ac; ac=an; bp=bc; bc=bn;
        scale=max([abs(ac),abs(ap),abs(bc),abs(bp)]);
        if scale>1e100, ac=ac/scale; ap=ap/scale; bc=bc/scale; bp=bp/scale; end
    end
    brem=(1+(d(2*M)-d(2*M+1))*z)/2;
    x=d(2*M+1)*z/brem;
    remainder=brem*x/(sqrt(1+x)+1); % stable sqrt(1+x)-1
    v=exp(sigma*t)*2/P*real((ac+remainder*ap)/(bc+remainder*bp));
end
