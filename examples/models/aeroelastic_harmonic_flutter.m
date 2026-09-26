function result=aeroelastic_harmonic_flutter(m)
%AEROELASTIC_HARMONIC_FLUTTER Independent harmonic neutral-boundary check.
%   For undamped structure, fix real k and solve K*q = U^2*B(k)*q.
%   Hankel functions are used directly (not the Laplace Bessel evaluator).
%   Search for positive real U^2. This finite k scan is not a global proof.
    if any(m.Cs(:)~=0), error('aeroelastic:Damping','This check assumes Cs=0.'); end
    grid=logspace(-2,1,1200); values=complex(zeros(2,numel(grid)));
    for j=1:numel(grid)
        v=speedsquared(grid(j));
        if j==1, [~,ix]=sort(real(v)); v=v(ix); else
            if sum(abs(v-values(:,j-1)))>sum(abs(flipud(v)-values(:,j-1))), v=flipud(v); end
        end
        values(:,j)=v;
    end
    found=zeros(0,3);
    for branch=1:2
        for j=1:numel(grid)-1
            if imag(values(branch,j))*imag(values(branch,j+1))>=0, continue; end
            anchor=(values(branch,j)+values(branch,j+1))/2;
            kr=fzero(@(k) imaginary_branch(k,anchor),grid(j:j+1));
            v=speedsquared(kr); [~,ix]=min(abs(v-anchor)); v=v(ix);
            if real(v)>0 && abs(imag(v))<1e-7*abs(v)
                U=sqrt(real(v)); found(end+1,:)=[U kr kr*U/m.b]; %#ok<AGROW>
            end
        end
    end
    assert(~isempty(found),'No harmonic neutral candidate found.');
    found=sortrows(found,1);
    result=struct('U',found(1,1),'k',found(1,2),'omega',found(1,3),'candidates',found);
    function y=imaginary_branch(k,anchor)
        v=speedsquared(k); [~,ix]=min(abs(v-anchor)); y=imag(v(ix));
    end
    function v=speedsquared(k)
        b=m.b; a=m.a; s=1i*k/b; % evaluate air loads at U=1
        h1=besselh(1,2,k); C=h1/(h1+1i*besselh(0,2,k));
        Q=zeros(2);
        for col=1:2
            q=zeros(2,1); q(col)=1; h=q(1); alpha=q(2);
            downwash=s*h+(1+b*(.5-a)*s)*alpha;
            L=pi*m.rho*b^2*(s^2*h+s*alpha-a*b*s^2*alpha) ...
                +2*pi*m.rho*b*C*downwash;
            moment=pi*m.rho*b^2*(a*b*s^2*h-b*(.5-a)*s*alpha ...
                -b^2*(1/8+a^2)*s^2*alpha)+2*pi*m.rho*b^2*(.5+a)*C*downwash;
            Q(:,col)=[-L;moment];
        end
        v=eig(m.K,(k/b)^2*m.M+Q);
    end
end
