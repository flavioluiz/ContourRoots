function [critical,info] = critical_delays(N,D,Tmax,varargin)
%CRITICAL_DELAYS Exact imaginary-axis crossings of D(s)+N(s)exp(-sT)=0.
%   CRITICAL = CRITICAL_DELAYS(N,D,TMAX) returns all positive frequencies
%   and delays 0 <= T <= TMAX for which s = +/-1i*omega is a root.
%
%   For real polynomial coefficients, an imaginary-axis crossing must obey
%       |D(i*omega)|^2 - |N(i*omega)|^2 = 0.
%   This condition is converted into a polynomial in x = omega^2. The phase
%   equation then generates every delay branch analytically. A two-variable
%   Newton correction is applied to (omega,T) before the result is returned.
%
%   Name-value options:
%     'Tmin'               lower delay bound (default 0)
%     'RootTolerance'      algebraic filtering tolerance (default 1e-9)
%     'ResidualTolerance'  accepted |F(i*omega,T)| (default 1e-9)
%
%   Output table columns:
%     Frequency, Delay, Branch, CrossingSpeed, Direction, Residual.
%   CrossingSpeed is Re(ds/dT). A positive value means a left-to-right
%   (destabilizing) crossing; a negative value means stabilizing.

    validateattributes(N,{'numeric'},{'vector','nonempty','finite','real'});
    validateattributes(D,{'numeric'},{'vector','nonempty','finite','real'});
    validateattributes(Tmax,{'numeric'},{'scalar','real','finite','nonnegative'});
    p = inputParser;
    addParameter(p,'Tmin',0,@(x) isscalar(x) && isreal(x) && x >= 0);
    addParameter(p,'RootTolerance',1e-9,@(x) isscalar(x) && x > 0);
    addParameter(p,'ResidualTolerance',1e-9,@(x) isscalar(x) && x > 0);
    parse(p,varargin{:});
    opt = p.Results;
    assert(opt.Tmin <= Tmax,'Tmin must not exceed Tmax.');

    N = trim_poly(N(:).');
    D = trim_poly(D(:).');
    maxDegree = max(numel(N),numel(D))-1;
    L = 2*maxDegree;
    dProd = pad_left(conv(D,poly_minus_argument(D)),L+1);
    nProd = pad_left(conv(N,poly_minus_argument(N)),L+1);
    gS = dProd-nProd;

    % gS(i*omega) is an even real polynomial. Set x=omega^2.
    gAsc = zeros(1,maxDegree+1);
    for k = 0:maxDegree
        power = 2*k;
        gAsc(k+1) = real(gS(L-power+1))*(-1)^k;
    end
    scale = max(1,norm(gAsc,inf));
    gAsc(abs(gAsc) < 100*eps(scale)) = 0;
    gX = trim_poly(fliplr(gAsc));

    if isscalar(gX)
        xRoots = [];
    else
        allX = roots(gX);
        keep = abs(imag(allX)) <= opt.RootTolerance*(1+abs(real(allX))) & ...
               real(allX) > opt.RootTolerance;
        xRoots = sort(real(allX(keep)));
        xRoots = unique_relative(xRoots,100*opt.RootTolerance);
    end
    frequencies = sqrt(xRoots(:));

    rows = zeros(0,6);
    directions = strings(0,1);
    for iw = 1:numel(frequencies)
        omega = frequencies(iw);
        Nw = polyval(N,1i*omega);
        Dw = polyval(D,1i*omega);
        if abs(Nw) < opt.RootTolerance*(1+norm(N,1))
            continue
        end
        ratio = -Dw/Nw;
        basePhase = mod(-angle(ratio),2*pi);
        kFirst = max(0,ceil((opt.Tmin*omega-basePhase)/(2*pi)-opt.RootTolerance));
        kLast = floor((Tmax*omega-basePhase)/(2*pi)+opt.RootTolerance);
        for branch = kFirst:kLast
            delay = (basePhase+2*pi*branch)/omega;
            if delay < opt.Tmin-opt.RootTolerance || delay > Tmax+opt.RootTolerance
                continue
            end
            [omegaRefined,delayRefined] = refine_crossing(N,D,omega,delay);
            s = 1i*omegaRefined;
            e = exp(-s*delayRefined);
            f = polyval(D,s)+e*polyval(N,s);
            dFs = polyval(polyder(D),s)+e*(polyval(polyder(N),s) ...
                - delayRefined*polyval(N,s));
            dsdt = s*e*polyval(N,s)/dFs;
            speed = real(dsdt);
            if speed > 100*eps*(1+abs(dsdt))
                label = "destabilizing";
            elseif speed < -100*eps*(1+abs(dsdt))
                label = "stabilizing";
            else
                label = "tangent/degenerate";
            end
            residual = abs(f);
            if residual <= opt.ResidualTolerance*(1+abs(polyval(D,s))) && ...
                    delayRefined >= opt.Tmin-opt.RootTolerance && ...
                    delayRefined <= Tmax+opt.RootTolerance
                rows(end+1,:) = [omegaRefined,delayRefined,branch,speed, ...
                    real(dsdt),residual]; %#ok<AGROW>
                directions(end+1,1) = label; %#ok<AGROW>
            end
        end
    end

    if isempty(rows)
        critical = table(zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
            strings(0,1),zeros(0,1),'VariableNames', ...
            {'Frequency','Delay','Branch','CrossingSpeed','Direction','Residual'});
    else
        [~,order] = sortrows(rows(:,[2 1]),[1 2]);
        rows = rows(order,:);
        directions = directions(order);
        critical = table(rows(:,1),rows(:,2),round(rows(:,3)),rows(:,4), ...
            directions,rows(:,6),'VariableNames', ...
            {'Frequency','Delay','Branch','CrossingSpeed','Direction','Residual'});
    end

    zeroRoot = abs(polyval(D,0)+polyval(N,0)) <= ...
        opt.RootTolerance*(1+abs(polyval(D,0))+abs(polyval(N,0)));
    info = struct('magnitudePolynomial',gX,'candidateFrequencies',frequencies, ...
        'zeroRootForAllDelays',zeroRoot,'delayRange',[opt.Tmin Tmax]);
end

function [omega,T] = refine_crossing(N,D,omega,T)
    dN = polyder(N); dD = polyder(D);
    for iteration = 1:12
        s = 1i*omega;
        e = exp(-s*T);
        f = polyval(D,s)+e*polyval(N,s);
        dFs = polyval(dD,s)+e*(polyval(dN,s)-T*polyval(N,s));
        dFw = 1i*dFs;
        dFT = -s*e*polyval(N,s);
        J = [real(dFw) real(dFT); imag(dFw) imag(dFT)];
        if rcond(J) < 1e-13, break; end
        delta = J\[real(f);imag(f)];
        omega = omega-delta(1);
        T = T-delta(2);
        if norm(delta) <= 1e-13*(1+norm([omega;T])), break; end
    end
end

function pm = poly_minus_argument(p)
    degree = numel(p)-1;
    pm = p.*((-1).^(degree:-1:0));
end

function p = pad_left(p,L)
    p = [zeros(1,L-numel(p)) p];
end

function p = trim_poly(p)
    scale = max(1,norm(p,inf));
    first = find(abs(p)>100*eps(scale),1,'first');
    if isempty(first), p=0; else, p=p(first:end); end
end

function y = unique_relative(x,tol)
    y = zeros(0,1);
    for k=1:numel(x)
        if isempty(y) || all(abs(x(k)-y)>tol*(1+abs(x(k))))
            y(end+1,1)=x(k); %#ok<AGROW>
        end
    end
end
