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
%   [CRITICAL,INFO] = CRITICAL_DELAYS(...) also returns INFO with fields
%     magnitudePolynomial   coefficients of |D|^2-|N|^2 in x = omega^2
%     candidateFrequencies  positive roots omega of that polynomial
%     zeroRootForAllDelays  true if s = 0 is a root for every delay
%     degenerate            true if |D(iw)| = |N(iw)| for every w (then no
%                           finite list of crossings exists; warning issued)
%     persistentFrequencies omega such that s = +-i*omega is a root for every
%                           delay (common factor of N and D; warning issued)
%   The coefficients are normalized internally, so multiplying N and D by
%   the same constant does not change the result.
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
    % The roots do not change when N and D are multiplied by the same
    % constant: normalize, so that every tolerance below is relative.
    coefScale = max([abs(N) abs(D)]);
    assert(coefScale > 0,'critical_delays:Zero','N and D must not both be zero.');
    N = N/coefScale; D = D/coefScale;
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
    % Coefficients that are round-off relative to the largest one are zero.
    % An exactly zero constant term then gives an exact root x = 0, which is
    % discarded; every remaining positive root is a genuine frequency, however
    % small (no absolute threshold on omega).
    gScale = norm(gAsc,inf);
    termScale = max([abs(dProd) abs(nProd)]);
    degenerate = gScale <= 1e3*eps*termScale;
    if degenerate
        gAsc(:) = 0;
    else
        gAsc(abs(gAsc) < 100*eps(gScale)) = 0;
    end
    gX = trim_poly(fliplr(gAsc));

    if isscalar(gX)
        xRoots = [];
    else
        allX = roots(gX);
        keep = abs(imag(allX)) <= opt.RootTolerance*abs(allX) & real(allX) > 0;
        xRoots = sort(real(allX(keep)));
        xRoots = unique_relative(xRoots,100*opt.RootTolerance);
    end
    frequencies = sqrt(xRoots(:));
    if degenerate
        warning('critical_delays:Degenerate', ['|D(i*w)| = |N(i*w)| for every w: ' ...
            'imaginary-axis roots are not isolated in (w,T), so no finite list of ' ...
            'critical delays exists. The result is empty and info.degenerate is true.']);
    end
    persistent_ = zeros(0,1);

    rows = zeros(0,6);
    directions = strings(0,1);
    for iw = 1:numel(frequencies)
        omega = frequencies(iw);
        Nw = polyval(N,1i*omega);
        Dw = polyval(D,1i*omega);
        % Size of the terms of N and D at this frequency, for relative tests.
        termN = polyval(abs(N),omega); termD = polyval(abs(D),omega);
        if abs(Nw) <= 1e3*eps*termN && abs(Dw) <= 1e3*eps*termD
            % Common factor on the imaginary axis: i*omega is a root for
            % EVERY delay. It is not a crossing; report it separately.
            persistent_(end+1,1) = omega; %#ok<AGROW>
            continue
        end
        ratio = -Dw/Nw;
        basePhase = mod(-angle(ratio),2*pi);
        kFirst = max(0,ceil((opt.Tmin*omega-basePhase)/(2*pi)-opt.RootTolerance));
        kLast = floor((Tmax*omega-basePhase)/(2*pi)+opt.RootTolerance);
        for branch = kFirst:kLast
            delay = (basePhase+2*pi*branch)/omega;
            tolT = opt.RootTolerance*(1+Tmax);
            if delay < opt.Tmin-tolT || delay > Tmax+tolT
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
            termScale = polyval(abs(D),omegaRefined)+polyval(abs(N),omegaRefined);
            if residual <= opt.ResidualTolerance*termScale && ...
                    delayRefined >= opt.Tmin-tolT && delayRefined <= Tmax+tolT
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

    % Common factors on the imaginary axis can also be multiple roots of the
    % magnitude polynomial, which roots() may return slightly complex; detect
    % them directly from the roots of D.
    rD = roots(D);
    onAxis = rD(abs(real(rD)) <= 1e-6*max(1,abs(rD)) & imag(rD) > 0);
    for r = onAxis(:).'
        omega = imag(r);
        if abs(polyval(N,1i*omega)) <= 1e-6*polyval(abs(N),omega) && ...
                all(abs(persistent_-omega) > 1e-6*omega)
            persistent_(end+1,1) = omega; %#ok<AGROW>
        end
    end
    D0 = polyval(D,0); N0 = polyval(N,0);
    zeroRoot = abs(D0+N0) <= 1e3*eps*(abs(D0)+abs(N0));
    if ~isempty(persistent_)
        warning('critical_delays:PersistentImaginaryRoots', ['s = +-%si is a root for ' ...
            'every delay (common factor of N and D on the imaginary axis). It is not ' ...
            'a crossing; see info.persistentFrequencies.'], num2str(persistent_.',6));
    end
    info = struct('magnitudePolynomial',gX,'candidateFrequencies',frequencies, ...
        'zeroRootForAllDelays',zeroRoot,'delayRange',[opt.Tmin Tmax], ...
        'degenerate',degenerate,'persistentFrequencies',persistent_);
end

function [omega,T] = refine_crossing(N,D,omega,T)
% Newton on [Re F; Im F] = 0 in the scaled unknowns (omega/omega0, T/T0),
% so that very small frequencies and very large delays are treated alike.
    dN = polyder(N); dD = polyder(D);
    so = omega; sT = max(T,1/omega);
    for iteration = 1:12
        s = 1i*omega;
        e = exp(-s*T);
        f = polyval(D,s)+e*polyval(N,s);
        dFs = polyval(dD,s)+e*(polyval(dN,s)-T*polyval(N,s));
        dFw = 1i*dFs*so;
        dFT = -s*e*polyval(N,s)*sT;
        J = [real(dFw) real(dFT); imag(dFw) imag(dFT)];
        if rcond(J) < 1e-13, break; end
        delta = J\[real(f);imag(f)];
        omega = omega-delta(1)*so;
        T = T-delta(2)*sT;
        if norm(delta) <= 1e-14, break; end
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
