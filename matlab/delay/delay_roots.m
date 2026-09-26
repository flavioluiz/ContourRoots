function [roots_delay, info] = delay_roots(N,D,T,xlimits,ylimits,nx,ny,varargin)
%DELAY_ROOTS Roots of D(s) + exp(-s*T) N(s) in a rectangle.
%   R = DELAY_ROOTS(N,D,T,XLIM,YLIM,NX,NY) launches damped complex
%   Newton iterations from a rectangular mesh. Poles of several diagonal
%   Pade approximations are also used as initial guesses.
%
%   [R,INFO] = DELAY_ROOTS(...,'Name',Value) accepts:
%     'SeedPoints'       extra complex initial guesses (default [])
%     'PadeSeedOrders'   Pade orders used as seeds (default [2 4 8 12])
%     'MaxIterations'    Newton iterations per seed (default 80)
%     'MaxRefinements'   mesh refinement passes (default 2)
%     'NewtonTolerance'  relative step tolerance (default 1e-11)
%     'ResidualTolerance' residual tolerance (default 1e-8)
%     'ClusterTolerance' relative clustering tolerance (default 2e-6)
%     'VerifyCount'      check the number of roots with the argument
%                        principle; INFO.countMatches is false when the
%                        count is inconclusive (NaN)
%                        (default true)
%     'Display'          print a root table (default false)
%     'Plot'             draw a pole map (default false)
%
%   INFO reports residuals, the number of Newton seeds and the root count
%   estimated independently with the argument principle.

    if nargin < 6 || isempty(nx), nx = 50; end
    if nargin < 7 || isempty(ny), ny = 80; end

    validateattributes(N,{'numeric'},{'vector','nonempty','finite'});
    validateattributes(D,{'numeric'},{'vector','nonempty','finite'});
    validateattributes(T,{'numeric'},{'scalar','real','finite','nonnegative'});
    validateattributes(xlimits,{'numeric'},{'vector','numel',2,'increasing'});
    validateattributes(ylimits,{'numeric'},{'vector','numel',2,'increasing'});
    validateattributes(nx,{'numeric'},{'scalar','integer','>=',2});
    validateattributes(ny,{'numeric'},{'scalar','integer','>=',2});

    p = inputParser;
    addParameter(p,'SeedPoints',[],@isnumeric);
    addParameter(p,'PadeSeedOrders',[2 4 8 12], ...
        @(x) isnumeric(x) && isvector(x) && all(x >= 1));
    addParameter(p,'MaxIterations',80,@(x) isscalar(x) && x >= 1);
    addParameter(p,'MaxRefinements',2,@(x) isscalar(x) && x >= 1);
    addParameter(p,'NewtonTolerance',1e-11,@(x) isscalar(x) && x > 0);
    addParameter(p,'ResidualTolerance',1e-8,@(x) isscalar(x) && x > 0);
    addParameter(p,'ClusterTolerance',2e-6,@(x) isscalar(x) && x > 0);
    addParameter(p,'VerifyCount',true,@(x) islogical(x) || ismember(x,[0 1]));
    addParameter(p,'Display',false,@(x) islogical(x) || ismember(x,[0 1]));
    addParameter(p,'Plot',false,@(x) islogical(x) || ismember(x,[0 1]));
    parse(p,varargin{:});
    opt = p.Results;

    N = trim_leading_zeros(N(:).');
    D = trim_leading_zeros(D(:).');
    dN = polyder(N);
    dD = polyder(D);
    F = @(s) polyval(D,s) + exp(-s*T).*polyval(N,s);
    dF = @(s) polyval(dD,s) + exp(-s*T).* ...
        (polyval(dN,s) - T*polyval(N,s));

    expected = NaN;
    winding = NaN;
    if opt.VerifyCount
        [expected,winding] = delay_root_count(N,D,T,xlimits,ylimits);
    end

    seedPoints = opt.SeedPoints(:);
    for q = unique(round(opt.PadeSeedOrders(:).'))
        rp = roots(pade_characteristic(N,D,T,q));
        keep = in_rectangle(rp,xlimits,ylimits,0.25);
        seedPoints = [seedPoints; rp(keep)]; %#ok<AGROW>
    end

    candidates = complex(zeros(0,1));
    attempts = 0;
    convergedCount = 0;
    refinementsUsed = 0;
    regionScale = max(diff(xlimits),diff(ylimits));
    escapeRadius = 6*(1 + max(abs([xlimits ylimits])));

    for level = 1:round(opt.MaxRefinements)
        refinementsUsed = level;
        nxl = ceil(nx*1.45^(level-1));
        nyl = ceil(ny*1.45^(level-1));
        [X,Y] = ndgrid(linspace(xlimits(1),xlimits(2),nxl), ...
                       linspace(ylimits(1),ylimits(2),nyl));
        starts = [X(:)+1i*Y(:); seedPoints];
        attempts = attempts + numel(starts);

        for j = 1:numel(starts)
            s = starts(j);
            converged = false;
            for k = 1:round(opt.MaxIterations)
                fs = F(s);
                dfs = dF(s);
                if ~isfinite_complex(fs) || ~isfinite_complex(dfs) || ...
                        abs(dfs) < 100*eps*(1+abs(fs))
                    break
                end
                step = fs/dfs;
                maxStep = max(1,0.35*regionScale);
                if abs(step) > maxStep
                    step = step*(maxStep/abs(step));
                end

                % Backtracking prevents Newton from crossing large
                % exponential-overflow regions unnecessarily.
                alpha = 1;
                snew = s - step;
                oldResidual = abs(fs);
                while alpha > 1/32
                    fnew = F(snew);
                    if isfinite_complex(fnew) && ...
                            (abs(fnew) <= 1.05*oldResidual || abs(step) < 1e-7)
                        break
                    end
                    alpha = alpha/2;
                    snew = s - alpha*step;
                end
                if ~isfinite_complex(snew) || abs(snew) > escapeRadius
                    break
                end
                if abs(snew-s) <= opt.NewtonTolerance*(1+abs(snew))
                    s = snew;
                    converged = abs(F(s)) <= opt.ResidualTolerance*(1+abs(polyval(D,s)));
                    break
                end
                s = snew;
            end

            if converged && in_rectangle(s,xlimits,ylimits,0)
                candidates(end+1,1) = s; %#ok<AGROW>
                convergedCount = convergedCount + 1;
            end
        end

        roots_delay = cluster_roots(candidates,opt.ClusterTolerance,F);
        if ~isnan(expected) && numel(roots_delay) >= expected
            break
        end
    end

    roots_delay = cluster_roots(candidates,opt.ClusterTolerance,F);
    [~,idx] = sortrows([-real(roots_delay), abs(imag(roots_delay))],[1 2]);
    roots_delay = roots_delay(idx);
    residuals = abs(F(roots_delay));

    info = struct('residuals',residuals,'seedCount',attempts, ...
        'convergedSeedCount',convergedCount,'refinementsUsed',refinementsUsed, ...
        'argumentPrincipleCount',expected,'windingNumber',winding, ...
        'countMatches',~opt.VerifyCount || (~isnan(expected) && numel(roots_delay)==expected), ...
        'characteristic',F);

    if opt.Display
        fprintf('\nRoots found: %d',numel(roots_delay));
        if ~isnan(expected), fprintf(' (argument-principle count: %d)',expected); end
        fprintf('\n\n       Re(s)              Im(s)             |F(s)|\n');
        fprintf('-----------------------------------------------------------\n');
        for k = 1:numel(roots_delay)
            fprintf('%14.8f   %14.8f   %12.3e\n',real(roots_delay(k)), ...
                imag(roots_delay(k)),residuals(k));
        end
    end

    if opt.Plot
        figure('Color','w');
        plot(real(roots_delay),imag(roots_delay),'x','MarkerSize',8,'LineWidth',1.6);
        grid on; xline(0,'--'); yline(0,'--');
        xlabel('Real(s)'); ylabel('Imag(s)');
        title(sprintf('D(s) + e^{-sT}N(s) = 0, T = %.4g',T));
        xlim(xlimits); ylim(ylimits);
    end
end

function tf = isfinite_complex(z)
    tf = isfinite(real(z)) && isfinite(imag(z));
end

function tf = in_rectangle(z,xl,yl,paddingFraction)
    px = paddingFraction*diff(xl);
    py = paddingFraction*diff(yl);
    tf = real(z) >= xl(1)-px & real(z) <= xl(2)+px & ...
         imag(z) >= yl(1)-py & imag(z) <= yl(2)+py;
end

function c = cluster_roots(values,tol,F)
    c = complex(zeros(0,1));
    if isempty(values), return; end
    [~,order] = sort(abs(F(values)),'ascend');
    values = values(order);
    for j = 1:numel(values)
        r = values(j);
        if isempty(c) || all(abs(r-c) > tol*(1+max(abs(r),abs(c))))
            c(end+1,1) = r; %#ok<AGROW>
        end
    end
end

function p = trim_leading_zeros(p)
    first = find(p~=0,1,'first');
    if isempty(first), p = 0; else, p = p(first:end); end
end
