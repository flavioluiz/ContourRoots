function [locations,info] = complex_spectrum(model,region,varargin)
%COMPLEX_SPECTRUM Numerical zeros or poles of a scalar nonrational function.
%   [R,INFO] = COMPLEX_SPECTRUM(F,[XMIN XMAX YMIN YMAX]) finds zeros.
%   ...,'Mode','poles' finds poles of a transfer function, after cancellations.
%
%   F: symbolic scalar/symfun, function handle, continuous SISO tf/ss/zpk,
%      or struct('Numerator',N,'Denominator',D) of analytic scalar functions.
%   For a symbolic quotient, numden separates the analytic factors. Numeric
%   coefficients in a pair are interpreted as descending polynomial powers.
%
%   Options:
%     Mode                'zeros' (default) or 'poles'
%     AssumeAnalytic      true asserts that the searched factors are analytic
%                         on a neighborhood of the CLOSED rectangle (false)
%     Derivative          derivative of a single characteristic handle ([])
%     SeedPoints          optional Newton seeds, e.g. Pade poles ([])
%     RootTolerance       local isolation radius / (1+abs(s)) (1e-7)
%     ContourPoints       initial samples per edge (32)
%     ContourRefinements  maximum contour doublings (7)
%     MaxDepth            subdivision depth (24)
%     MaxCells            visited-cell budget (5000)
%     MaxIterations       Newton iterations per seed (80)
%     GridSize            exploratory handle search mesh ([25 41])
%     Singularities       known branch/accumulation points; reject a region
%                         containing any of them ([])
%     Plot, Display       false by default
%
%   INFO.complete is a NUMERICAL completeness check under the analyticity
%   contract, not a rigorous mathematical certificate. Finite sampling cannot
%   certify an arbitrary black box. Opaque handles without AssumeAnalytic use
%   exploratory Newton + local winding checks and never claim completeness.
%   For pole-mode handles, AssumeAnalytic asserts analyticity of 1/G (hence
%   no zeros of G in the rectangle). A numerator/denominator pair is preferred.
%   INFO.cancelledLocations lists every candidate where some order was
%   removed by cancellation; INFO.cancellationOrders gives the orders
%   removed, and INFO.cancellationComplete is true where the candidate was
%   removed entirely (false for a partial cancellation, which remains in
%   the result with a reduced multiplicity).
%   INFO.multiplicity counts roots in a small isolating box; roots closer
%   than this resolution may be returned as one cluster. Boundary roots are
%   not silently included: use a slightly larger or shifted rectangle.
%
%   "Exact" means evaluation of the original nonrational function without
%   Pade or a finite modal truncation; locations are floating-point numbers.

    validateattributes(region,{'numeric'},{'real','finite','vector','numel',4});
    region=double(region(:).');
    assert(region(1)<region(2) && region(3)<region(4), ...
        'complex_spectrum:Region','Use [xmin xmax ymin ymax] with increasing bounds.');
    p=inputParser;
    addParameter(p,'Mode','zeros',@(x) any(strcmpi(x,{'zeros','poles'})));
    addParameter(p,'AssumeAnalytic',false,@(x) islogical(x) && isscalar(x));
    addParameter(p,'Derivative',[],@(x) isempty(x) || isa(x,'function_handle'));
    addParameter(p,'SeedPoints',[],@(x) isnumeric(x) && all(isfinite(x(:))));
    addParameter(p,'RootTolerance',1e-7,@positive);
    addParameter(p,'ContourPoints',32,@integerPositive);
    addParameter(p,'ContourRefinements',7,@(x) integerPositive(x) && x>=2);
    addParameter(p,'MaxDepth',24,@integerPositive);
    addParameter(p,'MaxCells',5000,@integerPositive);
    addParameter(p,'MaxIterations',80,@integerPositive);
    addParameter(p,'GridSize',[25 41],@(x) isnumeric(x) && numel(x)==2 && all(x>=2 & x==round(x)));
    addParameter(p,'Singularities',[],@(x) isnumeric(x) && all(isfinite(x(:))));
    addParameter(p,'Plot',false,@(x) islogical(x) && isscalar(x));
    addParameter(p,'Display',false,@(x) islogical(x) && isscalar(x));
    parse(p,varargin{:}); opt=p.Results;
    if any(real(opt.Singularities)>=region(1) & real(opt.Singularities)<=region(2) ...
        & imag(opt.Singularities)>=region(3) & imag(opt.Singularities)<=region(4))
        error('complex_spectrum:SingularityInRegion', ...
            'Region contains a declared branch or accumulation point. Choose a region excluding it.');
    end
    [target,dt,other,analytic,source]=spectrum_model(model,opt);
    [locations,m,scan]=spectrum_solve(target,dt,region,opt,~analytic);
    cancelled=complex(zeros(0,1)); cancellationOrders=zeros(0,1); fullyCancelled=false(0,1);
    uncertain=false(size(locations));
    if ~isempty(other)
        for k=1:numel(locations)
            z=locations(k); radius=scan.locationRadius(k)/sqrt(2);
            box=[real(z)-radius real(z)+radius imag(z)-radius imag(z)+radius];
            c=spectrum_count(other,box,opt);
            if ~c.ok || c.count<0
                uncertain(k)=true;
            else
                removed=min(m(k),c.count);
                if removed>0
                    cancelled(end+1,1)=z; cancellationOrders(end+1,1)=removed; %#ok<AGROW>
                    fullyCancelled(end+1,1)=m(k)<=c.count; %#ok<AGROW>
                end
                m(k)=max(0,m(k)-c.count);
            end
        end
    end
    keep=m>0;
    locations=locations(keep);
    info=scan;
    info.complete=scan.complete && ~any(uncertain);
    info.multiplicity=m(keep);
    info.locationRadius=scan.locationRadius(keep);
    info.residuals=scan.residuals(keep);
    info.targetCountBeforeCancellations=scan.count;
    info.count=sum(info.multiplicity);
    info.cancelledLocations=cancelled;
    info.cancellationOrders=cancellationOrders;
    info.cancellationComplete=fullyCancelled;
    info.cancellationUncertain=uncertain(keep);
    info.region=region; info.mode=lower(char(opt.Mode));
    info.analyticSource=source;
    info.certified=false;
    if info.complete
        info.status='numerically_complete';
    elseif ~analytic
        info.status='exploratory';
    else
        info.status='unresolved';
    end
    info.note=['Completeness is conditional on analyticity and resolved contour sampling. ', ...
        'Multiplicity may represent an unresolved cluster within locationRadius. ', ...
        'Numerical cancellation uses this same resolution; it is not an exact symbolic identity.'];
    if opt.Display
        disp(table(locations,info.multiplicity,info.residuals,info.locationRadius, ...
            'VariableNames',{'Location','Multiplicity','TargetResidual','LocationRadius'}));
        fprintf('Status: %s; analytic source: %s\n',info.status,source);
    end
    if opt.Plot
        figure('Color','w'); ax=axes; hold(ax,'on');
        plot(real(locations),imag(locations),'x','MarkerSize',9,'LineWidth',1.8);
        if ~isempty(cancelled)
            plot(real(cancelled),imag(cancelled),'o','Color',[.6 .6 .6]);
            legend('Retained','Cancelled at search resolution','Location','best');
        end
        xline(0,':'); yline(0,':'); grid on; ax.Box='on';
        xlim(region(1:2)); ylim(region(3:4)); xlabel('Re(s)'); ylabel('Im(s)');
        title(sprintf('%s: %s',info.mode,strrep(info.status,'_',' ')));
    end
end

function yes=positive(x)
    yes=isnumeric(x) && isscalar(x) && isreal(x) && isfinite(x) && x>0;
end
function yes=integerPositive(x)
    yes=positive(x) && x==round(x);
end
