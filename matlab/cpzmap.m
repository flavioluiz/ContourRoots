function [p,z,infoP,infoZ] = cpzmap(G,region,varargin)
%CPZMAP Pole-zero map of a nonrational transfer function.
%   CPZMAP(G,REGION) computes the poles (with CPOLES) and zeros (with
%   CZEROS) of G inside REGION = [xmin xmax ymin ymax] and plots them in
%   the complex plane: poles as x, zeros as o, as in PZMAP. Locations
%   removed by pole-zero cancellation are shown as small grey dots. The
%   search rectangle is drawn as a dashed box; roots outside it are not
%   searched and not shown.
%
%   [P,Z] = CPZMAP(G,REGION) returns the poles and zeros without plotting.
%   [P,Z,INFOP,INFOZ] = CPZMAP(...) also returns the diagnostics of each
%   search (see CROOTS and CPOLES).
%
%   CPZMAP(...,Name,Value) passes options to both searches, for example
%   'AssumeAnalytic',true. Use CPZMAP(...,'Plot',true) together with
%   output arguments to both plot and return results. CPZMAP(AX,...)
%   plots into the axes AX.
%
%   Example
%       G = ndpair(@(s) sinh(s/2), @(s) sinh(s));
%       cpzmap(G,[-1 1 -10 10],'AssumeAnalytic',true)
%
%   See also CPOLES, CZEROS, CROOTS, NDPAIR, PZMAP.

    ax = [];
    if nargin >= 1 && isscalar(G) && isgraphics(G,'axes')
        ax = G;
        if nargin < 3
            error('ContourRoots:cpzmap','Use CPZMAP(AX,G,REGION,...).');
        end
        G = region; region = varargin{1}; varargin(1) = [];
    end
    if nargin < 2
        error('ContourRoots:Region','CPZMAP needs a search rectangle [xmin xmax ymin ymax].');
    end
    k = find(strcmpi(varargin(1:2:end),'Plot'),1,'last');
    doPlot = nargout == 0;
    if ~isempty(k)
        doPlot = logical(varargin{2*k}); varargin(2*k-1:2*k) = [];
    end
    wantInfo = nargout >= 3;
    [p,infoP] = cr_search('cpzmap','poles',G,region,varargin,~wantInfo);
    [z,infoZ] = cr_search('cpzmap','zeros',G,region,varargin,~wantInfo);

    if doPlot
        if isempty(ax), ax = newplot; end
        holdState = ishold(ax); hold(ax,'on');
        restore = onCleanup(@() set_hold(ax,holdState));
        cancelled = [infoP.cancelledLocations; infoZ.cancelledLocations];
        frame = region([1 2 2 1 1]) + 1i*region([3 3 4 4 3]);
        plot(ax,real(frame),imag(frame),'--','Color',[.6 .6 .6],'HandleVisibility','off');
        xline(ax,0,':','Color',[.4 .4 .4],'HandleVisibility','off');
        yline(ax,0,':','Color',[.4 .4 .4],'HandleVisibility','off');
        hp = plot(ax,real(p),imag(p),'x','MarkerSize',10,'LineWidth',1.8, ...
            'Color',[0 0.4470 0.7410],'DisplayName','Poles');
        hz = plot(ax,real(z),imag(z),'o','MarkerSize',8,'LineWidth',1.5, ...
            'Color',[0.8500 0.3250 0.0980],'DisplayName','Zeros');
        handles = [hp hz];
        if ~isempty(cancelled)
            hc = plot(ax,real(cancelled),imag(cancelled),'.','MarkerSize',12, ...
                'Color',[.6 .6 .6],'DisplayName','Cancelled');
            handles(end+1) = hc;
        end
        grid(ax,'on'); box(ax,'on');
        pad = 0.04*[diff(region(1:2)) diff(region(3:4))];
        xlim(ax,region(1:2)+[-1 1]*pad(1)); ylim(ax,region(3:4)+[-1 1]*pad(2));
        xlabel(ax,'Real axis  Re(s)'); ylabel(ax,'Imaginary axis  Im(s)');
        status = infoP.status;
        if ~strcmp(infoZ.status,status), status = [status ' / ' infoZ.status]; end
        title(ax,sprintf('Pole-zero map (%s)',strrep(status,'_',' ')));
        legend(ax,handles,'Location','best');
    end
    if nargout == 0
        clear p z infoP infoZ
    end
end

function set_hold(ax,state)
    if isvalid(ax) && ~state, hold(ax,'off'); end
end
