function response_plot(kind,t,y,u,i,o)
    ax=o.Parent; if isempty(ax), f=figure('Color','w'); ax=axes(f); end
    wasHeld=ishold(ax); cleaner=onCleanup(@() restore_hold(ax,wasHeld));
    if ~wasHeld, cla(ax); end
    hold(ax,'on');
    plot(ax,t,y,'LineWidth',1.5,'DisplayName','Response');
    if strcmp(kind,'lsim'), plot(ax,t,u,'--','DisplayName','Input'); legend(ax,'show'); end
    bad=~i.resolvedMask;
    if any(bad), plot(ax,t(bad),y(bad),'rx','DisplayName','Unresolved'); end
    for j=1:numel(i.singularTerms)
        term=i.singularTerms(j);
        xline(ax,term.time,':',sprintf('Dirac weight %.4g',term.weight));
    end
    xlabel(ax,'Time'); ylabel(ax,'Response'); grid(ax,'on');
    title(ax,sprintf('%s response (%s, %s)',kind,i.method,i.status));
end
function restore_hold(ax,old)
    if isgraphics(ax), if old, hold(ax,'on'); else, hold(ax,'off'); end; end
end
