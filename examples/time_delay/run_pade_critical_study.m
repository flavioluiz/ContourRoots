function study = run_pade_critical_study(rootDir)
%RUN_PADE_CRITICAL_STUDY Compare exact and Pade critical delays.
    % Resolve paths from this file, so the study runs from any folder.
    here=fileparts(mfilename('fullpath')); repo=fileparts(fileparts(here));
    if exist('croots','file')~=2, run(fullfile(repo,'setup_contourroots.m')); end
    if nargin<1, rootDir=fullfile(repo,'output','time_delay'); end
    figDir=fullfile(rootDir,'figures');
    dataDir=fullfile(rootDir,'data');
    if ~exist(figDir,'dir'), mkdir(figDir); end
    if ~exist(dataDir,'dir'), mkdir(dataDir); end

    cases(1)=struct('id','P1','N',2,'D',[1 1],'name','P1: s+1+2e^{-sT}');
    cases(2)=struct('id','P3','N',3,'D',[1 .8 4],'name','P3: s^2+0.8s+4+3e^{-sT}');
    orders=[1 2 3 4 5 6 8 10 12];
    shown=[1 2 4 8 12];
    Tmax=10;
    exact=cell(2,1); predicted=cell(2,1);
    firstError=NaN(2,numel(orders)); coverage=NaN(2,numel(orders));

    for ic=1:2
        c=cases(ic);
        exact{ic}=critical_delays(c.N,c.D,Tmax);
        predicted{ic}=pade_critical_delays(c.N,c.D,Tmax,orders);
        writetable(exact{ic},fullfile(dataDir,[lower(c.id) '_exact_critical_to_T10.csv']));
        writetable(predicted{ic},fullfile(dataDir,[lower(c.id) '_pade_critical_to_T10.csv']));
        firstExact=exact{ic}.Delay(1);
        firstFrequency=exact{ic}.Frequency(1);
        for io=1:numel(orders)
            rows=predicted{ic}(predicted{ic}.Order==orders(io) & ...
                abs(predicted{ic}.Frequency-firstFrequency)<1e-7,:);
            if ~isempty(rows)
                firstError(ic,io)=abs(rows.Delay(1)-firstExact)/firstExact;
            end
            coverage(ic,io)=height(predicted{ic}(predicted{ic}.Order==orders(io),:))/height(exact{ic});
        end
    end

    % Phase mechanism for P1.
    f1=figure('Color','w','Position',[100 100 1000 520]); ax=axes(f1); hold(ax,'on');
    omega=exact{1}.Frequency(1); T=linspace(0,Tmax,2500);
    plot(ax,T,-omega*T,'k-','LineWidth',2.3,'DisplayName','Exact phase -\omegaT');
    colors=lines(4);
    for j=1:4
        n=[1 2 4 8]; n=n(j);
        [pn,qn]=pade_delay(1,n);
        phase=unwrap(angle(polyval(pn,1i*omega*T)./polyval(qn,1i*omega*T)));
        plot(ax,T,phase,'LineWidth',1.7,'Color',colors(j,:), ...
            'DisplayName',sprintf('Padé %d/%d',n,n));
    end
    target0=angle(-polyval(cases(1).D,1i*omega)/polyval(cases(1).N,1i*omega));
    if target0>0, target0=target0-2*pi; end
    for k=0:2
        y=target0-2*pi*k;
        yline(ax,y,':','Color',[.48 .48 .48],'HandleVisibility','off');
        text(ax,Tmax*.985,y+0.18,sprintf('branch k=%d',k), ...
            'HorizontalAlignment','right','Color',[.35 .35 .35]);
    end
    grid(ax,'on'); box(ax,'on'); xlim(ax,[0 Tmax]); ylim(ax,[-18.5 .5]);
    xlabel(ax,'Delay T'); ylabel(ax,'Unwrapped phase [rad]');
    title(ax,'Why Padé shifts or misses critical delays (P1)');
    legend(ax,'Location','southwest');
    save_figure(f1,figDir,'pade_critical_phase');

    % Accuracy and coverage.
    f2=figure('Color','w','Position',[120 120 1080 430]);
    tl=tiledlayout(f2,1,2,'TileSpacing','compact','Padding','compact');
    ax=nexttile(tl); hold(ax,'on');
    for ic=1:2
        semilogy(ax,orders,firstError(ic,:),'-o','LineWidth',1.9, ...
            'MarkerFaceColor','auto','DisplayName',cases(ic).id);
    end
    grid(ax,'on'); box(ax,'on'); set(ax,'YScale','log'); ylim(ax,[1e-15 1]);
    xlabel(ax,'Diagonal order n'); ylabel(ax,'Relative error in the first T_c');
    title(ax,'Accuracy of the first delay margin'); legend(ax,'Location','southwest');
    ax=nexttile(tl); hold(ax,'on');
    for ic=1:2
        plot(ax,orders,100*coverage(ic,:),'-s','LineWidth',1.9, ...
            'MarkerFaceColor','auto','DisplayName',sprintf('%s (%d exact events)', ...
            cases(ic).id,height(exact{ic})));
    end
    grid(ax,'on'); box(ax,'on'); ylim(ax,[0 108]); yticks(ax,0:20:100);
    xlabel(ax,'Diagonal order n'); ylabel(ax,'Events predicted up to T=10 [%]');
    title(ax,'Coverage of crossing branches'); legend(ax,'Location','southeast');
    title(tl,'Local accuracy versus global coverage');
    save_figure(f2,figDir,'pade_critical_accuracy');

    % Event maps show convergence and missing branches directly.
    f3=figure('Color','w','Position',[140 140 1080 450]);
    tl3=tiledlayout(f3,1,2,'TileSpacing','compact','Padding','compact');
    for ic=1:2
        ax=nexttile(tl3); hold(ax,'on');
        for j=1:numel(exact{ic}.Delay)
            xline(ax,exact{ic}.Delay(j),':','Color',[.25 .25 .25], ...
                'HandleVisibility','off');
        end
        for j=1:numel(shown)
            rows=predicted{ic}(predicted{ic}.Order==shown(j),:);
            scatter(ax,rows.Delay,shown(j)*ones(height(rows),1),55,'filled', ...
                'DisplayName',sprintf('Padé %d/%d',shown(j),shown(j)));
        end
        scatter(ax,exact{ic}.Delay,13*ones(height(exact{ic}),1),70,'kx', ...
            'LineWidth',1.8,'DisplayName','Exact');
        grid(ax,'on'); box(ax,'on'); xlim(ax,[0 Tmax]); ylim(ax,[0 14]);
        yticks(ax,[shown 13]); yticklabels(ax,{'1','2','4','8','12','exact'});
        xlabel(ax,'Critical delay'); ylabel(ax,'Padé order'); title(ax,cases(ic).name);
        if ic==2, legend(ax,'Location','eastoutside'); end
    end
    title(tl3,'Predicted events: vertical lines are exact critical delays');
    save_figure(f3,figDir,'pade_critical_events');

    write_table(dataDir,cases,orders,exact,predicted,firstError,coverage);
    study=struct('cases',cases,'orders',orders,'Tmax',Tmax,'exact',{exact}, ...
        'predicted',{predicted},'firstRelativeError',firstError,'coverage',coverage);
end

function write_table(dataDir,cases,orders,exact,predicted,firstError,coverage)
    fid=fopen(fullfile(dataDir,'pade_critical_tables.tex'),'w');
    assert(fid>0); c=onCleanup(@() fclose(fid)); %#ok<NASGU>
    shown=[1 2 4 8 12];
    fprintf(fid,'%% Generated by run_pade_critical_study.m\n');
    fprintf(fid,'\\begin{table}[htbp]\\centering\\small\n');
    fprintf(fid,'\\caption{First critical delay: exact value and Padé prediction.}\\label{tab:padeTc}\n');
    fprintf(fid,'\\begin{tabular}{lrrr}\\toprule\n');
    fprintf(fid,'%s\n','Case and order & predicted $T_c$ & absolute error & relative error \\ \midrule');
    for ic=1:2
        exactFirst=exact{ic}.Delay(1); freq=exact{ic}.Frequency(1);
        row=sprintf('%s exact & %.9f & -- & --',cases(ic).id,exactFirst);
        fprintf(fid,'%s %s\n',row,'\\');
        for n=shown
            rows=predicted{ic}(predicted{ic}.Order==n & abs(predicted{ic}.Frequency-freq)<1e-7,:);
            if isempty(rows), continue; end
            err=abs(rows.Delay(1)-exactFirst);
            row=sprintf('\\quad Padé %d/%d & %.9f & %.3e & %.3e',n,n,rows.Delay(1),err,err/exactFirst);
            fprintf(fid,'%s %s\n',row,'\\');
        end
        if ic==1, fprintf(fid,'\\addlinespace\n'); end
    end
    fprintf(fid,'\\bottomrule\\end{tabular}\\end{table}\n\n');
    fprintf(fid,'\\begin{table}[htbp]\\centering\\small\n');
    fprintf(fid,'\\caption{Coverage of the exact events up to $T=10$.}\\label{tab:padeCoverage}\n');
    fprintf(fid,'\\begin{tabular}{rrrrr}\\toprule\n');
    fprintf(fid,'%s\n','Order & events P1 & coverage P1 & events P3 & coverage P3 \\ \midrule');
    for io=1:numel(orders)
        row=sprintf('%d & %d/%d & %.1f\\%% & %d/%d & %.1f\\%%',orders(io), ...
            round(coverage(1,io)*height(exact{1})),height(exact{1}),100*coverage(1,io), ...
            round(coverage(2,io)*height(exact{2})),height(exact{2}),100*coverage(2,io));
        fprintf(fid,'%s %s\n',row,'\\');
    end
    fprintf(fid,'\\bottomrule\\end{tabular}\\end{table}\n');
end

function save_figure(fig,figDir,name)
    exportgraphics(fig,fullfile(figDir,[name '.pdf']),'ContentType','vector');
    exportgraphics(fig,fullfile(figDir,[name '.png']),'Resolution',220);
end
