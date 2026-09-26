function results = run_delay_study()
%RUN_DELAY_STUDY Reproduce all numerical experiments and report figures.
%   Outputs are written below output/data and output/figures.

    % Resolve paths from this file, so the study runs from any folder.
    here=fileparts(mfilename('fullpath')); repo=fileparts(fileparts(here));
    if exist('croots','file')~=2, run(fullfile(repo,'setup_contourroots.m')); end
    rootDir = fullfile(repo,'output','time_delay');
    figDir = fullfile(rootDir,'figures');
    dataDir = fullfile(rootDir,'data');
    if ~exist(figDir,'dir'), mkdir(figDir); end
    if ~exist(dataDir,'dir'), mkdir(dataDir); end

    oldDefaults = set_plot_defaults();
    cleanup = onCleanup(@() restore_plot_defaults(oldDefaults)); %#ok<NASGU>
    close all

    problems = define_problems();
    orders = [1 2 4 6 8 10 12];
    nominal = repmat(struct(),numel(problems),1);
    criticalResults = cell(numel(problems),1);

    fprintf('Computing nominal pole maps...\n');
    f1 = figure('Color','w','Position',[80 80 1250 440]);
    tl = tiledlayout(f1,1,numel(problems),'TileSpacing','compact','Padding','compact');
    for ip = 1:numel(problems)
        p = problems(ip);
        [rex,info] = delay_roots(p.N,p.D,p.T0,p.xlimits,p.ylimits,50,86, ...
            'PadeSeedOrders',[4 8 12 16],'MaxRefinements',2,'VerifyCount',true);
        nominal(ip).exact = rex;
        nominal(ip).info = info;
        nominal(ip).pade = cell(numel(orders),1);
        [criticalResults{ip},criticalInfo] = critical_delays( ...
            p.N,p.D,max(p.Ts),'Tmin',min(p.Ts));
        nominal(ip).criticalInfo = criticalInfo;
        writetable(criticalResults{ip}, ...
            fullfile(dataDir,[p.id '_critical_delays.csv']));
        ax = nexttile(tl); hold(ax,'on');
        plot(ax,real(rex),imag(rex),'x','Color',[0.04 0.25 0.52], ...
            'MarkerSize',7,'LineWidth',1.8,'DisplayName','Exact');
        shownOrders = [1 4 8];
        palette = [0.86 0.24 0.18; 0.96 0.58 0.12; 0.10 0.58 0.48];
        for io = 1:numel(orders)
            rp = roots(pade_characteristic(p.N,p.D,p.T0,orders(io)));
            nominal(ip).pade{io} = rp;
            jj = find(shownOrders==orders(io),1);
            if ~isempty(jj)
                keep = inside(rp,p.xlimits,p.ylimits);
                plot(ax,real(rp(keep)),imag(rp(keep)),'o','Color',palette(jj,:), ...
                    'MarkerSize',5.5,'LineWidth',1.1, ...
                    'DisplayName',sprintf('Padé %d/%d',orders(io),orders(io)));
            end
        end
        xline(ax,0,'--','Color',[0.35 0.35 0.35],'HandleVisibility','off');
        yline(ax,0,':','Color',[0.55 0.55 0.55],'HandleVisibility','off');
        grid(ax,'on'); box(ax,'on');
        xlim(ax,p.plotX); ylim(ax,p.plotY);
        xlabel(ax,'Re(s)'); ylabel(ax,'Im(s)');
        title(ax,sprintf('%s, T = %.2g',p.shortName,p.T0));
        if ip == numel(problems), legend(ax,'Location','eastoutside'); end

        write_nominal_csv(dataDir,p,rex,info,nominal(ip).pade,orders);
    end
    title(tl,'Exact roots and diagonal Padé poles');
    save_figure(f1,figDir,'pole_maps');

    fprintf('Sweeping delay values...\n');
    f2 = figure('Color','w','Position',[100 100 1250 440]);
    tl2 = tiledlayout(f2,1,numel(problems),'TileSpacing','compact','Padding','compact');
    sweepResults = cell(numel(problems),1);
    for ip = 1:numel(problems)
        p = problems(ip);
        Ts = p.Ts;
        alphaExact = NaN(size(Ts));
        alphaPade1 = NaN(size(Ts));
        alphaPade4 = NaN(size(Ts));
        alphaPade8 = NaN(size(Ts));
        seed = [];
        for it = 1:numel(Ts)
            rp8 = roots(pade_characteristic(p.N,p.D,Ts(it),8));
            [rr,~] = delay_roots(p.N,p.D,Ts(it),p.sweepX,p.sweepY,25,44, ...
                'SeedPoints',[seed;rp8],'PadeSeedOrders',[4 8 12], ...
                'MaxRefinements',1,'VerifyCount',false);
            if ~isempty(rr)
                alphaExact(it) = max(real(rr));
                seed = rr;
            end
            alphaPade1(it) = max(real(roots(pade_characteristic(p.N,p.D,Ts(it),1))));
            alphaPade4(it) = max(real(roots(pade_characteristic(p.N,p.D,Ts(it),4))));
            alphaPade8(it) = max(real(rp8));
        end
        tbl = table(Ts(:),alphaExact(:),alphaPade1(:),alphaPade4(:),alphaPade8(:), ...
            'VariableNames',{'T','Exact','Pade1','Pade4','Pade8'});
        writetable(tbl,fullfile(dataDir,[p.id '_delay_sweep.csv']));
        sweepResults{ip} = tbl;

        ax = nexttile(tl2); hold(ax,'on');
        plot(ax,Ts,alphaExact,'-','Color',[0.04 0.25 0.52], ...
            'LineWidth',2.3,'DisplayName','Exact');
        plot(ax,Ts,alphaPade1,'--','Color',[0.86 0.24 0.18], ...
            'LineWidth',1.5,'DisplayName','Padé 1/1');
        plot(ax,Ts,alphaPade4,'-.','Color',[0.96 0.58 0.12], ...
            'LineWidth',1.6,'DisplayName','Padé 4/4');
        plot(ax,Ts,alphaPade8,':','Color',[0.10 0.58 0.48], ...
            'LineWidth',2.0,'DisplayName','Padé 8/8');
        crit = criticalResults{ip};
        if ~isempty(crit)
            scatter(ax,crit.Delay,zeros(height(crit),1),52,'d','filled', ...
                'MarkerFaceColor',[0.48 0.16 0.62],'MarkerEdgeColor','w', ...
                'LineWidth',0.7,'DisplayName','Exact critical delay');
        end
        yline(ax,0,'k--','HandleVisibility','off');
        grid(ax,'on'); box(ax,'on');
        xlabel(ax,'Delay T'); ylabel(ax,'Spectral abscissa');
        title(ax,p.shortName);
        if ip == numel(problems), legend(ax,'Location','eastoutside'); end
    end
    title(tl2,'Migration of the dominant root with the delay');
    save_figure(f2,figDir,'spectral_abscissa');

    fprintf('Computing Pade convergence metrics...\n');
    f3 = figure('Color','w','Position',[130 130 900 580]);
    ax = axes(f3); hold(ax,'on');
    colors = lines(numel(problems));
    convergence = cell(numel(problems),1);
    for ip = 1:numel(problems)
        p = problems(ip);
        rex = nominal(ip).exact;
        % Compare the same dominant conjugate pair at every order.
        dominant = rex(1:min(2,numel(rex)));
        maxError = NaN(size(orders));
        rmsError = NaN(size(orders));
        spurious = NaN(size(orders));
        for io = 1:numel(orders)
            rp = nominal(ip).pade{io};
            comp = match_pade_roots(dominant,rp,numel(dominant));
            finiteError = comp.AbsoluteError(isfinite(comp.AbsoluteError));
            if ~isempty(finiteError)
                maxError(io) = max(finiteError);
                rmsError(io) = sqrt(mean(finiteError.^2));
            end
            nearWindow = inside(rp,p.plotX,p.plotY);
            spurious(io) = max(0,sum(nearWindow)-sum(isfinite(comp.AbsoluteError)));
        end
        convergence{ip} = table(orders(:),rmsError(:),maxError(:),spurious(:), ...
            'VariableNames',{'Order','RMSError','MaxError','UnmatchedPadePoles'});
        writetable(convergence{ip},fullfile(dataDir,[p.id '_pade_convergence.csv']));
        semilogy(ax,orders,rmsError,'-o','Color',colors(ip,:), ...
            'MarkerFaceColor',colors(ip,:),'LineWidth',1.8, ...
            'DisplayName',p.shortName);
    end
    grid(ax,'on'); box(ax,'on');
    set(ax,'YScale','log');
    ylim(ax,[1e-16 1]);
    xlabel(ax,'Order n of the [n/n] approximation');
    ylabel(ax,'RMS error of the dominant roots');
    title(ax,'Convergence of the dominant roots with the Padé order');
    legend(ax,'Location','southwest');
    save_figure(f3,figDir,'pade_convergence');

    f4 = figure('Color','w','Position',[160 160 1150 400]);
    tl4 = tiledlayout(f4,1,numel(problems),'TileSpacing','compact','Padding','compact');
    for ip=1:numel(problems)
        ax = nexttile(tl4); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
        crit = criticalResults{ip};
        if isempty(crit)
            text(ax,0.5,0.5,'No crossing in the range', ...
                'Units','normalized','HorizontalAlignment','center', ...
                'FontSize',12,'Color',[0.25 0.25 0.25]);
            xlim(ax,[min(problems(ip).Ts) max(problems(ip).Ts)]); ylim(ax,[0 1]);
        else
            destab = crit.CrossingSpeed>0;
            scatter(ax,crit.Delay(destab),crit.Frequency(destab),75,'^','filled', ...
                'MarkerFaceColor',[0.83 0.20 0.16],'DisplayName','Destabilizing');
            scatter(ax,crit.Delay(~destab),crit.Frequency(~destab),75,'v','filled', ...
                'MarkerFaceColor',[0.10 0.55 0.42],'DisplayName','Stabilizing');
            xlim(ax,[min(problems(ip).Ts) max(problems(ip).Ts)]);
            ylim(ax,[0 1.18*max(crit.Frequency)]);
        end
        xlabel(ax,'Critical delay T_c'); ylabel(ax,'Frequency \omega_c');
        title(ax,problems(ip).shortName);
        if ip==numel(problems) && ~isempty(crit), legend(ax,'Location','best'); end
    end
    title(tl4,'Exact imaginary-axis crossings');
    save_figure(f4,figDir,'critical_delays');

    write_latex_tables(rootDir,problems,nominal,orders,convergence,sweepResults,criticalResults);
    fprintf('Studying Pade predictions of critical delays...\n');
    padeCriticalStudy = run_pade_critical_study(rootDir);
    fprintf('Running Pade pitfall experiments...\n');
    run_pade_pitfalls(rootDir);
    results = struct('problems',problems,'nominal',nominal, ...
        'orders',orders,'sweeps',{sweepResults},'convergence',{convergence}, ...
        'criticalDelays',{criticalResults},'padeCriticalStudy',padeCriticalStudy);
    save(fullfile(dataDir,'study_results.mat'),'results');
    fprintf('Study complete. Figures: %s\n',figDir);
end

function problems = define_problems()
    problems(1) = struct('id','first_order_destabilizing', ...
        'shortName','P1: s+1+2e^{-sT}', ...
        'description','first order, destabilizing delayed feedback', ...
        'N',2,'D',[1 1],'T0',1,'xlimits',[-14 3],'ylimits',[-42 42], ...
        'plotX',[-7 2],'plotY',[-22 22],'sweepX',[-14 4],'sweepY',[-45 45], ...
        'Ts',linspace(0.1,2.4,24));
    problems(2) = struct('id','first_order_stable', ...
        'shortName','P2: s+2+e^{-sT}', ...
        'description','first order, stable for every delay', ...
        'N',1,'D',[1 2],'T0',1,'xlimits',[-14 2],'ylimits',[-42 42], ...
        'plotX',[-7 1],'plotY',[-22 22],'sweepX',[-14 3],'sweepY',[-45 45], ...
        'Ts',linspace(0.1,3.0,24));
    problems(3) = struct('id','second_order', ...
        'shortName','P3: s^2+0.8s+4+3e^{-sT}', ...
        'description','oscillatory second order with a delayed term', ...
        'N',3,'D',[1 0.8 4],'T0',0.8,'xlimits',[-14 3],'ylimits',[-45 45], ...
        'plotX',[-7 2],'plotY',[-24 24],'sweepX',[-14 4],'sweepY',[-48 48], ...
        'Ts',unique([linspace(0.08,3.2,40) linspace(2.62,2.9,29)]));
end

function tf = inside(z,xl,yl)
    tf = real(z)>=xl(1) & real(z)<=xl(2) & imag(z)>=yl(1) & imag(z)<=yl(2);
end

function write_nominal_csv(dataDir,p,rex,info,padeRoots,orders)
    rootTable = table(real(rex),imag(rex),info.residuals, ...
        'VariableNames',{'Real','Imaginary','Residual'});
    writetable(rootTable,fullfile(dataDir,[p.id '_exact_roots.csv']));
    for k = 1:numel(orders)
        rp = padeRoots{k};
        tbl = table(real(rp),imag(rp),'VariableNames',{'Real','Imaginary'});
        writetable(tbl,fullfile(dataDir,sprintf('%s_pade_%02d.csv',p.id,orders(k))));
    end
end

function write_latex_tables(rootDir,problems,nominal,orders,convergence,sweeps,criticalResults)
    % One file per table so that the report can place each one in context.
    fid = open_table(rootDir,'table_solver.tex');
    fprintf(fid,'\\begin{table}[htbp]\\centering\\small\n');
    fprintf(fid,'\\caption{Solver diagnostics at the nominal delay.}\\label{tab:solver}\n');
    fprintf(fid,'\\begin{tabular}{lrrrr}\\toprule\n');
    fprintf(fid,'%s\n','Case & Roots & Count & Max residual & Abscissa \\ \midrule');
    for i=1:numel(problems)
        infi = nominal(i).info; rr = nominal(i).exact;
        row = sprintf('P%d & %d & %d & %.2e & %.5f', ...
            i,numel(rr),infi.argumentPrincipleCount,max(infi.residuals),max(real(rr)));
        fprintf(fid,'%s %s\n',row,'\\');
    end
    fprintf(fid,'\\bottomrule\\end{tabular}\\end{table}\n\n');

    fclose(fid); fid = open_table(rootDir,'table_pade_rms.tex');
    fprintf(fid,'\\begin{table}[htbp]\\centering\\small\n');
    fprintf(fid,'\\caption{RMS error of the dominant pair at the nominal delay.}\\label{tab:pade}\n');
    fprintf(fid,'\\begin{tabular}{l%s}\\toprule\n',repmat('r',1,numel(orders)));
    fprintf(fid,'Case'); fprintf(fid,' & $n=%d$',orders);
    fprintf(fid,'%s\n',' \\ \midrule');
    for i=1:numel(problems)
        row = sprintf('P%d',i);
        for j=1:numel(orders)
            row = [row sprintf(' & %.2e',convergence{i}.RMSError(j))]; %#ok<AGROW>
        end
        fprintf(fid,'%s %s\n',row,'\\');
    end
    fprintf(fid,'\\bottomrule\\end{tabular}\\end{table}\n\n');

    fclose(fid); fid = open_table(rootDir,'table_sweep.tex');
    fprintf(fid,'\\begin{table}[htbp]\\centering\\small\n');
    fprintf(fid,'\\caption{Range of the exact spectral abscissa in the delay sweep.}\\label{tab:sweep}\n');
    fprintf(fid,'\\begin{tabular}{lr@{\\hspace{2.5em}}r@{\\hspace{2.5em}}r}\\toprule\n');
    fprintf(fid,'%s\n','Case & $T_{min}$--$T_{max}$ & $\min\alpha$ & $\max\alpha$ \\ \midrule');
    for i=1:numel(problems)
        t=sweeps{i};
        row = sprintf('P%d & %.2f--%.2f & %.4f & %.4f', ...
            i,min(t.T),max(t.T),min(t.Exact),max(t.Exact));
        fprintf(fid,'%s %s\n',row,'\\');
    end
    fprintf(fid,'\\bottomrule\\end{tabular}\\end{table}\n');

    fclose(fid); fid = open_table(rootDir,'table_critical.tex');
    fprintf(fid,'\n\\begin{table}[htbp]\\centering\\small\n');
    fprintf(fid,'\\caption{Imaginary-axis crossings in the studied delay range.}\\label{tab:critical}\n');
    fprintf(fid,'\\begin{tabular}{lrrrrl}\\toprule\n');
    fprintf(fid,'%s\n','Case & $T_c$ & $\omega_c$ & branch & $\operatorname{Re}(ds/dT)$ & direction \\ \midrule');
    for i=1:numel(problems)
        ctab=criticalResults{i};
        if isempty(ctab)
            fprintf(fid,'P%d & %s %s\n',i,'\multicolumn{5}{c}{no crossing}','\\');
        else
            for j=1:height(ctab)
                if ctab.CrossingSpeed(j)>0, direction='destabilizing'; else, direction='stabilizing'; end
                row=sprintf('P%d & %.8f & %.8f & %d & %+.3e & %s',i, ...
                    ctab.Delay(j),ctab.Frequency(j),ctab.Branch(j), ...
                    ctab.CrossingSpeed(j),direction);
                fprintf(fid,'%s %s\n',row,'\\');
            end
        end
    end
    fprintf(fid,'\\bottomrule\\end{tabular}\\end{table}\n');
    fclose(fid);
end

function fid = open_table(rootDir,name)
    out = fullfile(rootDir,'data',name);
    fid = fopen(out,'w');
    assert(fid>0,'Could not open %s.',out);
    fprintf(fid,'%% Generated by run_delay_study.m -- do not edit manually.\n');
end

function save_figure(fig,figDir,name)
    exportgraphics(fig,fullfile(figDir,[name '.pdf']),'ContentType','vector');
    exportgraphics(fig,fullfile(figDir,[name '.png']),'Resolution',220);
end

function old = set_plot_defaults()
    fields = {'defaultAxesFontName','defaultAxesFontSize','defaultAxesLineWidth', ...
              'defaultLineLineWidth','defaultFigureColor'};
    old = cell(size(fields));
    for k=1:numel(fields), old{k}=get(groot,fields{k}); end
    set(groot,'defaultAxesFontName','Helvetica','defaultAxesFontSize',18, ...
        'defaultAxesLineWidth',0.8,'defaultLineLineWidth',1.5, ...
        'defaultFigureColor','w');
    old = struct('fields',{fields},'values',{old});
end

function restore_plot_defaults(old)
    for k=1:numel(old.fields), set(groot,old.fields{k},old.values{k}); end
end
