function results=run_coupled_beam_study()
%RUN_COUPLED_BEAM_STUDY Exact characteristic poles versus m, k and c.
%   All values are nondimensional (L=EI=rho*A=1). Each sample repeats the
%   contour completeness check. Previous poles are only supplementary seeds.
    % Resolve paths from this file, so the study runs from any folder.
    here=fileparts(mfilename('fullpath')); repo=fileparts(fileparts(here));
    if exist('croots','file')~=2, run(fullfile(repo,'setup_contourroots.m')); end
    oldPath=addpath(fullfile(repo,'examples','models'));
    restorePath=onCleanup(@() path(oldPath));
    out=fullfile(repo,'output','coupled_beam'); if ~isfolder(out), mkdir(out); end
    region=[-60 1 -150 150];
    beta=fzero(@(b) cos(b)+1/cosh(b),[1 2]); omega1=beta^2;
    m0=.1; k0=m0*omega1^2; c0=.05;
    [D,G,meta]=coupled_beam_model('Mass',m0,'Stiffness',k0,'Damping',c0);
    [p,info]=characteristic_roots(D,region,'AssumeAnalytic',true);
    assert(info.complete);
    results=struct('region',region,'bareFirstFrequency',omega1, ...
        'baseline',struct('poles',p,'info',info,'parameters',meta.parameters));
    writePoles(fullfile(out,'baseline.csv'),p,meta);

    names={'Mass','Stiffness','Damping'};
    values={unique([linspace(.02,.5,25) m0]), ...
        unique([linspace(.05,5,35) k0]),unique([linspace(0,1.5,41) c0 .15 .3 .6 1])};
    for j=1:3
        rootsByValue=cell(numel(values{j}),1); alpha=zeros(size(values{j}));
        residual=alpha; count=alpha; seeds=p;
        rows=zeros(0,6);
        for k=1:numel(values{j})
            v=[m0 k0 c0]; v(j)=values{j}(k);
            [D,~,mi]=coupled_beam_model('Mass',v(1),'Stiffness',v(2),'Damping',v(3));
            [r,ii]=characteristic_roots(D,region,'AssumeAnalytic',true,'SeedPoints',seeds);
            assert(ii.complete,'Incomplete %s=%g',names{j},v(j));
            assert(max(real(r))<1e-7,'Unexpected unstable pole in passive model.');
            rootsByValue{k}=r; seeds=r;
            alpha(k)=max(real(r)); residual(k)=max(mi.normalizedResidual(r)); count(k)=sum(ii.multiplicity);
            assert(residual(k)<1e-7);
            rows=[rows; [repmat(v(j),numel(r),1) real(r) imag(r) ii.multiplicity ...
                mi.normalizedResidual(r) -real(r)./abs(r)]]; %#ok<AGROW>
        end
        results.(lower(names{j}))=struct('values',values{j},'poles',{rootsByValue}, ...
            'windowAbscissa',alpha,'normalizedResidual',residual,'count',count);
        writetable(array2table(rows,'VariableNames', ...
            {'Parameter','Real','Imag','Multiplicity','NormalizedResidual','DampingRatio'}), ...
            fullfile(out,[lower(names{j}) '_sweep.csv']));
        fprintf('%s: %d samples; count range %d..%d; max normalized residual %.3g\n', ...
            names{j},numel(values{j}),min(count),max(count),max(residual));
    end

    % Quantitative FE convergence at the baseline; no FE seeds in the search.
    meshes=[8 16 32 64]; error=zeros(size(meshes));
    for j=1:numel(meshes)
        pf=coupled_beam_fem(meta.parameters,meshes(j));
        error(j)=max(arrayfun(@(z) min(abs(pf-z)),p));
    end
    assert(all(diff(error)<0));
    results.fem=table(meshes.',error.','VariableNames',{'Elements','MaxPoleError'});
    writetable(results.fem,fullfile(out,'fem_convergence.csv'));

    % Selected damping values for an interpretable table; list every low pole.
    selected=[0 .05 .15 .3 .6 1 1.5]; rows=zeros(0,4);
    for c=selected
        [~,idx]=min(abs(results.damping.values-c)); r=results.damping.poles{idx};
        low=r(imag(r)>1e-7 & imag(r)<12 | abs(imag(r))<=1e-7);
        [~,order]=sort(abs(low)); low=low(order);
        rows=[rows; [repmat(c,numel(low),1) real(low) imag(low) -real(low)./abs(low)]]; %#ok<AGROW>
    end
    selectedTable=array2table(rows,'VariableNames',{'Damping','Real','Imag','DampingRatio'});
    writetable(selectedTable,fullfile(out,'selected_damping.csv')); disp(selectedTable);
    [bestAlpha,bestIdx]=min(results.damping.windowAbscissa);
    results.bestSample=struct('damping',results.damping.values(bestIdx),'windowAbscissa',bestAlpha);

    % Check whether enlarging the search window changes the reported abscissa.
    windowRows=zeros(3,4); selectedC=[c0 results.bestSample.damping 1.5];
    for j=1:3
        [Dc,~,~]=coupled_beam_model('Mass',m0,'Stiffness',k0,'Damping',selectedC(j));
        [re,ie]=characteristic_roots(Dc,[-90 1 -300 300],'AssumeAnalytic',true);
        assert(ie.complete);
        [~,idx]=min(abs(results.damping.values-selectedC(j)));
        smallAlpha=results.damping.windowAbscissa(idx); bigAlpha=max(real(re));
        assert(abs(smallAlpha-bigAlpha)<1e-7);
        windowRows(j,:)=[selectedC(j) numel(re) smallAlpha bigAlpha];
    end
    results.windowCheck=array2table(windowRows,'VariableNames', ...
        {'Damping','ExpandedWindowCount','OriginalAbscissa','ExpandedAbscissa'});
    writetable(results.windowCheck,fullfile(out,'window_check.csv'));

    colors=[.08 .35 .58;.85 .33 .12;.15 .55 .43];
    f=figure('Visible','off','Color','w','Position',[50 50 1500 480]);
    tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
    for j=1:3
        nexttile; hold on; set(gca,'FontSize',12);
        dat=results.(lower(names{j})); allLow=[];
        for k=1:numel(dat.values)
            r=dat.poles{k}; r=r(imag(r)>=-1e-7 & imag(r)<12);
            allLow=[allLow;r]; %#ok<AGROW>
            scatter(real(r),max(0,imag(r)),24,repmat(dat.values(k),numel(r),1),'filled');
        end
        colormap(parula); cb=colorbar; cb.Label.String=names{j};
        xline(0,':'); grid on; box on; xlabel('Re(s)'); ylabel('Im(s)');
        title([names{j} ': low-frequency poles']); ylim([-.25 12]);
        left=min(real(allLow)); xlim([min(-.05,1.08*left) max(.015,-.03*left)]);
    end
    exportgraphics(f,fullfile(out,'pole_sweeps.png'),'Resolution',190); close(f);

    f=figure('Visible','off','Color','w','Position',[50 50 1250 800]);
    tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
    for j=1:2
        nexttile; hold on; set(gca,'FontSize',12);
        dat=results.(lower(names{j})); w=zeros(numel(dat.values),2);
        for k=1:numel(dat.values)
            rp=dat.poles{k}; rp=sort(imag(rp(imag(rp)>1e-7))); w(k,:)=rp(1:2).';
        end
        plot(dat.values,w(:,1),'-','Color',colors(1,:),'LineWidth',1.8);
        plot(dat.values,w(:,2),'-','Color',colors(2,:),'LineWidth',1.8);
        if j==1, isolated=sqrt(k0./dat.values); else, isolated=sqrt(dat.values/m0); end
        plot(dat.values,isolated,'--','Color',[.5 .5 .5],'LineWidth',1);
        yline(omega1,':','Bare beam'); grid on; box on; xlabel(names{j}); ylabel('Im(s)');
        title('Two lowest oscillatory branches'); legend('Lower','Upper','Isolated oscillator','Location','best');
    end
    nexttile; set(gca,'FontSize',12);
    plot(results.damping.values,results.damping.windowAbscissa,'-','Color',colors(3,:),'LineWidth',1.8);
    hold on; plot(results.bestSample.damping,bestAlpha,'o','Color',colors(2,:),'MarkerFaceColor',colors(2,:));
    grid on; box on; xlabel('Damping c'); ylabel('max Re(s) in search window');
    title('Damping is not monotonically beneficial');
    nexttile; loglog(meshes,error,'o-','Color',colors(1,:),'LineWidth',1.8);
    grid on; box on; xlabel('Euler-Bernoulli finite elements'); ylabel('Max absolute pole error');
    title('Independent FEM convergence (10 poles)');
    exportgraphics(f,fullfile(out,'trends_and_validation.png'),'Resolution',190); close(f);

    % Mode illustrations for the conservative system, exact spatial solution.
    [D,~,~]=coupled_beam_model('Mass',m0,'Stiffness',k0,'Damping',0);
    [r,ii]=characteristic_roots(D,region,'AssumeAnalytic',true); assert(ii.complete);
    freq=sort(imag(r(imag(r)>1e-7))); x=linspace(0,1,300);
    f=figure('Visible','off','Color','w','Position',[50 50 1150 460]);
    tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
    modeRatios=zeros(2,1);
    for j=1:2
        b=sqrt(freq(j)); sigma=(cosh(b)+cos(b))/(sinh(b)+sin(b));
        shape=cosh(b*x)-cos(b*x)-sigma*(sinh(b*x)-sin(b*x)); shape=shape/shape(end);
        ratio=k0/(k0-m0*freq(j)^2); modeRatios(j)=ratio;
        nexttile; hold on; set(gca,'FontSize',12);
        plot(x,shape,'Color',colors(j,:),'LineWidth',2);
        plot([1 1.25],[1 ratio],':','Color',[.4 .4 .4],'LineWidth',1.5);
        plot(1.25,ratio,'s','MarkerSize',12,'MarkerFaceColor',colors(j,:),'Color',colors(j,:));
        yline(0,':'); grid on; box on; xlim([0 1.4]);
        xlabel('Beam coordinate x/L; mass at schematic x=1.25'); ylabel('Displacement normalized by tip');
        title(sprintf('Mode %d: omega=%.4f; z/y=%.3f',j,freq(j),ratio));
    end
    results.conservativeModes=struct('frequency',freq(1:2),'massToTipRatio',modeRatios);
    exportgraphics(f,fullfile(out,'coupled_modes.png'),'Resolution',190); close(f);

    % Machine-generated LaTeX tables keep the report tied to actual results.
    fid=fopen(fullfile(out,'tables.tex'),'w'); assert(fid>=0); cleaner=onCleanup(@() fclose(fid));
    fprintf(fid,'\\newcommand{\\BestDamping}{%.4g}\n\\newcommand{\\BestAlpha}{%.6f}\n',results.bestSample.damping,bestAlpha);
    fprintf(fid,'\\newcommand{\\FirstBareFrequency}{%.8f}\n',omega1);
    fprintf(fid,'\\newcommand{\\BaselineStiffness}{%.8f}\n',k0);
    fprintf(fid,'\\newcommand{\\BaselinePoleTable}{\\begin{tabular}{rrr}\\toprule Re(s)&Im(s)&$\\zeta$\\\\\\midrule\n');
    pp=p(imag(p)>0); [~,order]=sort(imag(pp)); pp=pp(order);
    for z=pp.', fprintf(fid,'%.8f & %.8f & %.6f\\\\\n',real(z),imag(z),-real(z)/abs(z)); end
    fprintf(fid,'\\bottomrule\\end{tabular}}\n');
    fprintf(fid,'\\newcommand{\\DampingPoleTable}{\\begin{tabular}{rrrr}\\toprule $c$&Re(s)&Im(s)&$\\zeta$\\\\\\midrule\n');
    for j=1:height(selectedTable)
        fprintf(fid,'%.2f & %.6f & %.6f & %.4f\\\\\n',selectedTable{j,:});
    end
    fprintf(fid,'\\bottomrule\\end{tabular}}\n');
    fprintf(fid,'\\newcommand{\\FEMTable}{\\begin{tabular}{rr}\\toprule Elements&Max. absolute error\\\\\\midrule\n');
    for j=1:numel(meshes), fprintf(fid,'%d & %.6g\\\\\n',meshes(j),error(j)); end
    fprintf(fid,'\\bottomrule\\end{tabular}}\n');
    save(fullfile(out,'results.mat'),'results');
    fprintf('Best sampled c=%.4g, window alpha=%.6g.\n',results.bestSample.damping,bestAlpha);
end

function writePoles(path,p,meta)
    writetable(table(real(p),imag(p),-real(p)./abs(p),meta.normalizedResidual(p), ...
        'VariableNames',{'Real','Imag','DampingRatio','NormalizedResidual'}),path);
end
