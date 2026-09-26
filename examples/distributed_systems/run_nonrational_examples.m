function results = run_nonrational_examples()
%RUN_NONRATIONAL_EXAMPLES Reproducible study of general nonrational spectra.
%   Writes tables, MATLAB data and PNG figures for the companion report.
    % Resolve paths from this file, so the study runs from any folder.
    here=fileparts(mfilename('fullpath')); repo=fileparts(fileparts(here));
    if exist('croots','file')~=2, run(fullfile(repo,'setup_contourroots.m')); end
    oldPath=addpath(fullfile(repo,'examples','models'));
    restorePath=onCleanup(@() path(oldPath));
    out=fullfile(repo,'output','distributed_systems');
    if ~isfolder(out), mkdir(out); end
    options={'AssumeAnalytic',true};
    models={ ...
        'heat_neumann', .5, 0, [-100 2 -3 3]; ...
        'heat_dirichlet', .5, 0, [-100 2 -3 3]; ...
        'heat_mixed', .5, 0, [-100 2 -3 3]; ...
        'wave', .5, 0, [-2 1 -30 30]; ...
        'wave_damped', .5, .5, [-2 1 -30 30]; ...
        'duct_constant', .5, 0, [-2 1 -20 20]; ...
        'duct_radiation', .5, 0, [-400 30 -2500 2500]; ...
        'beam', .5, .02, [-25 2 -160 160]};
    results=struct;
    for j=1:size(models,1)
        name=models{j,1}; kind=name;
        if strcmp(kind,'wave_damped'), kind='wave'; end
        [G,meta]=distributed_model(kind,'Sensor',models{j,2},'Damping',models{j,3});
        [p,info]=transfer_poles(G,models{j,4},options{:},'Singularities',meta.singularities);
        results.(name)=struct('poles',p,'info',info,'meta',meta);
        fprintf('%s: %d poles, %d cancelled locations, %s\n',name, ...
            numel(p),numel(info.cancelledLocations),info.status);
        assert(info.complete,'Unresolved example: %s',name);
        writetable(table(real(p),imag(p),info.multiplicity,info.residuals, ...
            'VariableNames',{'Real','Imag','Multiplicity','TargetResidual'}),fullfile(out,[name '.csv']));
    end

    % User example, evaluated without Symbolic Toolbox or Pade.
    Delta=@(s) 1+s+s.^2+(2*s+3).*exp(-s);
    [r,info]=characteristic_roots(Delta,[-8 2 -20 20],options{:});
    assert(info.complete);
    results.characteristic=struct('roots',r,'info',info);

    % Independent analytical/modal references, never used as Newton seeds.
    references={-(2*pi*(0:1)).^2, -(pi*[1 3]).^2, -(pi*((0:2)+.5)).^2};
    k=[-9:-1 1:9]; k=k(mod(k,4)~=0); references{4}=1i*pi*k;
    references{5}=[]; references{6}=log(.5)/2+1i*pi*(-6:6); references{7}=[];
    beta=zeros(4,1);
    for k=1:4
        beta(k)=fzero(@(b) cos(b)+1/cosh(b),[(k-1)*pi+.01 k*pi-.01]);
    end
    referenceBeam=[];
    for b=beta.'
        referenceBeam=[referenceBeam; roots([1 .02*b^4 b^4])]; %#ok<AGROW>
    end
    references{8}=referenceBeam(real(referenceBeam)>-25 & abs(imag(referenceBeam))<160);
    counts=zeros(8,1); cancellations=counts; maxerror=NaN(8,1);
    for j=1:8
        item=results.(models{j,1}); counts(j)=numel(item.poles);
        cancellations(j)=numel(item.info.cancelledLocations);
        ref=references{j};
        if ~isempty(ref)
            assert(numel(ref)==numel(item.poles));
            maxerror(j)=max(arrayfun(@(z) min(abs(item.poles-z)),ref));
            assert(maxerror(j)<1e-5);
        end
    end
    summary=table(string(models(:,1)),counts,cancellations,maxerror, ...
        'VariableNames',{'Model','Poles','CancelledLocations','MaxAnalyticError'});
    writetable(summary,fullfile(out,'summary.csv')); disp(summary);
    results.summary=summary;

    % Sensor-position experiment: cancelled transfer poles versus PDE modes.
    sensors=[.25 .5 .75 1]; sensorPoles=cell(3,numel(sensors));
    for j=1:3
        for k=1:numel(sensors)
            G=distributed_model(models{j,1},'Sensor',sensors(k));
            [sensorPoles{j,k},ii]=transfer_poles(G,[-100 2 -3 3],options{:});
            assert(ii.complete);
        end
    end
    results.sensorPoles=sensorPoles;
    colors=[.08 .35 .58; .84 .32 .12; .15 .55 .43; .55 .3 .65];
    f=figure('Visible','off','Color','w','Position',[100 100 1250 740]);
    tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
    titles={'Heat: Neumann','Heat: Dirichlet','Heat: mixed'};
    for j=1:3
        nexttile; hold on;
        for k=1:numel(sensors)
            rr=sensorPoles{j,k};
            scatter(real(rr),sensors(k)*ones(size(rr)),65,colors(k,:),'x','LineWidth',1.8);
        end
        xlim([-100 2]); ylim([.12 1.12]); yticks(sensors); grid on; box on;
        title(titles{j}); xlabel('Re(s)'); ylabel('Sensor x_0/L');
    end
    nexttile; hold on;
    plot(real(r),imag(r),'x','Color',colors(1,:),'MarkerSize',9,'LineWidth',1.8);
    xline(0,':'); yline(0,':'); grid on; box on;
    title('User characteristic: 7 roots in the rectangle'); xlabel('Re(s)'); ylabel('Im(s)');
    exportgraphics(f,fullfile(out,'heat_and_characteristic.png'),'Resolution',180); close(f);

    f=figure('Visible','off','Color','w','Position',[100 100 1250 780]);
    tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
    nexttile; hold on;
    a=results.wave; b=results.wave_damped;
    plot(real(a.poles),imag(a.poles),'o','Color',colors(1,:),'MarkerSize',7);
    plot(real(b.poles),imag(b.poles),'x','Color',colors(2,:),'MarkerSize',8,'LineWidth',1.5);
    plot(real(a.info.cancelledLocations),imag(a.info.cancelledLocations),'+','Color',[.6 .6 .6],'MarkerSize',9);
    grid on; box on; xlabel('Re(s)'); ylabel('Im(s)'); title('Wave: feedback damping and cancellations');
    legend('Undamped','Damping = 0.5','Cancelled modes','Location','best');
    xlim([-.22 .015]);
    nexttile; hold on;
    p=results.duct_constant.poles;
    plot(real(p),imag(p),'x','Color',colors(1,:),'MarkerSize',8,'LineWidth',1.5);
    xline(log(.5)/2,'--','Analytic line'); grid on; box on;
    xlim([-.55 -.15]); ylim([-20 20]);
    xlabel('Re(s)'); ylabel('Im(s)'); title('Duct: constant reflectance, normalized units');
    nexttile; hold on;
    p=results.duct_radiation.poles;
    plot(real(p),imag(p),'x','Color',colors(2,:),'MarkerSize',8,'LineWidth',1.5);
    xline(0,':'); grid on; box on; xlabel('Re(s) [1/s]'); ylabel('Im(s) [rad/s]');
    title('Duct: radiation impedance, article Fig. 5 parameters');
    nexttile; hold on;
    p=results.beam.poles;
    plot(real(p),imag(p),'x','Color',colors(3,:),'MarkerSize',9,'LineWidth',1.7);
    plot(real(references{8}),imag(references{8}),'o','Color',colors(1,:),'MarkerSize',10);
    grid on; box on; xlabel('Re(s)'); ylabel('Im(s)'); title('Beam: direct search and independent modal formula');
    legend('Nonrational search','Modal reference','Location','best');
    exportgraphics(f,fullfile(out,'wave_duct_beam.png'),'Resolution',180); close(f);

    % Accumulation is illustrated using the independent modal formula.
    k=(1:80).'; beta=(k-.5)*pi;
    for j=1:8
        beta(j)=fzero(@(b) cos(b)+1/cosh(b),[(j-1)*pi+.01 j*pi-.01]);
    end
    b4=beta.^4; discriminant=sqrt(complex(.02^2*b4.^2-4*b4));
    % Stable quadratic formula: avoid cancellation on the slow branch.
    fast=(-.02*b4-discriminant)/2; slow=b4./fast;
    results.beamAccumulation=struct('mode',k,'slow',slow,'fast',fast);
    f=figure('Visible','off','Color','w','Position',[100 100 1150 440]);
    tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
    nexttile; semilogy(k,abs(slow+50),'LineWidth',1.7,'Color',colors(1,:));
    grid on; xlabel('Mode index k'); ylabel('|s_{slow,k} + 50|'); title('Poles accumulate at s = -50');
    nexttile; plot(real(slow),imag(slow),'o','Color',colors(3,:),'MarkerSize',4);
    xline(-50,'--','Accumulation'); grid on; xlabel('Re(s)'); ylabel('Im(s)');
    title('Every neighborhood of -50 contains infinitely many poles');
    exportgraphics(f,fullfile(out,'beam_accumulation.png'),'Resolution',180); close(f);
    save(fullfile(out,'results.mat'),'results');
end
