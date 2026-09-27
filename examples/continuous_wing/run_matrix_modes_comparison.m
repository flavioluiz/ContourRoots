function results=run_matrix_modes_comparison(outputDir)
%RUN_MATRIX_MODES_COMPARISON M3: matrix modes versus scalar determinant roots.
%   Fixed-size implicit matrices only. Hybrid port matrices are meromorphic
%   coordinate charts and must not be passed to an analytic mode search.
    repo=fileparts(fileparts(fileparts(mfilename('fullpath'))));
    addpath(fullfile(repo,'matlab'),fullfile(repo,'examples','models'));
    if nargin<1, outputDir=fullfile(repo,'output','matrix_modes'); end
    if ~isfolder(outputDir), mkdir(outputDir); end
    rows=cell(6,12); runs=cell(1,6); row=0;
    for kind={'dry','goland'}
        wing=wing_model(kind{1},1);
        if strcmp(kind{1},'dry'), U=0; box=[-2 2 1 320]; else, U=150; box=[-40 15 1 320]; end
        reference=[];
        for pieces=[1 2 4]
            fprintf('CMODES: %s, %d exact pieces...\n',kind{1},pieces);
            H=@(s) wing_matrix(s,U,wing,pieces);
            clock=tic; [r,i]=cmodes(wing_cdyn(U,wing,pieces),box,'AssumeAnalytic',true);
            matrixSeconds=toc(clock);
            clock=tic; [d,j]=croots(@(s) det(full(H(s))),box,'AssumeAnalytic',true);
            determinantSeconds=toc(clock);
            assert(i.complete&&j.complete&&numel(r)==numel(d));
            [~,ix]=sort(imag(r)); r=r(ix); [~,ix]=sort(imag(d)); d=d(ix);
            errorDet=max(abs(r-d)./(1+abs(d))); assert(errorDet<1e-8);
            if isempty(reference), reference=r; end
            errorPieces=max(abs(r-reference)./(1+abs(reference))); assert(errorPieces<1e-8);
            oracleError=NaN;
            if strcmp(kind{1},'dry')
                e=wing.strips; bending=[1.875104068711961 4.694091132974174].^2*sqrt(e.EI/e.mu)/wing.L^2;
                torsion=[1 3]*pi/(2*wing.L)*sqrt(e.GJ/e.Ialpha);
                exact=1i*sort([bending torsion]).';
                oracleError=max(abs(r-exact)./(1+abs(exact))); assert(oracleError<1e-8);
            end
            pivots=[i.history.minPivotRatio]; pivots=pivots(isfinite(pivots));
            row=row+1;
            rows(row,:)={kind{1},pieces,6*(pieces+1),numel(r),errorDet,errorPieces, ...
                oracleError,max(i.residuals),min(pivots),i.evaluations,matrixSeconds,determinantSeconds};
            runs{row}=struct('model',kind{1},'pieces',pieces,'modes',r,'determinantRoots',d, ...
                'matrixInfo',i,'scalarInfo',j);
        end
    end
    results.wing=cell2table(rows,'VariableNames',{'Model','Pieces','Dimension','Count', ...
        'RelativeDetError','RelativeSubdivisionError','RelativeClosedFormError', ...
        'ScaledResidual','MinContourPivotRatio','MatrixEvaluations','MatrixSeconds','DeterminantSeconds'});
    results.runs=runs;
    stress={};
    for scale=[1e-200 1 1e200]
        H=@(s) scale*[0 s+2;s+1 0]; box=[-3 0 -1 1];
        [r,i]=cmodes(H,box,'AssumeAnalytic',true);
        assert(i.complete&&max(abs(r-[-1;-2]))<1e-10);
        scalarComplete=false;
        try
            [~,j]=croots(@(s) det(H(s)),box,'AssumeAnalytic',true);
            scalarStatus=j.status; scalarComplete=j.complete;
        catch err
            scalarStatus=err.identifier;
        end
        stress(end+1,:)={scale,det(H(.5)),i.count,i.complete,scalarComplete,scalarStatus}; %#ok<AGROW>
    end
    results.scaling=cell2table(stress,'VariableNames',{'Scale','RawDetAtReference','MatrixCount', ...
        'MatrixComplete','DeterminantComplete','DeterminantStatus'});
    results.metadata=struct('matlab',version,'toolbox',contourroots_version, ...
        'description','Fixed-dimension exact strip matrices; counts concern characteristic modes, not transfer poles.');
    save(fullfile(outputDir,'matrix_modes_comparison.mat'),'results');
    writetable(results.wing,fullfile(outputDir,'wing_comparison.csv'));
    writetable(results.scaling,fullfile(outputDir,'scaling_comparison.csv'));
    fig=figure('Visible','off','Position',[100 100 1100 440]); cleanup=onCleanup(@() close(fig));
    tiledlayout(1,2,'Padding','compact');
    nexttile; hold on
    for k=[1 4]
        r=runs{k}.modes; plot(real(r),imag(r),'x','MarkerSize',9,'LineWidth',1.8);
        d=runs{k}.determinantRoots; plot(real(d),imag(d),'o','MarkerSize',9,'LineWidth',1);
    end
    grid on; xline(0,':'); xlabel('Re(s) [1/s]'); ylabel('Im(s) [rad/s]');
    title('Continuous wing: characteristic modes');
    legend('dry: cmodes','dry: determinant','150 m/s: cmodes','150 m/s: determinant','Location','best');
    nexttile; hold on
    for k=[1 4]
        ix=k:k+2;
        semilogy(results.wing.Pieces(ix),max(results.wing.RelativeDetError(ix),eps),'o-','LineWidth',1.5);
    end
    set(gca,'YScale','log'); grid on; xlabel('Exact pieces per span'); ylabel('Relative location discrepancy');
    title('cmodes versus croots(det H)'); legend('dry','150 m/s','Location','best');
    exportgraphics(fig,fullfile(outputDir,'matrix_modes_comparison.png'),'Resolution',180);
    % Flutter mode shape: the null vector of K at the flutter root, propagated
    % exactly along the span (wing_mode_shape). No modal basis is involved.
    wing=wing_model('goland'); Uf=136.983977449;         % Tutorial 11
    [lf,lfi]=cmodes(@(s) full(wing_matrix(s,Uf,wing)),[-2 2 60 80],'AssumeAnalytic',true);
    assert(lfi.complete&&numel(lf)==1);
    y=linspace(0,wing.L,121); [w,alpha,z]=wing_mode_shape(lf,Uf,wing,lfi.rightVectors{1},y);
    b=wing.strips(1).b; c=w(end); w=w/c; ba=b*alpha/c;   % normalize: tip deflection = 1
    results.flutterMode=struct('lambda',lf,'U',Uf,'y',y,'w',w,'balpha',ba, ...
        'tipRatio',abs(ba(end)),'tipPhaseDeg',angle(ba(end))*180/pi, ...
        'tipLoadResidual',norm(z(4:6,end))/norm(z(4:6,1)),'info',lfi);
    save(fullfile(outputDir,'matrix_modes_comparison.mat'),'results');
    fig2=figure('Visible','off','Position',[100 100 1100 420]); cleanup2=onCleanup(@() close(fig2));
    tiledlayout(1,2,'Padding','compact');
    nexttile; plot(y,abs(w),'LineWidth',1.8); hold on; plot(y,abs(ba),'LineWidth',1.8); grid on
    xlim([0 wing.L]); xlabel('span position y [m]'); ylabel('amplitude (tip deflection = 1)');
    legend('bending |w|','twist b|\alpha|','Location','northwest');
    title(sprintf('Flutter mode at U = %.2f m/s, \\omega = %.2f rad/s',Uf,imag(lf)));
    nexttile; th=linspace(0,360,361);
    plot(th,real(w(end)*exp(1i*th*pi/180)),'LineWidth',1.8); hold on
    plot(th,real(ba(end)*exp(1i*th*pi/180)),'LineWidth',1.8); grid on; xlim([0 360]); xticks(0:90:360)
    xlabel('phase in one cycle [deg]'); ylabel('tip motion');
    legend('bending w','twist b\alpha','Location','southwest');
    title(sprintf('Tip: twist lags bending by %.0f deg',-results.flutterMode.tipPhaseDeg));
    exportgraphics(fig2,fullfile(outputDir,'wing_flutter_mode.png'),'Resolution',160);
    disp(results.wing); disp(results.scaling);
end
