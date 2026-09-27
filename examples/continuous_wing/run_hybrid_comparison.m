function results=run_hybrid_comparison(outputDir,channels)
%RUN_HYBRID_COMPARISON M1/M2/M4 acceptance: ports, transfers and time response.
%   Saves comparison data and a figure. Not a generic network assembler.
    root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
    addpath(fullfile(root,'matlab'),fullfile(root,'examples','models'));
    if nargin<1, outputDir=fullfile(root,'output','hybrid_comparison'); end
    if nargin<2, channels=[1 3]; end % all nine frequency channels are always checked
    if ~isfolder(outputDir), mkdir(outputDir); end
    samples=[0 2+20i 3+70i 5+200i 5+500i]; rows=[];
    for kind={'dry','goland','tapered'}
        for n=[1 2 4 8]
            wing=wing_model(kind{1},n); U=120; worst=0;
            for s=samples
                A=wing_transfer_matrix(s,U,wing,'implicit'); B=wing_transfer_matrix(s,U,wing,'hybrid');
                worst=max(worst,norm(A-B,'fro')/norm(A,'fro'));
            end
            assert(worst<1e-8);
            rows=[rows;{kind{1},n,worst}]; %#ok<AGROW>
        end
    end
    results.transfer=cell2table(rows,'VariableNames',{'Model','Strips','RelativeError'});
    wing=wing_model('goland',2); U=120;
    A=wing_response_model(U,wing,'implicit',channels);
    B=wing_response_model(U,wing,'hybrid',channels);
    t=(0:.01:.2).'; allInputs=[1+.2*sin(10*t),.1*t,.3*cos(5*t)];
    u=allInputs(:,channels); nc=numel(channels);
    opts={'SingularityBound',5,'AbsTol',1e-3,'RelTol',1e-3};
    fprintf('Preparing all implicit channels...\n');
    [KA,prepA]=ckernel(A,t,opts{:});
    fprintf('Preparing all hybrid channels...\n');
    [KB,prepB]=ckernel(B,t,opts{:});
    out={'AbsTol',.01,'RelTol',2e-3};
    [ya,~,ia]=clsim(KA,u,t,out{:}); [yb,~,ib]=clsim(KB,u,t,out{:});
    assert(ia.converged&&ib.converged&&ia.evaluations==0&&ib.evaluations==0);
    difference=max(abs(ya-yb),[],1); assert(all(difference<.01));
    % Independent calls through the old scalar engine, not the matrix bank.
    local=zeros(size(ya));
    fprintf('Checking independent scalar superposition...\n');
    for j=1:nc, for i=1:nc
        [v,~,iv]=clsim(@(s) 1e6*wing_transfer(s,U,wing,channels(i),channels(j)),u(:,j),t, ...
            'SingularityBound',5,'AbsTol',.001,'RelTol',1e-3);
        assert(iv.converged); local(:,i)=local(:,i)+v;
    end, end
    superpositionDifference=max(abs(ya-local),[],1); assert(all(superpositionDifference<.01));
    % A second input uses all stored channels, including a previously small one.
    other=zeros(size(u)); other(:,end)=1;
    [second,~,reuse]=clsim(KA,other,t,out{:});
    assert(reuse.converged&&reuse.evaluations==0);
    results.time=struct('t',t,'channels',channels,'inputs',u,'implicit',ya,'hybrid',yb,'scalarSuperposition',local, ...
        'difference',difference,'superpositionDifference',superpositionDifference,'secondInput',second);
    results.preparation=struct('implicit',prepA,'hybrid',prepB);
    results.metadata=struct('matlab',version,'toolbox',contourroots_version, ...
        'boundAssumption',5,'kernelAbsTol',1e-3,'kernelRelTol',1e-3,'outputAbsTol',.01,'outputRelTol',2e-3);
    save(fullfile(outputDir,'hybrid_comparison.mat'),'results');
    writetable(results.transfer,fullfile(outputDir,'transfer_comparison.csv'));
    fig=figure('Visible','off','Position',[100 100 1000 700]); cleanup=onCleanup(@() close(fig));
    tiledlayout(nc,2); labels={'deflection [mm]','slope [mrad]','twist [mrad]'};
    for i=1:nc
        nexttile; plot(t,ya(:,i),'LineWidth',1.5); hold on; plot(t,yb(:,i),'--','LineWidth',1.2);
        grid on; ylabel(labels{channels(i)}); xlabel('t [s]'); if i==1, legend('implicit','hybrid'); end
        nexttile; plot(t,ya(:,i)-yb(:,i),'LineWidth',1.3); grid on;
        ylabel('assembly difference'); xlabel('t [s]');
    end
    exportgraphics(fig,fullfile(outputDir,'hybrid_comparison.png'),'Resolution',180);
    % Superposition on a finer grid: what a MIMO response is made of.
    % One bank, three inputs (both loads, force only, torque only).
    ts=(0:.005:.4).'; us=[1+.2*sin(10*ts), .3*cos(5*ts)];
    P=wing_response_model(U,wing,'implicit',[1 3]);
    fprintf('Preparing the fine-grid bank for the superposition figure...\n');
    [KS,prepS]=ckernel(P,ts,opts{:});
    [yboth,~,iboth]=clsim(KS,us,ts,out{:});
    [yforce,~,iforce]=clsim(KS,[us(:,1) 0*ts],ts,out{:});
    [ytorque,~,itorque]=clsim(KS,[0*ts us(:,2)],ts,out{:});
    assert(iboth.converged&&iforce.converged&&itorque.converged);
    results.superposition=struct('t',ts,'inputs',us,'both',yboth,'force',yforce, ...
        'torque',ytorque,'residual',max(abs(yboth-yforce-ytorque),[],'all'), ...
        'preparation',prepS);
    save(fullfile(outputDir,'hybrid_comparison.mat'),'results');
    fig2=figure('Visible','off','Position',[100 100 1000 420]); cleanup2=onCleanup(@() close(fig2));
    tiledlayout(1,2); names={'tip deflection [mm]','tip twist [mrad]'};
    for i=1:2
        nexttile; plot(ts,yforce(:,i),'-','Color',[0 .35 .65],'LineWidth',1.1); hold on
        plot(ts,ytorque(:,i),'-','Color',[.85 .35 .05],'LineWidth',1.1);
        plot(ts,yboth(:,i),'k-','LineWidth',1.8); grid on
        xlabel('t [s]'); ylabel(names{i});
        if i==1, legend('from the force','from the torque','both (clsim)','Location','southeast'); end
    end
    exportgraphics(fig2,fullfile(outputDir,'hybrid_superposition.png'),'Resolution',160);
    fprintf('Maximum transfer discrepancy: %.3e\n',max(results.transfer.RelativeError));
    fprintf('Time discrepancy (selected outputs):'); fprintf(' %.3e',difference); fprintf('\n');
    fprintf('Scalar-superposition discrepancy:'); fprintf(' %.3e',superpositionDifference); fprintf('\n');
end
