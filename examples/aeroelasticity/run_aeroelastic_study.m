function results=run_aeroelastic_study(outputDir)
%RUN_AEROELASTIC_STUDY Reproducible 2-DOF exact-Theodorsen study.
%   No optional toolbox needed. Writes figures, CSV and MAT to outputDir.
%   Default: output/aeroelasticity. All searches exclude the branch cut.
    % Resolve paths from this file, so the study runs from any folder.
    here=fileparts(mfilename('fullpath')); root=fileparts(fileparts(here));
    if exist('croots','file')~=2, run(fullfile(root,'setup_contourroots.m')); end
    oldPath=addpath(fullfile(root,'examples','models'));
    restorePath=onCleanup(@() path(oldPath));
    if nargin==0, outputDir=fullfile(root,'output','aeroelasticity'); end
    if ~isfolder(outputDir), mkdir(outputDir); end
    results=struct;
    methods={'exact','jones'};
    colors=[0 .35 .65;.85 .35 .05;.35 .55 .25];
    for kind={'nasa','dlr'}
        name=kind{1}; m=aeroelastic_section(name);
        if strcmp(name,'nasa'), speeds=0:5:250; else, speeds=0:5:300; end
        n=numel(speeds); poles=complex(zeros(n,2,2)); pk=complex(zeros(n,2));
        counts=NaN(n,2); residual=zeros(n,2); candidates=cell(2,1);
        for im=1:2
            method=methods{im}; prev=[];
            for j=1:n
                U=speeds(j);
                if im==1
                    [r,info]=croots(@(s) aeroelastic_delta(s,U,m,method), ...
                        [-180 100 .1 230],'AssumeAnalytic',true,'SeedPoints',prev);
                    assert(info.complete && info.count==2 && numel(r)==2, ...
                        'aeroelastic:Incomplete','%s at U=%g: %s, count=%g',name,U,info.status,info.count);
                    counts(j,im)=info.count;
                else
                    allr=aeroelastic_rational_roots(U,m,'jones');
                    r=allr(imag(allr)>1e-6); assert(numel(r)==2);
                end
                if isempty(prev), [~,ix]=sort(imag(r)); r=r(ix); else, r=match_two(r,prev); end
                poles(j,:,im)=r.'; prev=r;
                for k=1:2
                    d=aeroelastic_matrix(r(k),U,m,method); sig=svd(d);
                    residual(j,im)=max(residual(j,im),sig(end)/sig(1));
                end
                if im==1
                    if j==1, pk(j,:)=r.'; else
                        [pkr,~]=aeroelastic_pk_roots(U,m,pk(j-1,:).');
                        pk(j,:)=match_two(pkr,pk(j-1,:).').';
                    end
                end
            end
            candidates{im}=aeroelastic_flutter(m,method);
            fprintf('%s %s: Uf=%.10g, omega=%.10g, k=%.10g\n', ...
                name,method,candidates{im}.U,candidates{im}.omega,candidates{im}.k);
        end
        quasiAlpha=arrayfun(@(u) max(real(aeroelastic_rational_roots(u,m,'quasisteady'))),speeds);
        if max(real(aeroelastic_rational_roots(.01,m,'quasisteady')))>1e-8
            fq=struct('U',NaN,'omega',NaN,'k',NaN);
            quasiNote='Already unstable at U=0.01; no positive onset reported';
        else
            fq=aeroelastic_flutter(m,'quasisteady'); quasiNote='Positive-speed flutter candidate';
        end
        % Recheck several full regions with doubled contour resolution.
        for j=unique(round(linspace(1,n,6)))
            [r,i]=croots(@(s) aeroelastic_delta(s,speeds(j),m),[-220 130 .05 300], ...
                'AssumeAnalytic',true,'ContourPoints',64,'RootTolerance',1e-8);
            assert(i.complete && i.count==2 && numel(r)==2);
            for z=poles(j,:,1), assert(min(abs(r-z))<1e-6); end
        end
        f=candidates{1};
        % Stability count: roots with Re(s)>0 in a right-half-plane window,
        % which contains the positive real axis (no branch cut there).
        rhp=[1e-3 60 -250 250]; unstable=NaN(n,1); rhpStatus=strings(n,1);
        for j=2:n    % U=0 is conservative: its roots lie on the imaginary axis
            [r,i]=croots(@(s) aeroelastic_delta(s,speeds(j),m),rhp, ...
                'AssumeAnalytic',true,'ContourRefinements',10,'Warn',false);
            rhpStatus(j)=i.status;
            if i.complete, unstable(j)=i.count; end
            expected=2*(speeds(j)>f.U);
            % Inconclusive is acceptable only when a mode is within 0.1 of the
            % imaginary axis (very low speed, or close to flutter).
            nearAxis=min(abs(real(poles(j,:,1))))<.1;
            assert(unstable(j)==expected || (~i.complete && nearAxis), ...
                'aeroelastic:StabilityCount','%s at U=%g: count %g, status %s', ...
                name,speeds(j),unstable(j),i.status);
        end
        % Above the static divergence speed a REAL unstable root appears; it
        % lies on the positive real axis, inside the right-half-plane window.
        divergence=sqrt(m.K(2,2)/(2*pi*m.rho*m.b^2*(.5+m.a)));
        [rd,id]=croots(@(s) aeroelastic_delta(s,1.05*divergence,m),[1e-3 200 -300 300], ...
            'AssumeAnalytic',true,'ContourRefinements',10);
        assert(id.complete && any(abs(imag(rd))<1e-8 & real(rd)>0), ...
            'aeroelastic:Divergence','No real unstable root above divergence.');
        realRoot=real(rd(abs(imag(rd))<1e-8));
        harmonic=aeroelastic_harmonic_flutter(m);
        assert(abs(harmonic.U-f.U)<1e-6 && abs(harmonic.k-f.k)<1e-8);
        assert(abs(f.U-m.reference.U)<m.reference.velocityTolerance);
        if isfinite(m.reference.k), assert(abs(f.k-m.reference.k)<m.reference.frequencyTolerance); end
        fpk=aeroelastic_flutter(m,'pk'); assert(abs(fpk.U-f.U)<1e-6);
        tab=table(string(methods(:)),cellfun(@(x)x.U,candidates), ...
            cellfun(@(x)x.omega,candidates),cellfun(@(x)x.k,candidates), ...
            'VariableNames',{'Method','FlutterSpeed','AngularFrequency','ReducedFrequency'});
        tab.RelativeSpeedError=tab.FlutterSpeed/f.U-1;
        tab=[tab;table("quasisteady",fq.U,fq.omega,fq.k,fq.U/f.U-1, ...
            'VariableNames',tab.Properties.VariableNames)];
        writetable(tab,fullfile(outputDir,[name '_flutter.csv']));
        trajectory=table(speeds.',real(poles(:,1,1)),imag(poles(:,1,1)), ...
            real(poles(:,2,1)),imag(poles(:,2,1)),real(pk(:,1)),imag(pk(:,1)), ...
            real(pk(:,2)),imag(pk(:,2)),counts(:,1),residual(:,1), ...
            'VariableNames',{'U','Sigma1','Omega1','Sigma2','Omega2','PKSigma1','PKOmega1', ...
            'PKSigma2','PKOmega2','UpperCount','MaxModalResidual'});
        trajectory.QuasiSteadyMaxGrowth=quasiAlpha.';
        trajectory.UnstableRootsRHP=unstable;
        trajectory.RHPStatus=rhpStatus;
        writetable(trajectory,fullfile(outputDir,[name '_trajectories.csv']));
        results.(name)=struct('model',m,'speeds',speeds,'poles',poles,'pk',pk, ...
            'counts',counts,'modalResiduals',residual,'flutter',tab,'divergence',divergence, ...
            'quasiSteadyMaxGrowth',quasiAlpha,'quasiSteadyNote',quasiNote, ...
            'harmonicFlutter',harmonic,'unstableCount',unstable, ...
            'divergenceCheck',struct('U',1.05*divergence,'realRoot',realRoot));
        fprintf('%s: divergence %.6g; at 1.05*Ud a real root s = %.6g\n',name,divergence,realRoot);
        fig=figure('Visible','off','Color','w','Position',[100 100 1100 740]);
        tl=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
        ax=nexttile(tl); hold(ax,'on');
        plot(real(poles(:,:,1)),imag(poles(:,:,1)),'-','Color',colors(1,:),'LineWidth',1.7);
        plot(real(pk),imag(pk),'--','Color',colors(2,:),'LineWidth',1.2);
        xline(0,':'); plot(0,f.omega,'kp','MarkerFaceColor',[1 .75 0],'MarkerSize',11);
        grid on; xlabel('Growth rate [1/s]'); ylabel('Frequency [rad/s]');
        title('Root loci: full lines exact; dashed p-k');
        ax=nexttile(tl); hold(ax,'on');
        plot(speeds,imag(poles(:,:,1)),'-','Color',colors(1,:),'LineWidth',1.7);
        plot(speeds,imag(pk),'--','Color',colors(2,:),'LineWidth',1.2);
        xline(f.U,':'); grid on; xlabel(['U [' m.velocityUnits ']']); ylabel('Frequency [rad/s]');
        title('Both oscillatory branches');
        ax=nexttile(tl); hold(ax,'on'); handles=gobjects(3,1);
        for im=1:2
            handles(im)=plot(speeds,max(real(poles(:,:,im)),[],2),'Color',colors(im,:),'LineWidth',1.6);
        end
        handles(3)=plot(speeds,quasiAlpha,'Color',colors(3,:),'LineWidth',1.6);
        yline(0,':'); xline(f.U,':'); grid on; xlabel(['U [' m.velocityUnits ']']); ylabel('Max Re(s) in window [1/s]');
        legend(handles,{'Exact','Jones two-lag','Quasi-steady C=1'},'Location','best');
        title('Aerodynamic approximation changes flutter');
        ax=nexttile(tl); hold(ax,'on');
        plot(speeds,real(pk-poles(:,:,1)),'LineWidth',1.5); yline(0,':'); xline(f.U,':');
        grid on; xlabel(['U [' m.velocityUnits ']']); ylabel('p-k minus exact growth rate [1/s]');
        title('Same neutral boundary; different damping away from it');
        legend('Lower-frequency branch at U=0','Higher-frequency branch at U=0','Location','best');
        title(tl,sprintf('%s benchmark | exact Uf = %.6g, kf = %.6g',upper(name),f.U,f.k));
        set(findall(fig,'Type','axes'),'FontSize',11);
        exportgraphics(fig,fullfile(outputDir,[name '_study.png']),'Resolution',170); close(fig);
    end
    k=logspace(-3,1.5,350); c=theodorsen_laplace(1i*k);
    j=1-.165*(1i*k)./(1i*k+.0455)-.335*(1i*k)./(1i*k+.3);
    fig=figure('Visible','off','Color','w','Position',[100 100 1100 360]);
    tl=tiledlayout(fig,1,3,'TileSpacing','compact');
    nexttile; semilogx(k,real(c),k,-imag(c),k,real(j),'--',k,-imag(j),'--','LineWidth',1.5);
    grid on; xlabel('k'); title('Theodorsen and Jones'); legend('Re C','-Im C','Re Jones','-Im Jones');
    nexttile; loglog(k,abs(j-c),'LineWidth',1.7); grid on; xlabel('k'); ylabel('|C_J-C|'); title('Aerodynamic error');
    nexttile; hold on;
    patch([.06 .95 .95 .06],[-1.6 -1.6 1.6 1.6],[1 .9 .88],'EdgeColor',[.8 .3 .1],'LineWidth',1.8);
    rectangle('Position',[-1.7 .12 2.3 1.6],'EdgeColor',colors(1,:),'LineWidth',1.8);
    plot([-2 0],[0 0],'r','LineWidth',3); plot(0,0,'ro','MarkerFaceColor','r');
    text(-1.95,-.25,'branch cut','Color','r','FontSize',10);
    text(-1.6,1.55,'modes','Color',colors(1,:),'FontSize',10);
    text(.12,-1.4,'Re(s)>0','Color',[.8 .3 .1],'FontSize',10);
    xline(0,':'); yline(0,':'); xlim([-2 1]); ylim([-1.8 1.9]); grid on;
    xlabel('Re(s)'); ylabel('Im(s)'); title('Valid search regions');
    set(findall(fig,'Type','axes'),'FontSize',11);
    exportgraphics(fig,fullfile(outputDir,'theodorsen_and_domain.png'),'Resolution',170); close(fig);
    % Stiffness-ratio study, keeping mass, inertia, geometry and omega_alpha fixed.
    m=aeroelastic_section('nasa'); ratios=(.3:.05:.8).'; uf=zeros(size(ratios)); kf=uf;
    wa=sqrt(m.K(2,2)/m.M(2,2));
    for j=1:numel(ratios)
        p=m; p.K(1,1)=p.M(1,1)*(wa*ratios(j))^2;
        f=aeroelastic_flutter(p); uf(j)=f.U; kf(j)=f.k;
        for factor=[.98 1.02]
            [r,i]=croots(@(s) aeroelastic_delta(s,f.U*factor,p),[-220 130 .05 300],'AssumeAnalytic',true);
            assert(i.complete && i.count==2 && sign(max(real(r)))==sign(factor-1));
        end
    end
    results.parameterStudy=table(ratios,uf,kf,'VariableNames',{'FrequencyRatio','FlutterSpeed','ReducedFrequency'});
    writetable(results.parameterStudy,fullfile(outputDir,'nasa_frequency_ratio.csv'));
    fig=figure('Visible','off','Color','w','Position',[100 100 850 350]);
    tiledlayout(1,2); nexttile; plot(ratios,uf,'o-','LineWidth',1.6); grid on;
    xlabel('Uncoupled omega_h / omega_alpha'); ylabel('Flutter U [ft/s]');
    nexttile; plot(ratios,kf,'o-','LineWidth',1.6); grid on;
    xlabel('Uncoupled omega_h / omega_alpha'); ylabel('Flutter reduced frequency');
    set(findall(fig,'Type','axes'),'FontSize',11,'XLim',[.28 .82]);
    exportgraphics(fig,fullfile(outputDir,'nasa_parameter_study.png'),'Resolution',170); close(fig);
    save(fullfile(outputDir,'aeroelastic_results.mat'),'results');
end
function r=match_two(r,previous)
    if sum(abs(r-previous))>sum(abs(flipud(r)-previous)), r=flipud(r); end
end
