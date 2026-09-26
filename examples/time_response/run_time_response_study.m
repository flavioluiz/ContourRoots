function results=run_time_response_study(outputDir)
%RUN_TIME_RESPONSE_STUDY Reproducible nonrational transient checks.
%   Base MATLAB only; CSV/MAT/PNG outputs, no release or solver changes.
    % Resolve paths from this file, so the study runs from any folder.
    here=fileparts(mfilename('fullpath')); root=fileparts(fileparts(here));
    if exist('croots','file')~=2, run(fullfile(root,'setup_contourroots.m')); end
    oldPath=addpath(fullfile(root,'examples','models'));
    restorePath=onCleanup(@() path(oldPath));
    if nargin==0, outputDir=fullfile(root,'output','time_response'); end
    if ~isfolder(outputDir), mkdir(outputDir); end
    results=struct;
    t=(0:.025:6).';
    cases={@delay_feedback,@diffusion,@heat_transfer,@beam_transfer};
    names={'Delay feedback','Diffusion','Finite rod','Coupled beam'};
    results.summary=table;
    fig=figure('Visible','off','Color','w','Position',[100 100 1100 740]);
    tiledlayout(fig,2,2,'TileSpacing','compact');
    for j=1:4
        [G,reference,aux]=cases{j}(t);
        [y,~,i]=cstep(G,t,'SingularityBound',0);
        assert(i.converged,'%s: %s',names{j},i.stopReason);
        err=max(abs(y-reference));
        target=[2e-5 3e-6 3e-6 2e-4]; assert(err<target(j));
        ax=nexttile; plot(ax,t,y,'LineWidth',1.6); hold(ax,'on');
        plot(ax,t,reference,'k--','LineWidth',1); grid(ax,'on');
        xlabel(ax,'Time'); ylabel(ax,'Step response'); title(ax,sprintf('%s | max error %.2g',names{j},err));
        legend(ax,'Transfer inversion','Independent reference','Location','best');
        results.(sprintf('case%d',j))=struct('time',t,'response',y,'reference',reference,'info',i,'auxiliary',aux);
        row=table(string(names{j}),err,i.points,i.evaluations,'VariableNames',{'Case','MaxAbsoluteError','MaxPoints','Evaluations'});
        results.summary=[results.summary;row]; %#ok<AGROW>
    end
    set(findall(fig,'Type','axes'),'FontSize',11);
    exportgraphics(fig,fullfile(outputDir,'time_response_validation.png'),'Resolution',160); close(fig);
    % Aeroelastic transfer = pitch response to a pitch moment, not det(D).
    m=aeroelastic_section('nasa'); ta=(0:.002:1).';
    fig=figure('Visible','off','Color','w','Position',[100 100 1000 410]); tiledlayout(1,2);
    for j=1:2
        U=[160 180]; U=U(j); G=@(s) aero_channel(s,U,m);
        % Bound 5 is an explicit study assumption, not inferred from a finite pole plot.
        [y,~,i]=cstep(G,ta,'SingularityBound',5,'AbsTol',1e-10,'RelTol',1e-4);
        pick=unique([1 round(linspace(2,numel(ta),18))]);
        [check,~,ic]=cstep(G,ta(pick),'SingularityBound',5,'Method','quadrature','AbsTol',1e-10,'RelTol',1e-4);
        assert(i.converged && ic.converged); err=max(abs(y(pick)-check));
        fprintf('Aero U=%g: cross-error %.6e; FFT estimate %.6e; quadrature estimate %.6e\n',U,err,max(i.errorEstimate),max(ic.errorEstimate));
        snapshot=struct('time',ta,'y',y,'pick',pick,'check',check,'fft',i,'quadrature',ic);
        save(fullfile(outputDir,sprintf('aero_validation_%g.mat',U)),'snapshot');
        assert(err<2e-8,'Aeroelastic inversion discrepancy %.6e at U=%g.',err,U);
        nexttile; plot(ta,y,'LineWidth',1.5); hold on; plot(ta(pick),check,'ko'); grid on;
        xlabel('Time [s]'); ylabel('Pitch / unit moment'); title(sprintf('NASA U=%g ft/s',U));
        legend('Shifted FFT','Independent inversion','Location','best');
        results.(sprintf('aero%d',j))=struct('time',ta,'response',y,'info',i,'checkInfo',ic,'crossError',err);
    end
    set(findall(fig,'Type','axes'),'FontSize',11);
    exportgraphics(fig,fullfile(outputDir,'time_response_aeroelasticity.png'),'Resolution',160); close(fig);
    % Two refinement axes and a deliberately wrong imaginary-axis inverse.
    tf=(0:.02:4).';
    [y,~,i]=cstep(@(s) exp(-.713*s)./(s+1),tf,'SingularityBound',0);
    exact=(tf>=.713).*(1-exp(-max(0,tf-.713)));
    fig=figure('Visible','off','Color','w','Position',[100 100 1000 410]); tiledlayout(1,2);
    nexttile; semilogy(i.history(:,2),max(i.history(:,3:5),1e-16),'o-','LineWidth',1.5); grid on;
    xlabel('Internal FFT points'); ylabel('Refinement difference'); legend('Period','Bandwidth','Shift');
    title('Separate convergence checks');
    [g,~,ig]=cimpulse(ndpair(1,[1 -.3]),tf);
    N=8192; dt=.01; w=2*pi/(N*dt)*[0:N/2 -N/2+1:-1].';
    bad=real(ifft(1./(1i*w-.3)))/dt; % intentionally invalid causal inversion
    nexttile; plot(tf,g,tf,exp(.3*tf),'--',(0:400)'*dt,bad(1:401),':','LineWidth',1.5);
    grid on; legend('Shifted FFT','Analytical causal','Wrong imaginary-axis result','Location','best');
    xlabel('Time'); title('Unstable plant: contour choice matters');
    set(findall(fig,'Type','axes'),'FontSize',11);
    exportgraphics(fig,fullfile(outputDir,'time_response_convergence.png'),'Resolution',160); close(fig);
    results.refinement=struct('history',i.history,'delayError',max(abs(y-exact)),'unstableInfo',ig);
    writetable(results.summary,fullfile(outputDir,'validation.csv'));
    save(fullfile(outputDir,'time_response_results.mat'),'results');
    disp(results.summary);
end
function [G,y,aux]=delay_feedback(t)
    a=.5; T=1; G=@(s) 1./(s+1+a*exp(-T*s)); y=zeros(size(t));
    % Expansion in delayed powers of 1/(s+1), integrated analytically.
    for k=0:floor(max(t)/T)
        q=t-k*T; active=q>=0;
        y(active)=y(active)+(-a)^k*gammainc(q(active),k+1,'lower');
    end
    aux='Finite method-of-steps series (not root or rational fitting).';
end
function [G,y,aux]=diffusion(t)
    a=.7; G=@(s) exp(-a*sqrt(s)); y=zeros(size(t));
    y(2:end)=erfc(a./(2*sqrt(t(2:end)))); aux='Analytical heat kernel.';
end
function [G,y,aux]=heat_transfer(t)
    x=.5; G=distributed_model('heat_dirichlet','Sensor',x);
    n=(1:150).'; coeff=2*(-1).^n.*sin(n*pi*x)./(n*pi);
    y=x+sum(coeff.*exp(-(n*pi).^2*t.'),1).'; y(1)=0;
    aux='Separation-of-variables rod solution, 150 modes; exact initial value.';
end
function [G,y,aux]=beam_transfer(t)
    [~,trans,meta]=coupled_beam_model; G=trans.TipForceToTip;
    solutions=zeros(numel(t),2);
    for j=1:2
        elements=[16 32]; [~,M,C,K]=coupled_beam_fem(meta.parameters,elements(j));
        n=size(M,1); force=zeros(n,1); force(n-2)=1;
        A=[zeros(n) eye(n);-M\K -M\C]; B=[zeros(n,1);M\force];
        output=[force.' zeros(1,n)]; eq=-A\B;
        P=expm(A*(t(2)-t(1))); state=zeros(2*n,1);
        for k=2:numel(t), state=eq+P*(state-eq); solutions(k,j)=output*state; end
    end
    y=solutions(:,2); aux=max(abs(diff(solutions,1,2)));
    assert(aux<2e-4,'FEM transient has not converged.');
end
function g=aero_channel(s,U,m)
    g=arrayfun(@(z) component(z,U,m),s);
end
function g=component(s,U,m)
    q=aeroelastic_matrix(s,U,m)\[0;1]; g=q(2);
end
