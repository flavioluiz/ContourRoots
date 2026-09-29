function report=benchmark_shared_grid(mode,outputDir,repetitions,caseNames)
%BENCHMARK_SHARED_GRID Compare actual matrix work and frozen-bank outputs.
%   MODE='legacy' also runs on 0.8.0 (no SharedGrid option); 'shared' runs
%   0.9.0. Use the same MATLAB thread setting for both. No RNG is changed.
%   Optional CASENAMES resumes a subset, retaining other saved cases.
%   Dense case: de Hoog, t=0:.5:1, 16-node solve warmup; other cases: FFT.
%   Nonrational hybrid internal factorizations cannot be counted here.
    if nargin<1, mode='legacy'; end
    if nargin<2, outputDir=fullfile('output','shared_grid'); end
    if nargin<3, repetitions=3; end
    if nargin<4, caseNames={'small','hybrid','dense300'}; end
    mode=validatestring(mode,{'legacy','shared'});
    validateattributes(repetitions,{'numeric'},{'scalar','integer','positive','finite'});
    root=fileparts(fileparts(mfilename('fullpath'))); old=path;
    cleanup=onCleanup(@() path(old));
    addpath(fullfile(root,'matlab'),fullfile(root,'examples','models'));
    if ~isfolder(outputDir), mkdir(outputDir); end
    report=struct('version',contourroots_version,'matlab',version,'computer',computer, ...
        'threads',maxNumCompThreads,'mode',mode,'repetitions',repetitions,'cases',struct);
    file=fullfile(outputDir,[mode '.mat']);
    if isfile(file), prior=load(file,'report'); report.cases=prior.report.cases; end
    extra={}; if strcmp(mode,'shared'), extra={'SharedGrid',true}; end
    names=cellstr(caseNames);
    for test=1:numel(names)
        seen=complex(zeros(0,1));
        switch names{test}
            case 'small'
                M=cdyn(@small,eye(2),eye(2),zeros(2),'Dimensions',[2 2 2]);
                t=(0:.1:1).'; opts={'SingularityBound',0}; U=[1+t sin(t)];
            case 'hybrid'
                wing=wing_model('goland');
                M=cmimo(@hybrid,'Size',[2 2]); t=(0:.01:.2).';
                opts={'SingularityBound',5,'AbsTol',1e-3,'RelTol',1e-3}; U=[1+.2*sin(10*t) .3*cos(5*t)];
            case 'cheapChannels'
                M=cmimo({ndpair(1,[1 1]),2;0,ndpair(2,[1 2])});
                t=(0:.1:.6).'; opts={}; U=[1+t sin(t)];
            case 'dense300'
                n=300; Q=sqrt(2/(n+1))*sin(pi*(1:n)'*(1:n)/(n+1));
                A=Q*diag(linspace(.5,20,n))*Q'; Id=eye(n);
                B=[sin((1:n)'*.7) cos((1:n)'*.3)]/sqrt(n);
                C=cos((1:30)'*(1:n)*.13)/sqrt(n); C=C.*logspace(-2,2,30)';
                M=cdyn(@dense,B,C,zeros(30,2),'Dimensions',[n 30 2]);
                t=(0:.5:1).'; opts={'Method','dehoog','SingularityBound',0,'AbsTol',logspace(-6,-4,30)','RelTol',1e-5}; U=[1+t sin(t)];
            otherwise
                error('benchmark_shared_grid:Case','Unknown case %s.',names{test});
        end
        fprintf('%s %s warmup\n',mode,names{test});
        if strcmp(names{test},'dense300'), ceval(M,1+1i*(0:15));
        else, ckernel(M,t,opts{:},extra{:}); end
        elapsed=zeros(1,repetitions);
        for repeat=1:repetitions
            seen=complex(zeros(0,1)); timer=tic;
            [K,info]=ckernel(M,t,opts{:},extra{:}); elapsed(repeat)=toc(timer);
            fprintf('%s %s run %d: %.3f s, %d evaluations\n',mode,names{test},repeat,elapsed(repeat),info.evaluations);
        end
        [y,~,out]=clsim(K,U,t,'AbsTol',.02,'RelTol',2e-3);
        assert(out.converged); distinct=numel(unique(seen));
        result=struct('seconds',elapsed,'medianSeconds',median(elapsed), ...
            'distinctNodes',distinct,'info',info,'cacheHitRate',info.cacheHits/(info.evaluations+info.cacheHits), ...
            'time',t,'options',{opts},'output',y,'maxOutputDifference',NaN);
        if strcmp(mode,'shared')
            before=load(fullfile(outputDir,'legacy.mat'),'report');
            if isfield(before.report.cases,names{test})
                result.maxOutputDifference=max(abs(y-before.report.cases.(names{test}).output),[],'all');
            end
        end
        if strcmp(names{test},'dense300')
            reference=foh_reference(linspace(.5,20,n)',C*Q,Q'*B,U,t);
            result.maxReferenceError=max(abs(y-reference),[],'all');
        end
        report.cases.(names{test})=result;
        save(fullfile(outputDir,[mode '.mat']),'report');
        fid=fopen(fullfile(outputDir,[mode '.json']),'w');
        fprintf(fid,'%s',jsonencode(report,PrettyPrint=true)); fclose(fid);
    end
    function H=small(s), seen(end+1,1)=s; H=diag([s+1 s+2]); end
    function H=dense(s)
        seen(end+1,1)=s; H=s*Id+A;
        if mod(numel(seen),10000)==0, fprintf('%s dense300: %d node calls\n',mode,numel(seen)); end
    end
    function G=hybrid(s)
        seen(end+1,1)=s; G=1e6*wing_transfer_matrix(s,120,wing,'hybrid'); G=G([1 3],[1 3]);
    end
end

function y=foh_reference(rates,C,B,u,t)
% Exact modal state update for the same piecewise-linear sampled input.
    state=zeros(numel(rates),1); y=zeros(numel(t),size(C,1));
    for k=2:numel(t)
        dt=t(k)-t(k-1); decay=exp(-rates*dt);
        first=-expm1(-rates*dt)./rates;
        second=(dt-first)./rates;
        state=decay.*state+first.*(B*u(k-1,:)')+second.*(B*(u(k,:)-u(k-1,:))'/dt);
        y(k,:)=(C*state)';
    end
end
