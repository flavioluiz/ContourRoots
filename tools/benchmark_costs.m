function report=benchmark_costs(root,outputDir,mode)
%BENCHMARK_COSTS Reproducible wall-time/profile study, not a release gate.
%   Use an immutable source copy for before/after comparisons. MODE is
%   'core', 'pipeline' or 'package'. Packaging writes only below OUTPUTDIR.
%   No artifact is published or installed. Timings depend on warm-up/load.
    if nargin<3, mode='core'; end
    if ~isfolder(outputDir), mkdir(outputDir); end
    oldPath=path; oldFolder=pwd; oldVisible=get(groot,'defaultFigureVisible');
    beforeFigures=findall(groot,'Type','figure');
    cleanup=onCleanup(@() restore_state(oldPath,oldFolder,oldVisible,beforeFigures));
    addpath(fullfile(root,'matlab'),fullfile(root,'matlab','delay'), ...
        fullfile(root,'examples','models'),fullfile(root,'tools'));
    cd(root); set(groot,'defaultFigureVisible','off');
    report=struct('matlab',version,'computer',computer,'sourceRoot',root, ...
        'mode',mode,'timestamp',char(datetime('now')),'stages',struct);
    switch mode
        case 'core'
            wing=wing_model('goland'); H=@(s) wing_matrix(s,150,wing,4);
            jobs={@() cmodes(H,[-30 15 1 320],'AssumeAnalytic',true), ...
                @() croots(@(s) wing_delta(s,150,wing),[-30 15 1 320],'AssumeAnalytic',true), ...
                @() ckernel(wing_response_model(120,wing,'hybrid',[1 3]),(0:.01:.2).', ...
                    'SingularityBound',5,'AbsTol',1e-3,'RelTol',1e-3)};
            names={'matrixModes','scalarModes','hybridKernel'};
            for j=1:numel(jobs)
                fprintf('BENCHMARK %s: warming up...\n',names{j}); jobs{j}();
                elapsed=zeros(1,3); result=[]; diagnostics=[];
                for k=1:3
                    timer=tic; [result,diagnostics]=jobs{j}(); elapsed(k)=toc(timer);
                    fprintf('BENCHMARK %s run %d: %.3f s\n',names{j},k,elapsed(k));
                end
                report.stages.(names{j})=struct('seconds',elapsed,'medianSeconds',median(elapsed));
                save(fullfile(outputDir,[names{j} '.mat']),'result','diagnostics');
                checkpoint();
            end
            profile clear; profile on
            jobs{3}();
            profile off; profiling=profile('info'); save(fullfile(outputDir,'kernel_profile.mat'),'profiling');
            [~,order]=sort([profiling.FunctionTable.TotalTime],'descend');
            for k=order(1:min(15,numel(order)))
                p=profiling.FunctionTable(k); fprintf('PROFILE %.3f s %d calls %s\n',p.TotalTime,p.NumCalls,p.FunctionName);
            end
        case 'pipeline'
            addpath(fullfile(root,'tests'));
            measure('unitRegression',@() run_tests('unit','regression'));
            measure('documentation',@() run_tests('docs'));
            scripts={'time_delay/run_delay_study.m','distributed_systems/run_nonrational_examples.m', ...
                'coupled_beam/run_coupled_beam_study.m'};
            names={'delayStudy','distributedStudy','coupledBeamStudy'};
            for k=1:numel(scripts)
                file=fullfile(root,'examples',scripts{k}); measure(names{k},@() run_script(file));
            end
            folders={'aeroelasticity','time_response','continuous_wing'};
            for k=1:numel(folders), addpath(fullfile(root,'examples',folders{k})); end
            measure('aeroelasticStudy',@() run_aeroelastic_study);
            measure('timeResponseStudy',@() run_time_response_study);
            measure('continuousWingStudy',@() run_continuous_wing_study);
            measure('hybridStudy',@() run_hybrid_comparison);
            measure('matrixModesStudy',@() run_matrix_modes_comparison);
        case 'package'
            % Reproduce the pre-optimization packaging stages explicitly.
            stage=fullfile(outputDir,'stage','ContourRoots');
            if isfolder(stage), error('benchmark_costs:Stage','Choose a fresh output directory.'); end
            mkdir(stage);
            items={'matlab','examples','docs','setup_contourroots.m','README.md','LICENSE', ...
                'THIRD_PARTY_NOTICES.md','CITATION.cff','CHANGELOG.md','CONTRIBUTING.md','ROADMAP.md'};
            timer=tic;
            for k=1:numel(items), copyfile(fullfile(root,items{k}),fullfile(stage,items{k})); end
            dev=fullfile(stage,'docs','development'); if isfolder(dev), rmdir(dev,'s'); end
            for pattern={'.DS_Store','*.asv'}
                d=dir(fullfile(stage,'**',pattern{1}));
                for k=1:numel(d), delete(fullfile(d(k).folder,d(k).name)); end
            end
            finish('staging',timer);
            timer=tic; zip(fullfile(outputDir,'ContourRoots.zip'),items,stage); finish('zip',timer);
            fprintf('BENCHMARK ToolboxOptions constructor...\n'); timer=tic;
            opts=matlab.addons.toolbox.ToolboxOptions(stage,'6f1b2e7c-3c1a-4a55-9b7e-0c4b8f2d9a31');
            finish('toolboxOptions',timer);
            timer=tic;
            opts.ToolboxName='ContourRoots'; opts.ToolboxVersion=contourroots_version;
            opts.AuthorName='Flávio Luiz Cardoso-Ribeiro';
            opts.Summary='Performance study of ContourRoots packaging.';
            opts.ToolboxMatlabPath={fullfile(stage,'matlab'),fullfile(stage,'matlab','delay')};
            opts.OutputFile=fullfile(outputDir,'ContourRoots.mltbx');
            finish('metadata',timer);
            fprintf('BENCHMARK packageToolbox...\n'); timer=tic;
            profile clear; profile on
            matlab.addons.toolbox.packageToolbox(opts);
            profile off; finish('mltbx',timer);
            profiling=profile('info'); save(fullfile(outputDir,'packaging_profile.mat'),'profiling');
            report.files=numel(opts.ToolboxFiles);
        otherwise
            error('benchmark_costs:Mode','Use core, pipeline or package.');
    end
    checkpoint();
    function measure(name,job)
        fprintf('BENCHMARK starting %s...\n',name); timer=tic; job(); finish(name,timer);
    end
    function finish(name,timer)
        seconds=toc(timer); report.stages.(name)=struct('seconds',seconds);
        fprintf('BENCHMARK completed %s: %.3f s\n',name,seconds); checkpoint();
    end
    function checkpoint()
        save(fullfile(outputDir,[mode '_timings.mat']),'report');
        fid=fopen(fullfile(outputDir,[mode '_timings.json']),'w');
        if fid<0, error('benchmark_costs:Write','Cannot write timing report.'); end
        closer=onCleanup(@() fclose(fid)); fprintf(fid,'%s\n',jsonencode(report,'PrettyPrint',true));
    end
end
function run_script(file), run(file); end
function restore_state(oldPath,oldFolder,oldVisible,beforeFigures)
    profile off; close(setdiff(findall(groot,'Type','figure'),beforeFigures));
    path(oldPath); cd(oldFolder); set(groot,'defaultFigureVisible',oldVisible);
end
