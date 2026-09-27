function report=compare_performance_runs(beforeDir,afterDir)
%COMPARE_PERFORMANCE_RUNS Numerical equality and warm timing ratios.
%   Inputs are CORE output folders from BENCHMARK_COSTS. No timing threshold
%   determines correctness; compare spectra and simulated frozen kernels.
    root=fileparts(fileparts(mfilename('fullpath')));
    oldPath=addpath(fullfile(root,'matlab')); cleanupPath=onCleanup(@() path(oldPath));
    a=load(fullfile(beforeDir,'core_timings.mat')); b=load(fullfile(afterDir,'core_timings.mat'));
    report=struct('timings',struct,'equivalence',struct);
    names=fieldnames(a.report.stages);
    for k=1:numel(names)
        key=names{k}; before=a.report.stages.(key).medianSeconds; after=b.report.stages.(key).medianSeconds;
        report.timings.(key)=struct('beforeSeconds',before,'afterSeconds',after,'speedup',before/after);
        x=load(fullfile(beforeDir,[key '.mat'])); y=load(fullfile(afterDir,[key '.mat']));
        if strcmp(key,'hybridKernel')
            t=x.result.Time; u=[1+.2*sin(10*t),.3*cos(5*t)];
            [p,~,ip]=clsim(x.result,u,t,'AbsTol',.01,'RelTol',2e-3);
            [q,~,iq]=clsim(y.result,u,t,'AbsTol',.01,'RelTol',2e-3);
            discrepancy=max(abs(p-q),[],'all');
            assert(ip.converged&&iq.converged&&discrepancy<1e-10);
            assert(isequal(ip.resolvedMask,iq.resolvedMask));
            assert(isequal(ip.errorEstimate,iq.errorEstimate));
            assert(x.diagnostics.evaluations==y.diagnostics.evaluations);
            assert(x.diagnostics.cacheHits==y.diagnostics.cacheHits);
        else
            assert(x.diagnostics.complete&&y.diagnostics.complete);
            assert(numel(x.result)==numel(y.result)); ref=x.result;
            discrepancy=0;
            for z=y.result(:).'
                [distance,j]=min(abs(z-ref)); discrepancy=max(discrepancy,distance/(1+abs(ref(j)))); ref(j)=[];
            end
            assert(discrepancy<1e-10);
            if strcmp(key,'matrixModes'), assert(x.diagnostics.count==y.diagnostics.count); end
        end
        report.equivalence.(key)=struct('maxDiscrepancy',discrepancy,'passed',true);
        fprintf('%s: %.3f -> %.3f s (%.2fx), discrepancy %.3g\n',key,before,after,before/after,discrepancy);
    end
    fid=fopen(fullfile(afterDir,'comparison.json'),'w'); assert(fid>=0);
    cleanup=onCleanup(@() fclose(fid)); fprintf(fid,'%s\n',jsonencode(report,'PrettyPrint',true));
end
