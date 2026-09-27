function plan = buildfile
%BUILDFILE Development tasks for ContourRoots (run with BUILDTOOL).
%   buildtool            same as "buildtool test"
%   buildtool test       unit and regression tests (fast)
%   buildtool docs       run every MATLAB block of README.md and docs/
%   buildtool examples   run the full example studies (slow, writes output/)
%   buildtool manual     run the studies, then compile the PDF manual (LaTeX)
%   buildtool release    tests, docs, studies, manual, then dist/ContourRoots.zip/.mltbx
%   buildtool package    archives ONLY; no tests/studies/manual regeneration
%   buildtool clean      delete dist/ and LaTeX build files
%
%   Users of the toolbox do not need any of this: see README.md.

    plan = buildplan(localfunctions);
    plan.DefaultTasks = "test";
    % A release regenerates every study, so that the manual never uses
    % results computed by an older version of the code.
    plan("manual").Dependencies = "examples";
    plan("release").Dependencies = ["test" "docs" "manual"];
end

function testTask(~)
% Run the unit and regression tests.
    addpath(fullfile(repo_root(),'tests'));
    timed('tests',@() run_tests('unit','regression'));
end

function docsTask(~)
% Execute the README, the Markdown documentation and the quick-start examples.
    addpath(fullfile(repo_root(),'tests'));
    timed('docs',@() run_tests('docs'));
end

function examplesTask(~)
% Run the full example studies (several minutes); results go to output/.
    root = repo_root();
    run(fullfile(root,'setup_contourroots.m'));
    studies = {fullfile('time_delay','run_delay_study.m'), ...
        fullfile('distributed_systems','run_nonrational_examples.m'), ...
        fullfile('coupled_beam','run_coupled_beam_study.m')};
    for k = 1:numel(studies)
        fprintf('Running %s ...\n', studies{k});
        file=fullfile(root,'examples',studies{k});
        [~,name]=fileparts(file); timed(name,@() run_study(file));
    end
    addpath(fullfile(root,'examples','aeroelasticity'));
    timed('aeroelastic_study',@() run_aeroelastic_study);
    addpath(fullfile(root,'examples','time_response'));
    timed('time_response_study',@() run_time_response_study);
    addpath(fullfile(root,'examples','continuous_wing'));
    timed('continuous_wing_study',@() run_continuous_wing_study);
    timed('hybrid_comparison',@() run_hybrid_comparison);
    timed('matrix_modes_comparison',@() run_matrix_modes_comparison);
end

function manualTask(~)
% Compile manual/ContourRoots_manual.tex and copy the PDF to docs/.
    addpath(fullfile(repo_root(),'tools'));
    timed('manual',@() build_manual());
end

function releaseTask(~)
% Build dist/ContourRoots.zip and dist/ContourRoots.mltbx.
    addpath(fullfile(repo_root(),'tools'));
    timed('packaging',@() build_release());
end

function packageTask(~)
% Archives only. This is NOT a substitute for the release validation gates.
    fprintf('PACKAGING ONLY: no tests, studies or manual regeneration. Use release for validation.\n');
    addpath(fullfile(repo_root(),'tools'));
    timed('packaging',@() build_release());
end

function timed(name,operation)
    addpath(fullfile(repo_root(),'tools'));
    time_build_step(name,operation);
end
function run_study(file), run(file); end

function cleanTask(~)
% Remove release builds and LaTeX intermediate files.
    root = repo_root();
    if isfolder(fullfile(root,'dist')), rmdir(fullfile(root,'dist'),'s'); end
    if isfolder(fullfile(root,'output','manual')), rmdir(fullfile(root,'output','manual'),'s'); end
end

function root = repo_root()
    root = fileparts(mfilename('fullpath'));
end
