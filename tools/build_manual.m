function build_manual()
%BUILD_MANUAL Compile the PDF manual and copy it to docs/.
%   Runs the example studies if their results are missing, regenerates the
%   documentation figures, compiles manual/ContourRoots_manual.tex with
%   latexmk (a TeX distribution must be installed) and copies the PDF to
%   docs/ContourRoots_manual.pdf.

    root = fileparts(fileparts(mfilename('fullpath')));
    needed = {fullfile(root,'output','time_delay','data','table_critical.tex'), ...
        fullfile(root,'output','distributed_systems','summary.csv'), ...
        fullfile(root,'output','coupled_beam','tables.tex')};
    if ~all(cellfun(@isfile,needed))
        fprintf('Study results missing: running the examples first.\n');
        run(fullfile(root,'setup_contourroots.m'));
        run(fullfile(root,'examples','time_delay','run_delay_study.m'));
        run(fullfile(root,'examples','distributed_systems','run_nonrational_examples.m'));
        run(fullfile(root,'examples','coupled_beam','run_coupled_beam_study.m'));
    end
    if ~isfile(fullfile(root,'output','aeroelasticity','nasa_study.png'))
        addpath(fullfile(root,'examples','aeroelasticity')); run_aeroelastic_study;
    end
    if ~isfile(fullfile(root,'output','time_response','time_response_results.mat'))
        addpath(fullfile(root,'examples','time_response')); run_time_response_study;
    end
    if ~isfile(fullfile(root,'output','continuous_wing','continuous_wing_results.mat'))
        addpath(fullfile(root,'examples','continuous_wing')); run_continuous_wing_study;
    end
    make_doc_figures();

    tex = fullfile(root,'manual','ContourRoots_manual.tex');
    texPath = 'PATH="$PATH:/Library/TeX/texbin:/usr/local/bin:/opt/homebrew/bin"; ';
    if ispc, texPath = ''; end
    cmd = sprintf('%slatexmk -pdf -cd -interaction=nonstopmode -halt-on-error "%s"', texPath, tex);
    [status,out] = system(cmd);
    if status ~= 0
        fprintf('%s\n', out);
        error('build_manual:LaTeX','latexmk failed; see manual/ContourRoots_manual.log.');
    end
    log = fileread(fullfile(root,'manual','ContourRoots_manual.log'));
    if contains(log,'undefined') || contains(log,'Citation') && contains(log,'undefined')
        warning('build_manual:Undefined','The manual has undefined references or citations.');
    end
    copyfile(fullfile(root,'manual','ContourRoots_manual.pdf'), ...
        fullfile(root,'docs','ContourRoots_manual.pdf'));
    fprintf('Manual written to %s\n', fullfile(root,'docs','ContourRoots_manual.pdf'));
end
