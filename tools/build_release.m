function report=build_release(outputDir)
%BUILD_RELEASE Build dist/ContourRoots.zip and dist/ContourRoots.mltbx.
%   The release contains what users need: the toolbox (matlab/), examples,
%   documentation (including the PDF manual), setup_contourroots.m, README,
%   LICENSE, CITATION.cff, CHANGELOG.md, CONTRIBUTING.md and ROADMAP.md
%   (all linked from the README). Tests, tools, manual sources and
%   development notes are not included.
%
%   The file names do not contain the version on purpose: the README links
%   to https://github.com/<owner>/ContourRoots/releases/latest/download/<name>.
%
%   The .mltbx installer needs MATLAB R2023a or later to be built.
%   BUILD_RELEASE(OUTPUTDIR) uses a separate output directory. REPORT gives
%   per-stage timings. This function ONLY packages; BUILDTOOL RELEASE runs
%   the validation/study/manual gates first. No files are published online.

    root = fileparts(fileparts(mfilename('fullpath')));
    if nargin<1, outputDir=fullfile(root,'dist'); end
    dist=char(outputDir);
    if ~isfolder(dist), mkdir(dist); end
    % Normalize before packageToolbox receives a potentially relative path.
    [ok,attributes]=fileattrib(dist); assert(ok,'build_release:Output','Cannot access output directory.');
    dist=attributes.Name;
    totalTimer=tic; report=struct('matlab',version,'version',contourroots_version_of(root), ...
        'scope','Packaging only; buildtool release supplies validation gates.','stages',struct);
    manual = fullfile(root,'docs','ContourRoots_manual.pdf');
    assert(isfile(manual),'build_release:Manual', ...
        'docs/ContourRoots_manual.pdf is missing; run "buildtool manual" first.');

    % Stage a clean copy of the released files.
    fprintf('Packaging: staging released files...\n'); timer=tic;
    stage = fullfile(tempname,'ContourRoots');
    mkdir(stage);
    cleaner = onCleanup(@() rmdir(fileparts(stage),'s'));
    items = {'matlab','examples','docs','setup_contourroots.m','README.md', ...
        'LICENSE','THIRD_PARTY_NOTICES.md','CITATION.cff','CHANGELOG.md','CONTRIBUTING.md','ROADMAP.md'};
    for k = 1:numel(items)
        copyfile(fullfile(root,items{k}), fullfile(stage,items{k}));
    end
    dev = fullfile(stage,'docs','development');
    if isfolder(dev), rmdir(dev,'s'); end
    remove_matching(stage,{'.DS_Store','*.asv'});
    elapsed('staging',timer);

    % ZIP: files at the top level, so that unzip(...,"ContourRoots") and the
    % operating-system "extract" both give ContourRoots/setup_contourroots.m.
    % Prepare both archives before replacing existing release artifacts.
    fprintf('Packaging: ZIP archive...\n'); timer=tic;
    zipFile = fullfile(fileparts(stage),'ContourRoots.zip');
    zip(zipFile, items, stage);
    elapsed('zip',timer);

    % Add-On installer.
    if verLessThan('matlab','9.14')   % R2023a
        warning('build_release:mltbx','MATLAB R2023a or later is needed to build the .mltbx.');
        movefile(zipFile,fullfile(dist,'ContourRoots.zip'),'f');
        finish();
        return
    end
    fprintf('Packaging: MATLAB installer options...\n'); timer=tic;
    uuid = '6f1b2e7c-3c1a-4a55-9b7e-0c4b8f2d9a31';   % fixed across versions
    opts = matlab.addons.toolbox.ToolboxOptions(stage, uuid);
    opts.ToolboxName = 'ContourRoots';
    opts.ToolboxVersion = contourroots_version_of(root);
    opts.AuthorName = 'Flávio Luiz Cardoso-Ribeiro';
    opts.Summary = 'Poles, zeros and roots of nonrational functions in a rectangle of the complex plane.';
    opts.Description = ['Finds roots of scalar analytic functions and poles and zeros of ' ...
        'nonrational transfer functions (time delays, PDE models) with the argument ' ...
        'principle, without Padé approximation or modal truncation. Includes ' ...
        'time-delay stability tools, examples and a PDF manual.'];
    opts.ToolboxMatlabPath = {fullfile(stage,'matlab'),fullfile(stage,'matlab','delay')};
    opts.OutputFile = fullfile(fileparts(stage),'ContourRoots.mltbx');
    report.fileCount=numel(opts.ToolboxFiles); elapsed('options',timer);
    fprintf('Packaging: MATLAB installer archive...\n'); timer=tic;
    matlab.addons.toolbox.packageToolbox(opts);
    elapsed('mltbx',timer);
    timer=tic;
    movefile(zipFile,fullfile(dist,'ContourRoots.zip'),'f');
    movefile(opts.OutputFile,fullfile(dist,'ContourRoots.mltbx'),'f');
    elapsed('copyArtifacts',timer); finish();
    function elapsed(name,timer)
        report.stages.(name)=toc(timer);
        fprintf('Packaging: %s completed in %.3f s.\n',name,report.stages.(name));
    end
    function finish()
        report.totalSeconds=toc(totalTimer);
        fid=fopen(fullfile(dist,'packaging_timings.json'),'w');
        if fid<0
            warning('build_release:Timings','Archives were built, but the timing report cannot be written.');
        else
            closer=onCleanup(@() fclose(fid)); fprintf(fid,'%s\n',jsonencode(report,'PrettyPrint',true));
        end
        fprintf('Packages written to %s (%.3f s; validation not included).\n',dist,report.totalSeconds);
    end
end

function v = contourroots_version_of(root)
    text = fileread(fullfile(root,'matlab','contourroots_version.m'));
    v = char(regexp(text,'''(\d+\.\d+\.\d+)''','tokens','once'));
end

function remove_matching(folder,patterns)
    for p = patterns
        d = dir(fullfile(folder,'**',p{1}));
        for k = 1:numel(d), delete(fullfile(d(k).folder,d(k).name)); end
    end
end
