function build_release()
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

    root = fileparts(fileparts(mfilename('fullpath')));
    dist = fullfile(root,'dist');
    if ~isfolder(dist), mkdir(dist); end
    manual = fullfile(root,'docs','ContourRoots_manual.pdf');
    assert(isfile(manual),'build_release:Manual', ...
        'docs/ContourRoots_manual.pdf is missing; run "buildtool manual" first.');

    % Stage a clean copy of the released files.
    stage = fullfile(tempname,'ContourRoots');
    mkdir(stage);
    cleaner = onCleanup(@() rmdir(fileparts(stage),'s'));
    items = {'matlab','examples','docs','setup_contourroots.m','README.md', ...
        'LICENSE','CITATION.cff','CHANGELOG.md','CONTRIBUTING.md','ROADMAP.md'};
    for k = 1:numel(items)
        copyfile(fullfile(root,items{k}), fullfile(stage,items{k}));
    end
    dev = fullfile(stage,'docs','development');
    if isfolder(dev), rmdir(dev,'s'); end
    remove_matching(stage,{'.DS_Store','*.asv'});

    % ZIP: files at the top level, so that unzip(...,"ContourRoots") and the
    % operating-system "extract" both give ContourRoots/setup_contourroots.m.
    zipFile = fullfile(dist,'ContourRoots.zip');
    if isfile(zipFile), delete(zipFile); end
    zip(zipFile, items, stage);
    fprintf('Wrote %s\n', zipFile);

    % Add-On installer.
    if verLessThan('matlab','9.14')   % R2023a
        warning('build_release:mltbx','MATLAB R2023a or later is needed to build the .mltbx.');
        return
    end
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
    opts.ToolboxMatlabPath = {fullfile(stage,'matlab'), fullfile(stage,'matlab','delay')};
    opts.OutputFile = fullfile(dist,'ContourRoots.mltbx');
    matlab.addons.toolbox.packageToolbox(opts);
    fprintf('Wrote %s\n', opts.OutputFile);
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
