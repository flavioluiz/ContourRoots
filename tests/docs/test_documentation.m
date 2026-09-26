function tests = test_documentation
%TEST_DOCUMENTATION Run the README, the Markdown docs and the short examples.
%   Every ```matlab block of README.md and docs/**/*.md is executed (except
%   blocks marked <!-- no-test -->), as well as the quick-start scripts.
    tests = functiontests(localfunctions);
end

function setupOnce(tc)
    repo = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    import matlab.unittest.fixtures.PathFixture
    tc.applyFixture(PathFixture(fullfile(repo,'matlab')));
    tc.applyFixture(PathFixture(fullfile(repo,'matlab','delay')));
    tc.applyFixture(PathFixture(fullfile(repo,'examples','models')));
    tc.applyFixture(PathFixture(fullfile(repo,'tools')));
    tc.TestData.repo = repo;
    old = get(groot,'defaultFigureVisible');
    set(groot,'defaultFigureVisible','off');
    tc.addTeardown(@() set(groot,'defaultFigureVisible',old));
end

function testReadmeAndDocs(tc)
    repo = tc.TestData.repo;
    files = [{fullfile(repo,'README.md')}; markdown_files(fullfile(repo,'docs'))];
    report = check_doc_snippets(files);
    tc.verifyFalse(any(report.Status == "failed"));
    tc.verifyGreaterThan(sum(report.Status == "passed"), 0);
end

function testQuickstartScripts(tc)
    repo = tc.TestData.repo;
    scripts = [dir(fullfile(repo,'examples','quickstart','*.m')); ...
               dir(fullfile(repo,'examples','coupled_beam','beam_direct_nd.m'))];
    tc.assertNotEmpty(scripts);
    for k = 1:numel(scripts)
        file = fullfile(scripts(k).folder,scripts(k).name);
        run_isolated(file);
        close(findall(groot,'Type','figure'));
    end
end

function files = markdown_files(folder)
    d = dir(fullfile(folder,'**','*.md'));
    d = d(~contains({d.folder},'development'));   % internal planning notes
    files = fullfile({d.folder},{d.name}).';
end

function run_isolated(file)
    % Each script gets a fresh workspace (this function's).
    evalc('run(file)');
end
