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
               dir(fullfile(repo,'examples','time_response','first_step.m')); ...
               dir(fullfile(repo,'examples','time_response','arbitrary_input.m')); ...
               dir(fullfile(repo,'examples','time_response','delays_and_diffusion.m')); ...
               dir(fullfile(repo,'examples','time_response','unstable_response.m')); ...
               dir(fullfile(repo,'examples','time_response','reuse_kernels.m')); ...
               dir(fullfile(repo,'examples','coupled_beam','beam_direct_nd.m')); ...
               dir(fullfile(repo,'examples','aeroelasticity','flutter_quickstart.m'))];
    tc.assertNotEmpty(scripts);
    before = findall(groot,'Type','figure');
    for k = 1:numel(scripts)
        file = fullfile(scripts(k).folder,scripts(k).name);
        run_isolated(file);
        close(setdiff(findall(groot,'Type','figure'),before));
    end
end

function testCodeBlocksMatchExampleFiles(tc)
    % A block preceded by <!-- file: path --> must equal that file.
    repo = tc.TestData.repo;
    files = markdown_files(fullfile(repo,'docs'));
    checked = 0;
    for i = 1:numel(files)
        lines = splitlines(string(fileread(files{i})));
        for k = find(startsWith(strtrim(lines),'<!-- file:')).'
            target = strtrim(extractBetween(lines(k),'<!-- file:','-->'));
            tc.assertEqual(strtrim(lines(k+1)),"```matlab", ...
                sprintf('%s: marker not followed by a matlab block',files{i}));
            stop = k + 1 + find(strtrim(lines(k+2:end)) == "```",1);
            block = strjoin(lines(k+2:stop-1),newline);
            source = strip(string(fileread(fullfile(repo,target))),'right');
            tc.verifyEqual(block,source,sprintf('%s is out of sync with %s', ...
                files{i},target));
            checked = checked + 1;
        end
    end
    tc.verifyGreaterThan(checked,0);
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
