function results = run_tests(varargin)
%RUN_TESTS Run the ContourRoots test suite.
%   RUN_TESTS runs every test below this folder (unit, regression and
%   documentation tests) and throws an error if any test fails. Tests that
%   need an unavailable optional toolbox are reported as "Incomplete"
%   (skipped), not as failures.
%
%   RUN_TESTS('unit') or RUN_TESTS('regression','docs') runs selected
%   subfolders only.
%
%   RESULTS = RUN_TESTS(...) also returns the matlab.unittest results.

    here = fileparts(mfilename('fullpath'));
    if nargin == 0
        folders = {'unit','regression','docs'};
    else
        folders = varargin;
    end
    suites = cellfun(@(f) testsuite(fullfile(here,f)), folders, 'UniformOutput', false);
    suite = [suites{:}];
    results = run(suite);
    disp(table(results));
    skipped = sum([results.Incomplete]);
    if skipped > 0
        fprintf('%d test(s) skipped because an optional toolbox is not available.\n', skipped);
    end
    failed = results([results.Failed]);
    if ~isempty(failed)
        error('run_tests:Failed','%d test(s) failed: %s', numel(failed), ...
            strjoin({failed.Name}, ', '));
    end
    fprintf('All %d ContourRoots tests passed.\n', numel(results) - skipped);
    if nargout == 0, clear results; end
end
