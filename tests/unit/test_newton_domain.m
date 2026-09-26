function tests = test_newton_domain
%TEST_NEWTON_DOMAIN The solver never evaluates the function outside the region.
    tests = functiontests(localfunctions);
end

function setupOnce(tc)
    repo = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    tc.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(repo,'matlab')));
end

function testFunctionUndefinedOutsideRegion(tc)
    % sqrt has a cut on the negative real axis; this guard makes any
    % evaluation outside the rectangle an error.
    region = [0.01 3 -2 2];
    F = @(s) guarded(s,region);
    [r,info] = croots(F,region,'AssumeAnalytic',true);
    tc.verifyTrue(info.complete);
    tc.verifyEqual(r,0.25,'AbsTol',1e-8);
    [r,info] = croots(@(s) guarded(s,[0.01 3 -2 2]).*(s-2+1i),region,'AssumeAnalytic',true);
    tc.verifyTrue(info.complete); tc.verifyEqual(numel(r),2);
end

function v = guarded(s,region)
    if any(real(s(:))<region(1) | real(s(:))>region(2) | imag(s(:))<region(3) | imag(s(:))>region(4))
        error('test:OutsideRegion','Evaluated outside the search region.');
    end
    v = sqrt(s) - 0.5;
end
