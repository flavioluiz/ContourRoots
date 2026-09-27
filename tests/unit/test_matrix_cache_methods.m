function tests=test_matrix_cache_methods, tests=functiontests(localfunctions); end
function setupOnce(tc)
    root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
    tc.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,'matlab')));
end
function testOpaqueNonFFTRoutes(tc)
    % The channels representation bypasses the shared matrix cache. Exercise
    % the opaque representation explicitly with both non-FFT inversions.
    M=cmimo(@(s) [1/(s+1)^2 .5/(s+2)^2;0 1/(s+3)^2],'Size',[2 2]);
    t=[0;.2;.7]; ref=zeros(3,2,2);
    ref(:,1,1)=1-(1+t).*exp(-t);
    ref(:,1,2)=.5*(1-(1+2*t).*exp(-2*t))/4;
    ref(:,2,2)=(1-(1+3*t).*exp(-3*t))/9;
    for method={'dehoog','quadrature'}
        [y,~,i]=cstep(M,t,'Method',method{1},'SingularityBound',0);
        tc.verifyTrue(i.converged); tc.verifyEqual(y,ref,'AbsTol',2e-6);
        tc.verifyGreaterThan(i.cacheHits,0);
    end
end
function testBatchedAndScalarEvaluatorsAgree(tc)
    weights=[1 2;3 4]; t=(0:.1:.5).';
    scalar=cmimo(@(s) weights/(s+1),'Size',[2 2]);
    batched=cmimo(@(s) weights.*reshape(1./(s+1),1,1,[]),'Size',[2 2],'Batched',true);
    [a,~,ia]=cstep(scalar,t,'SingularityBound',0);
    [b,~,ib]=cstep(batched,t,'SingularityBound',0);
    tc.verifyTrue(ia.converged&&ib.converged); tc.verifyEqual(a,b,'AbsTol',1e-12);
    tc.verifyEqual(ia.evaluations,ib.evaluations); tc.verifyEqual(ia.cacheHits,ib.cacheHits);
    tc.verifyLessThan(ib.factorEvaluations,ia.factorEvaluations);
end
