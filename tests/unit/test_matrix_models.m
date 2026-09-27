function tests=test_matrix_models, tests=functiontests(localfunctions); end
function setupOnce(tc)
    root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
    tc.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,'matlab')));
end
function testNoProbeAndBlockSolve(tc)
    calls=0; M=cdyn(@H,eye(2),[1 2],zeros(1,2),'Dimensions',[2 1 2]);
    tc.verifyEqual(calls,0);
    [G,i]=ceval(M,[1+2i 3-1i]);
    tc.verifySize(G,[1 2 2]); tc.verifyEqual(calls,2);
    tc.verifyEqual(i.factorizations,2); tc.verifyEqual(i.rhsColumns,4);
    tc.verifyEqual(G(:,:,1),[1/((1+2i)+1) 2/((1+2i)+2)],'AbsTol',1e-14);
    function v=H(s), calls=calls+1; v=diag([s+1 s+2]); end
end
function testBatchAndShapes(tc)
    A=cmimo(@(s) reshape(1./(s+1),1,1,[]),'Size',[1 1],'Batched',true);
    B=cmimo(@(s) 1/(s+1),'Size',[1 1]);
    tc.verifyEqual(ceval(A,[1 2 3]),ceval(B,[1 2 3]));
    tc.verifyError(@() ceval(cmimo(@(s) eye(2),'Size',[1 2]),1),'ContourRoots:MatrixShape');
    tc.verifyError(@() cmimo(@(s) 1/s),'ContourRoots:MatrixShape');
    tc.verifyError(@() cdyn(eye(2),ones(3,1),ones(1,2),0),'ContourRoots:MatrixShape');
    tc.verifyError(@() cmimo({[1 2]}),'ContourRoots:MatrixRepresentation');
    tc.verifyEqual(ceval(cmimo([1 2;3 4]),1),[1 2;3 4]);
    tc.verifyEqual(ceval(cdyn(speye(2),eye(2),eye(2),zeros(2)),[1 2]),repmat(eye(2),1,1,2));
end
function testDomainBudgetAndChannels(tc)
    M=cmimo({ndpair(1,[1 1]),0;@(s) exp(-s),2},'DomainCheck',@(s) real(s)>0);
    tc.verifyEqual(ceval(M,1),[.5 0;exp(-1) 2],'AbsTol',1e-14);
    f=cchannel(M,1,1); tc.verifyEqual(f([1 2]),[.5 1/3]);
    tc.verifyError(@() ceval(M,-1),'ContourRoots:MatrixDomain');
    tc.verifyError(@() ceval(M,[1 2],'MaxEvaluations',1),'ContourRoots:MatrixBudget');
    tc.verifyError(@() ceval(M,1,'MaxMemoryMB',1e-6),'ContourRoots:MatrixBudget');
    tc.verifyError(@() ceval(cmimo(@(s) NaN,'Size',[1 1]),1),'ContourRoots:MatrixShape');
end
function testOptionalSymbolic(tc)
    tc.assumeTrue(~isempty(ver('symbolic')));
    s=sym('s'); M=cmimo([1/(s+1) exp(-s);0 2]);
    tc.verifyEqual(ceval(M,1),[.5 exp(-1);0 2],'AbsTol',1e-14);
    tc.verifyEqual(ceval(cdyn(s+1,1,1,0),1),.5);
end
function testOptionalLTIAndDescriptor(tc)
    tc.assumeTrue(~isempty(ver('control')));
    sys=ss(diag([-1 -2]),eye(2),[1 2],zeros(1,2)); M=cmimo(sys);
    tc.verifyEqual(ceval(M,2),evalfr(sys,2),'AbsTol',1e-14);
    tc.verifyError(@() cmimo(dss(eye(2),eye(2),eye(2),zeros(2),diag([1 0]))),'ContourRoots:MatrixRepresentation');
end
function testBatchedStructuredAndFreshParameters(tc)
    s=[1+2i 3+4i];
    M=cdyn(@(z) reshape(z+1,1,1,[]),1,2,0,'Dimensions',[1 1 1],'Batched',true);
    [g,i]=ceval(M,s); tc.verifyEqual(reshape(g,1,[]),2./(s+1));
    tc.verifyEqual(i.factorEvaluations,1); tc.verifyEqual(i.factorizations,2);
    rate=1; N=cmimo(@value,'Size',[1 1]); a=ceval(N,1); rate=2; b=ceval(N,1);
    tc.verifyEqual(a,.5); tc.verifyEqual(b,1/3);
    function v=value(z), v=1/(z+rate); end
end
