function tests=test_wing_hybrid, tests=functiontests(localfunctions); end
function setupOnce(tc)
    root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
    tc.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,'matlab')));
    tc.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,'examples','models')));
end
function testElementJoinAndInterface(tc)
    for kind={'dry','goland'}
        w=wing_model(kind{1}); half=w.strips; half.length=w.L/2;
        for s=[0 2+40i 5+200i]
            [A,ia]=wing_hybrid(s,120,half,w.scale); [H,ih]=wing_hybrid(s,120,w.strips,w.scale);
            [J,ij]=wing_hybrid_join(A,A); tc.verifyTrue(ia.available&&ih.available&&ij.available);
            tc.verifyLessThan(norm(J-H,'fro')/norm(H,'fro'),1e-10);
            u=[.001;.002;.003;1;2;3]; qI=ij.interfaceMap*u;
            fLB=A(1:3,:)*[qI;u(4:6)]; left=A*[u(1:3);-fLB]; right=A*[qI;u(4:6)];
            tc.verifyLessThan(norm(left(4:6)-qI),1e-10);
            tc.verifyLessThan(norm([left(1:3);right(4:6)]-J*u)/max(1,norm(J*u)),1e-10);
        end
    end
end
function testTransferAndStructuredEquivalence(tc)
    for kind={'dry','goland','tapered'}
        w=wing_model(kind{1},4); M=wing_cdyn(120,w,2);
        for s=[0 3+70i 5+300i]
            A=wing_transfer_matrix(s,120,w,'implicit'); B=wing_transfer_matrix(s,120,w,'hybrid');
            tc.verifyLessThan(norm(A-B,'fro')/norm(A,'fro'),1e-9);
            tc.verifyLessThan(norm(A-ceval(M,s),'fro')/norm(A,'fro'),1e-9);
            for j=1:3, for i=1:3
                tc.verifyLessThan(abs(A(i,j)-wing_transfer(s,120,w,i,j)),1e-10*max(1,abs(A(i,j))));
            end, end
        end
    end
end
function testStaticAndUnavailablePivot(tc)
    w=wing_model('dry'); [H,info]=wing_hybrid(0,0,w.strips,w.scale); tc.verifyTrue(info.available);
    e=w.strips; L=w.L; ref=[L^3/(3*e.EI) L^2/(2*e.EI) 0;L^2/(2*e.EI) L/e.EI 0;0 0 L/e.GJ];
    tc.verifyEqual(H(4:6,4:6),ref,'RelTol',1e-12,'AbsTol',1e-15);
    omega=1.875104068711961^2*sqrt(e.EI/e.mu)/L^2;
    [H,info]=wing_hybrid(1i*omega,0,e,w.scale);
    tc.verifyFalse(info.available); tc.verifyTrue(all(isnan(H(:))));
    tc.verifyTrue(all(isfinite(info.implicit(:))));
end
function testHierarchicalVersusSequential(tc)
    w=wing_model('goland'); e=w.strips; e.length=e.length/16;
    [A,info]=wing_hybrid(5+2e4i,120,e,w.scale); tc.verifyTrue(info.available);
    for n=[1 2 3 7 16]
        sequential=A;
        for k=2:n, [sequential,info]=wing_hybrid_join(sequential,A); tc.verifyTrue(info.available); end
        fast=wing_hybrid_repeat(A,n);
        tc.verifyLessThan(norm(fast-sequential,'fro')/norm(sequential,'fro'),1e-10);
    end
end
