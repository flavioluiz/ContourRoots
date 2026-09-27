function M=wing_cdyn(U,wing,substeps)
%WING_CDYN Fixed-size structured tip model, useful for inspecting assembly.
%   SUBSTEPS fixes the number of exact numerical pieces per physical strip.
%   Select enough for the frequency range and check subdivision convergence.
%   For unbounded/adaptive inversion grids use CMIMO with
%   WING_TRANSFER_MATRIX: its internal solve size adapts while G stays 3x3.
    if nargin<3, substeps=1; end
    validateattributes(substeps,{'numeric'},{'scalar','integer','positive','finite'});
    n=6*(numel(wing.strips)*substeps+1); B=zeros(n,3); C=zeros(3,n);
    B(end-2:end,:)=diag(1./wing.scale(4:6));
    C(:,end-5:end-3)=diag(wing.scale(1:3));
    M=cdyn(@(s) wing_matrix(s,U,wing,substeps),B,C,zeros(3), ...
        'Dimensions',[n 3 3],'InputNames',{'force','moment','torque'}, ...
        'OutputNames',{'deflection','slope','twist'},'InputUnits',{'N','N m','N m'}, ...
        'OutputUnits',{'m','rad','rad'},'Feedthrough',zeros(3));
end
