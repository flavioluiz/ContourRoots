function M=wing_response_model(U,wing,assembly,channels)
%WING_RESPONSE_MODEL Select tip channels and convert to mm/mrad per kN/kNm.
%   Defaults to force/torque -> deflection/twist, channels [1 3]. The full
%   6-state strip is retained; only external inputs/outputs are selected.
%   Channels [1 2 3] expose all nine transfers (costlier kernel preparation).
    if nargin<3, assembly='implicit'; end
    if nargin<4, channels=[1 3]; end
    validateattributes(channels,{'numeric'},{'vector','nonempty','integer','>=',1,'<=',3});
    if numel(unique(channels))~=numel(channels), error('wing:Channels','Use distinct channel indices.'); end
    inputs={'force','moment','torque'}; outputs={'deflection','slope','twist'};
    iu={'kN','kN m','kN m'}; ou={'mm','mrad','mrad'};
    M=cmimo(@evaluate,'Size',[numel(channels) numel(channels)], ...
        'InputNames',inputs(channels),'OutputNames',outputs(channels), ...
        'InputUnits',iu(channels),'OutputUnits',ou(channels),'Feedthrough',zeros(numel(channels)));
    function G=evaluate(s)
        allChannels=wing_transfer_matrix(s,U,wing,assembly);
        G=1e6*allChannels(channels,channels);
    end
end
