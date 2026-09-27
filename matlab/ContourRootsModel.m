classdef (Sealed) ContourRootsModel
%CONTOURROOTSMODEL Explicit matrix transfer or structured dynamic model.
%   Construct with CMIMO or CDYN; evaluate with CEVAL. Dimensions are never
%   inferred by probing a handle. This value object holds evaluators, not a
%   frozen numerical snapshot: use CKERNEL to freeze a time-response model.
    properties (SetAccess=private)
        Representation
        Size
        Dimensions
        Factors
        Options
    end
    methods
        function obj=ContourRootsModel(kind,varargin)
            [obj.Representation,obj.Size,obj.Dimensions,obj.Factors,obj.Options]=matrix_model(kind,varargin{:});
        end
    end
end
