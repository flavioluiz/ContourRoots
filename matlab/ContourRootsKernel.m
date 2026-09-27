classdef (Sealed) ContourRootsKernel
%CONTOURROOTSKERNEL Read-only numerical snapshot created with CKERNEL.
%   Time, Interpolation and Info are inspectable, not assignable. The
%   originating evaluator is not retained: changing a captured parameter
%   cannot silently change the prepared system. Prepare a new K instead.
    properties (Access=private)
        Data
    end
    properties (Dependent, SetAccess=private)
        Time
        Interpolation
        Info
    end
    methods
        function obj=ContourRootsKernel(G,t,varargin)
            obj.Data=response_prepare(G,t,varargin);
        end
        function value=get.Time(obj), value=obj.Data.time; end
        function value=get.Interpolation(obj), value=obj.Data.options.Interpolation; end
        function value=get.Info(obj), value=obj.Data.info; end
        function [y,tOut,info]=clsim(obj,u,t,varargin)
            if ~isscalar(obj)
                error('ContourRoots:KernelModel','Use one prepared kernel object per simulation.');
            end
            [y,tOut,info]=response_reuse(obj.Data,u,t,varargin,nargout==0);
        end
    end
end
