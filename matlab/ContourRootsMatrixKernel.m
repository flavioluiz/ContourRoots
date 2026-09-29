classdef (Sealed) ContourRootsMatrixKernel
%CONTOURROOTSMATRIXKERNEL Immutable numeric MIMO kernel bank from CKERNEL.
%   Includes every channel, even if its input is zero in the first run.
%   SharedGrid=true opts into common fft/dehoog preparation (see CKERNEL).
%   No model handles or shared evaluation cache are retained. Scalar
%   ContourRootsKernel snapshots remain unchanged and load as before.
    properties (Access=private)
        Data
    end
    properties (Dependent, SetAccess=private)
        Time
        Interpolation
        Size
        Info
    end
    methods
        function obj=ContourRootsMatrixKernel(M,t,varargin)
            obj.Data=matrix_response_prepare(M,t,varargin);
        end
        function v=get.Time(obj), v=obj.Data.time; end
        function v=get.Interpolation(obj), v=obj.Data.interpolation; end
        function v=get.Size(obj), v=obj.Data.size; end
        function v=get.Info(obj), v=obj.Data.info; end
        function [y,t,info]=clsim(obj,u,t,varargin)
            if ~isscalar(obj), error('ContourRoots:KernelModel','Use one matrix kernel bank.'); end
            [y,t,info]=matrix_response_run('lsim',[],u,t,varargin,nargout==0,obj.Data);
        end
    end
end
