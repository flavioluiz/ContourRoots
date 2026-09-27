function [v,i]=response_get_kernel(m,t,power,o,data)
% Same convolution path for freshly inverted and explicitly reused kernels.
    if isempty(data)
        [v,i]=response_kernels(m,t,power,o);
    else
        if power==1
            v=data.step; i=data.stepInfo;
        else
            v=data.ramp; i=data.rampInfo;
        end
        i.evaluations=0; % No evaluation in THIS call; history is preparation.
    end
end
