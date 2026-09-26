function [y,i]=response_invert(f,t,o,bound,initial)
    i=response_info(t,o); y=NaN(size(t)); pos=t>0;
    if any(pos)
        if strcmp(o.Method,'fft')
            [v,j]=response_fft(f,t(pos),o,bound);
        elseif strcmp(o.Method,'dehoog')
            [v,j]=response_dehoog(f,t(pos),o,bound);
        else
            [v,j]=response_quadrature(f,t(pos),o,bound);
        end
        y(pos)=v; i=j; i.time=t;
        for field={'errorEstimate','resolvedMask','periodError','bandwidthError','shiftError'}
            key=field{1}; q=i.(key);
            if strcmp(key,'resolvedMask'), out=false(size(t)); else, out=inf(size(t)); end
            out(pos)=q; i.(key)=out;
        end
    end
    y(~pos)=initial;
    if isfinite(initial)
        for field={'errorEstimate','periodError','bandwidthError','shiftError'}
            i.(field{1})(~pos)=0;
        end
        i.resolvedMask(~pos)=true;
    end
    i.converged=all(i.resolvedMask);
    if i.converged, i.status='converged'; else, i.status='unresolved'; end
    if any(~pos)&&~isfinite(initial)
        i.warnings{end+1}='Unknown value at t=0; provide an analytical InitialValue or request positive times.';
        i.stopReason='Unknown zero-time limit; inspect positive-time convergence separately.';
    end
end
