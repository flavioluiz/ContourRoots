function c = spectrum_count(f,box,opt)
% Numerical argument principle; convergence is evidence, not an interval proof.
    previous = NaN; stable = 0; oldZ=[]; oldV=[];
    c = struct('count',NaN,'ok',false,'center',NaN,'samples',0);
    for level = 0:opt.ContourRefinements
        n = opt.ContourPoints*2^level;
        t = (0:n-1)/n;
        z = [box(1)+diff(box(1:2))*t+1i*box(3), ...
            box(2)+1i*(box(3)+diff(box(3:4))*t), ...
            box(2)-diff(box(1:2))*t+1i*box(4), ...
            box(1)+1i*(box(4)-diff(box(3:4))*t)];
        % Dyadic refinement retains every previous node. Reuse only exact
        % coordinate matches, within this contour call (never persistently).
        if ~isempty(oldZ) && isequal(z(1:2:end),oldZ)
            v=zeros(size(z)); v(1:2:end)=oldV;
            v(2:2:end)=spectrum_values(f,z(2:2:end));
        else
            v = spectrum_values(f,z);
        end
        oldZ=z; oldV=v;
        c.samples = numel(z);
        if any(~isfinite(v) | v==0), return; end
        phase = angle(v); logmag = log(abs(v));
        dphase = angle(exp(1i*(phase([2:end 1])-phase)));
        dmag = logmag([2:end 1])-logmag;
        winding = sum(dphase)/(2*pi);
        count = round(winding);
        resolved = max(abs(dphase)) < pi/3 && max(abs(dmag)) < 2;
        if count==previous && resolved && abs(winding-count)<1e-7
            stable = stable+1;
        else
            stable = 0;
        end
        previous = count;
        if stable >= 2
            c.count = count; c.ok = true;
            % First logarithmic moment, used only as a Newton seed.
            midpoint = (z+z([2:end 1]))/2;
            if count~=0
                c.center = sum(midpoint.*(dmag+1i*dphase))/(2i*pi*count);
            end
            return
        end
    end
end
