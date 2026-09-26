function sigma=response_line(o,bound,period)
    if ~isempty(o.Abscissa), sigma=o.Abscissa;
    elseif o.AssumeStable, sigma=0;
    else
        if isempty(bound)||~isfinite(bound), bound=0; end
        sigma=max(0,bound)+max(12,-log(min(o.AbsTol,o.RelTol)*.01))/period;
    end
end
