function [locations,info] = cr_search(caller,mode,model,region,args,warnIfIncomplete)
% Shared entry point of CROOTS, CPOLES and CZEROS.
%   Converts convenience inputs, calls COMPLEX_SPECTRUM and, when requested,
%   warns if the result is not numerically complete. The numerical work is
%   done entirely by COMPLEX_SPECTRUM.

    if nargin < 2 || isempty(region)
        error('ContourRoots:Region', ['%s needs a search rectangle, e.g. ' ...
            '%s(F,[-5 1 -10 10]) for -5<=Re(s)<=1, -10<=Im(s)<=10.'], ...
            upper(caller), caller);
    end
    if ~(isnumeric(region) && numel(region)==4)
        error('ContourRoots:Region', ['The search region must be ' ...
            '[xmin xmax ymin ymax], e.g. [-5 1 -10 10].']);
    end
    % Name-value pairs owned by the wrappers are removed before forwarding.
    [warnOption,args] = take_option(args,'Warn',true);
    if any(strcmpi(args(1:2:end),'Mode'))
        error('ContourRoots:Mode','%s sets Mode itself; use CROOTS, CPOLES or CZEROS.', ...
            upper(caller));
    end
    if isnumeric(model)
        % A coefficient vector (descending powers), as accepted by ROOTS.
        validateattributes(model,{'numeric'},{'vector','nonempty','finite'},caller,'F');
        model = ndpair(model,1);
    end
    [locations,info] = complex_spectrum(model,region,args{:},'Mode',mode);
    if warnIfIncomplete && warnOption && ~info.complete
        switch info.status
            case 'exploratory'
                msg = ['The search is exploratory: completeness is not checked for an ' ...
                    'opaque function handle. If the function is analytic in the region, ' ...
                    'pass ''AssumeAnalytic'',true (see "help %s").'];
            otherwise
                msg = ['The search is not numerically complete (status "%s"). Some ' ...
                    'roots may be missing. Call [r,info] = %s(...) and inspect ' ...
                    'info.unresolvedBoxes; try a slightly shifted or larger rectangle.'];
        end
        if strcmp(info.status,'exploratory')
            warning('ContourRoots:Exploratory',msg,caller);
        else
            warning('ContourRoots:Incomplete',msg,info.status,caller);
        end
    end
end

function [value,args] = take_option(args,name,default)
    value = default;
    k = find(strcmpi(args(1:2:end),name),1,'last');
    if ~isempty(k)
        value = args{2*k};
        args(2*k-1:2*k) = [];
        validateattributes(value,{'logical','numeric'},{'scalar'},'',name);
        value = logical(value);
    end
end
