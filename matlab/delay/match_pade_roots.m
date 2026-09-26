function comparison = match_pade_roots(exactRoots,padeRoots,nDominant)
%MATCH_PADE_ROOTS Match each dominant exact root to its nearest Pade root.

    if nargin < 3, nDominant = numel(exactRoots); end
    nDominant = min(nDominant,numel(exactRoots));
    exactRoots = exactRoots(1:nDominant);
    matched = NaN(size(exactRoots));
    errors = Inf(size(exactRoots));
    available = true(size(padeRoots));
    for k = 1:numel(exactRoots)
        idx = find(available);
        if isempty(idx), break; end
        [errors(k),j] = min(abs(padeRoots(idx)-exactRoots(k)));
        matched(k) = padeRoots(idx(j));
        available(idx(j)) = false;
    end
    comparison = table(exactRoots,matched,errors, ...
        'VariableNames',{'ExactRoot','PadeRoot','AbsoluteError'});
end
