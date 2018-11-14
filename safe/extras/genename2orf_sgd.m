function S = genename2orf_sgd(genenames)

load sgd_141010;

genenames = lower(strtrim(genenames));
[genenames,ia,ic] = unique(genenames);

s = cell(size(genenames));

% Check if any of the genenames are already in ORF format
expr = 'Y[A-P][RL][0-9]{3}[CW](-[ABC])*';
inds = find(~cellfun(@isempty, regexpi(genenames, expr)));
s(inds) = upper(genenames(inds));

% Now, match the rest to the SGD genenames
inds = find(cellfun(@isempty,s));
[~,ind1,ind2] = intersect(genenames(inds), sgd.genenames);
s(inds(ind1)) = sgd.orfs(ind2);

% Print out things that rest that couldn't be matched
inds = find(cellfun(@isempty, s));
if ~isempty(inds)
    warning('These genenames could not be matched:')
    for i = 1 : length(inds)
        fprintf('%s\n', genenames{inds(i)});
    end
    fprintf('\n');
end

s(inds) = genenames(inds);

S = s(ic);