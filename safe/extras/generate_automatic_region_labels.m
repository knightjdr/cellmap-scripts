function layout = generate_automatic_region_labels(layout)

SGN = size(layout.pval, 3);


for sgn = 1 : SGN

    ucolors = unique(layout.groupColor(:,sgn));
    ucolors(ucolors==1) = [];
    ucolors_label = cell(size(ucolors));
    words_to_exclude = {'to','or','and','the','a','an','via','of','from','into','in','by','process'};
    for i = 1 : length(ucolors)
        inds = find(layout.groupColor(:,sgn) == ucolors(i));
        tmp = [];
        for j = 1 : length(inds)
            tmp = [tmp regexp(layout.group_names{inds(j)},' ','split')];
        end
        tmp = tmp(~ismember(tmp, words_to_exclude));
        tmp = cellfun(@(x) regexprep(x,',',''), tmp,'UniformOutput',0);
        [ulabels,~,ic] = unique(tmp);
        ulabels_n = histc(ic,1:numel(ulabels));
        [~,ix] = sort(ulabels_n,'descend');
        
        % Join labels into a single string
        %         ucolors_label{i} = strjoin(ulabels(ix(1:min(length(ix),5))),' ');     % Only for Matlab 2013 and higher.
        c = ulabels(ix(1:min(length(ix),5)));
        j = cell(2,numel(c));
        j(1,:) = reshape(c,1,numel(c));
        j(2,1:numel(c)-1) = {' '};
        ucolors_label{i} = [j{:}];
        
    end
    
    layout.regionId{sgn} = ucolors;
    layout.regionName{sgn} = ucolors_label;
    
end