function layout = minimize_colors(layout)

SGN = size(layout.pval, 3);

layout.labelColor = nan(length(layout.label), SGN);
layout.labelOpacity = nan(length(layout.label), SGN);
layout.labelOpacity01 = nan(length(layout.label), SGN);

for sgn = 1 : SGN

    ucolors = unique(layout.groupColor(:,sgn));
    ucolors(ucolors == 0) = [];
    
    layout.cumOpacityByColor{sgn} = nan(length(layout.label), length(ucolors));
    layout.cumOpacity01ByColor{sgn} = nan(length(layout.label), length(ucolors));
    
    if isempty(ucolors)
        break;
    end

    for i = 1 : length(ucolors)
        inds = find(layout.groupColor(:,sgn) == ucolors(i));
        layout.cumOpacityByColor{sgn}(:,i) = sum(layout.opacity(:,inds,sgn),2);
        layout.cumOpacity01ByColor{sgn}(:,i) = sum(layout.opacity_01(:,inds,sgn),2);
    end

    for i = 1 : length(layout.label)
        [r,~] = tiedrank(layout.cumOpacity01ByColor{sgn}(i,:));
        inds = find(r == max(r));
        [~,ixColor] = max(layout.cumOpacityByColor{sgn}(i,inds));
        
        layout.labelColor(i,sgn) = ucolors(inds(ixColor));
        layout.labelOpacity(i,sgn) = max(layout.opacity(i, layout.groupColor(:,sgn) == ucolors(inds(ixColor)), sgn));
    end

    layout.labelOpacity01(:,sgn) = layout.labelOpacity(:,sgn) > layout.thresholdOpacity;

    if layout.groupMinimize == 1
        
        % If a certain color has less than a minimum number of nodes at opacity > significance, eliminate that color
        ucolors = unique(layout.labelColor(:,sgn));
        ucolors_size = histc(layout.labelColor(layout.labelOpacity01(:,sgn) > 0,sgn), ucolors);
        layout.labelColor(ismember(layout.labelColor(:,sgn), ucolors(ucolors_size < layout.groupEnrichmentMinSize)), sgn) = 0;
        layout.labelColor(layout.labelOpacity01(:,sgn) ~= 1,sgn) = 0;
        
        % Re-number the colors, so that they are consecutive. This also sets no color (0) to 1, if any.
        [ucolors, ~, ib] = unique(layout.labelColor(:,sgn));
        layout.labelColor(:,sgn) = ib;
        
        % Don't forget to update the colors of the attributes themselves too
        [~, Locb] = ismember(layout.groupColor(:,sgn),ucolors); 
        Locb(Locb==0) = 1;
        layout.groupColor(:,sgn) = Locb;
        
        % And also the columns in the summary
        ucolors(ucolors==0) = [];
        layout.cumOpacityByColor{sgn} = layout.cumOpacityByColor{sgn}(:,ucolors);
        layout.cumOpacity01ByColor{sgn} = layout.cumOpacity01ByColor{sgn}(:,ucolors);
        
    elseif layout.groupMinimize == 2    % "Legacy" procedure
        
        inds = find(layout.groupIsTop(:,sgn));
        [layout.labelOpacity(:,sgn),ix] = max(layout.opacity(:,inds,sgn),[],2);
        layout.labelColor(:,sgn) = layout.groupColor(inds(ix),sgn);

        % If a certain color has less than 10 genes at opacity > significance, ignore that color
        ucolors = unique(layout.labelColor(:,sgn));
        ucolors_size = histc(layout.labelColor(layout.labelOpacity(:,sgn) > layout.thresholdOpacity, sgn), ucolors);
        layout.labelColor(ismember(layout.labelColor(:,sgn), ucolors(ucolors_size < 10)),sgn) = 0;

        % Minimize the colors, such that they are consecutive [this also sets
        % no color (0) = 1].
        [ucolors, ~, ib] = unique(layout.labelColor(:,sgn));
        ucolors_min = [1 : length(ucolors)]';
        layout.labelColor(:,sgn) = ucolors_min(ib);

        % Don't forget to update the colors of the terms themselves too
        [ucolors_term, ~, ib_term] = unique(layout.groupColor(:,sgn));
        ucolors_term_new = ones(size(ucolors_term));
        [~,ind1,ind2] = intersect(ucolors,ucolors_term);
        ucolors_term_new(ind2) = ucolors_min(ind1);
        layout.groupColor(:,sgn) = ucolors_term_new(ib_term);
        
        % And also the columns in the summary
        ucolors(ucolors==0) = [];
        layout.cumOpacityByColor{sgn} = layout.cumOpacityByColor{sgn}(:,ucolors);
        layout.cumOpacity01ByColor{sgn} = layout.cumOpacity01ByColor{sgn}(:,ucolors);

    else
        
        layout.groupColor = layout.groupColor+1;
        layout.labelColor = layout.labelColor+1;
        
    end
    
end

