function layout = safe(layout)

%% 2. Load network

layout = load_network(layout, 'ImportEdges', 1, 'IsDirected', 0);

%% 3. Apply network layout, if necessary

if ~isempty(layout.layoutAlgorithm)
    layout = generate_network_map(layout);
end

%% 4. Plot network

if layout.plotNetwork
    layout = plot_original_network(layout);
end

%% 5. Translate node labels, if necessary
% The label_orf field is used for enrichment analysis. Multiple nodes can have the same label_orf.

inds = find(cellfun(@isempty, layout.label_orf));

if ~isempty(inds)
    layout.label_orf(inds) = genename2orf(layout.label(inds));  % Yeast ORFs only
    tmp = split_by_delimiter('_', layout.label_orf(inds));
    layout.label_orf(inds) = tmp(:,1);
end

%% 6. Load functional annotations

layout = load_annotation_file(layout);

%% 7. Calculate distances between node pairs

layout = calculate_node_distances(layout);
    
%% 8. Calculate enrichments

layout = compute_enrichment(layout);

%% 9. Plot sample enrichment landscapes

if layout.plotSampleGroups
    
    % Plot 10 random groups
    plot_sample_groups(layout)
    
    % Plot 1 random group
    plot_sample_groups(layout,'Random',1);
    
    % To plot a specific set of groups, indicate their IDs
    % plot_sample_groups(layout,'Range',[1:2]);  
    
    % OPTIONAL. Save the current figure as PDF.
    export_figures(layout);
    
end

% OPTIONAL. Save all the enrichment scores to a text file.
write_neighborhood_scores(layout);

%% 10. Combine functional attributes into domains

% Identify region-specific vs multi-regional functional attributes
layout = compute_unimodality(layout);

% Collapse together attributes that have similar enrichment landscapes.
layout = cluster_groups(layout);

% Eliminate the least represented colors & re-sort the rest 
layout = minimize_colors(layout);

%% 11. Generate automatic domain labels

layout = generate_automatic_region_labels(layout);

%% 12. Plot annotated network

if strcmp('both', layout.annotationsign)
    layout = plot_composite_map(layout, 'AnnotationSign', 'highest');
    if layout.showLabels
        plot_labels(layout, 'AnnotationSign', 'highest');
    end
    layout = plot_composite_map(layout, 'AnnotationSign', 'lowest');
    if layout.showLabels
        plot_labels(layout, 'AnnotationSign', 'lowest');
    end
else
    layout = plot_composite_map(layout, 'AnnotationSign', layout.annotationsign);
    if layout.showLabels
        plot_labels(layout, 'AnnotationSign', layout.annotationsign);
    end
end

%% 13. Print results

print_output_files(layout);
export_figures(layout);

%% 14. Save SAFE session

layout = orderfields(layout);

save([layout.outputdir 'safe_session.mat'],'layout');




