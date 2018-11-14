function print_output_files(layout, varargin)

%% Check inputs

if ~isfield(layout, 'outputdir')
    p = mfilename('fullpath');
    ind = regexp(p, '/io/load_network');
    path_to_results_folder = p(1:ind);
    
    [pathstr, ~, ~] = fileparts(path_to_results_folder);
    layout.outputdir = [pathstr '/safe-' datestr(now,'yyyy-mm-dd-HH-MM-SS') '/'];  
end

OUTPUTDIR = layout.outputdir;

if ~exist(OUTPUTDIR, 'dir')
    mkdir(OUTPUTDIR);
end

PROGRESSBAR = 0;
if ~isempty(find(strcmpi('ProgressBar', varargin)))
    PROGRESSBAR = varargin{find(strcmpi('ProgressBar', varargin))+1};
end

%%

SGN = size(layout.pval,3);
nfiles = 3;
annotationsigns = {'highest'; 'lowest'};

if PROGRESSBAR
    w = waitbar(0, 'Exporting the results files...');
else
    fprintf('\nExporting the results files...\n');
end

for sgn = 1 : SGN
    
    %% Attribute file with node domain IDs and node neighborhood scores. Neighborhood score = opacity = -log10 maximum enrichment.

    i = 1;
    if PROGRESSBAR   
        w = waitbar( ((sgn-1)*nfiles+i) / (SGN*nfiles), w);
    end

    fid = fopen([OUTPUTDIR 'node_properties_annotation-' annotationsigns{sgn} '.txt'],'w');
    
    fprintf(fid, '## \n');
    fprintf(fid, '## This file lists the properties of all nodes in the network.\n');
    fprintf(fid, '## \n\n');
    
    columns = {'Node label', 'Node label ORF', 'Domain (predominant)', 'Neighborhood score [max=1, min=0] (predominant)', 'Total number of enriched domains','Number of enriched attributes per domain'};
    columns_field = {'label', 'label_orf', 'labelColor', 'labelOpacity'};
    columns_format = {'%s', '%s', '%d', '%.3f'};
    
    for k = 1 : length(columns)
        if k > 1, fprintf(fid, '\t'); end
        fprintf(fid, '%s', columns{k});
    end
    fprintf(fid, '\n');
    
    for i = 1 : length(layout.label)
        for k = 1 : length(columns_field)
            if k > 1, fprintf(fid, '\t'); end
            if strcmp('%s', columns_format{k})
                item = layout.(columns_field{k}){i};
            else
                item = layout.(columns_field{k})(i,sgn);
            end
            fprintf(fid, columns_format{k}, item);
        end
        
        fprintf(fid, '\t%d', length(find(layout.cumOpacity01ByColor{sgn}(i,:)>0)));
        
        fprintf(fid, '\t');
        for k = 1 : length(layout.regionId{sgn})
            if k > 1, fprintf(fid, ','); end
            fprintf(fid, '%d', layout.cumOpacity01ByColor{sgn}(i,k));
        end
        fprintf(fid, '\n');
    end
    fclose(fid);

    
    %% List of functional attributes organized by domains.
    
    i = 2;
    if PROGRESSBAR   
        w = waitbar( ((sgn-1)*nfiles+i) / (SGN*nfiles), w);
    end
    
    ucolors = layout.regionId{sgn};
    fid = fopen([OUTPUTDIR 'attribute_properties_annotation-' annotationsigns{sgn} '.txt'],'w');
    
    fprintf(fid, '## \n');
    fprintf(fid, '## This file lists the properties of all attributes used to annotate the network.\n');
    fprintf(fid, '## \n\n');
    
    fprintf(fid, 'Attribute Id\tAttribute name\tDomain Id\n');
    for i = 1 : length(ucolors)
        inds = find(layout.groupColor(:,sgn) == ucolors(i));
        for j = 1 : length(inds)
            fprintf(fid, '%d\t%s\t%d\n', layout.group_ids(inds(j)), layout.group_names{inds(j)}, ucolors(i));
        end
    end
    fclose(fid);

    %% List of RGB colors for all the domains
    
    i = 3;
    if PROGRESSBAR   
        w = waitbar( ((sgn-1)*nfiles+i) / (SGN*nfiles), w);
    end
    
    ucolors = layout.regionId{sgn};
    map_colors = layout.mapColors(2:end,:);
    fid = fopen([OUTPUTDIR 'domain_properties_annotation-' annotationsigns{sgn} '.txt'],'w');
    
    fprintf(fid, '## \n');
    fprintf(fid, '## This file lists the properties of all functional domains detected in the network.\n');
    fprintf(fid, '## \n\n');
    
    fprintf(fid, 'Domain number\tDomain name\tRGB\n');
    for i = 1 : length(ucolors)
        fprintf(fid, '%d\t%s\t%d\t%d\t%d\n', ucolors(i), layout.regionName{sgn}{i}, round(map_colors(i,1)*255), round(map_colors(i,2)*255), round(map_colors(i,3)*255));
    end
    fclose(fid);
    
end

if PROGRESSBAR
    delete(w);
end