function layout = load_annotation_file(layout, varargin)

%% Check inputs

if mod(nargin,2) ~= 1 
    error('Check your inputs. This function requires a set of property name/property value pairs.');  
end

PATHTOANNOTATION = '';
if isfield(layout, 'annotationfile')
    PATHTOANNOTATION = layout.annotationfile;
end

PROGRESSBAR = 0;
if ~isempty(find(strcmpi('ProgressBar', varargin)))
    PROGRESSBAR = varargin{find(strcmpi('ProgressBar', varargin))+1};
end

datadir = 'data/';

%%

if isempty(PATHTOANNOTATION) || strcmp(PATHTOANNOTATION, 'GO Biological Process (S. cerevisiae)') 
    PATHTOANNOTATION = [datadir 'go_bp_140819.mat'];
elseif strcmp(PATHTOANNOTATION, 'Hoepfner~Movva, 2014')
    PATHTOANNOTATION = [datadir 'hoepfner_movva_2014_hop_known.txt'];
end

[~,~,ext] = fileparts(PATHTOANNOTATION);

if strcmp('.mat',ext)
    
    if PROGRESSBAR
        w = waitbar(PROGRESSBAR, 'Loading the Gene Ontology Biological Process standard...');
    else
        fprintf('\nLoading the Gene Ontology Biological Process standard...\n');
    end

    % Example: Gene Ontology (biological process)
    load(PATHTOANNOTATION);

    layout.group_ids = go.term_ids;    % unique group identifiers
    layout.group_names = go.term_names;    % group names (don't have to be unique)
    layout.label2group = zeros(length(layout.label), length(layout.group_ids)); % matrix containing label-to-group mappings

    for i = 1 : length(layout.group_ids)
        inds = ismember(layout.label_orf, go.orfs(go.term2orf(i,:)>0));
        layout.label2group(inds,i) = 1;
    end
    
    layout.group_global_size = sum(go.term2orf,2);
    layout.group_total_gene_number = length(go.orfs);
    
    if PROGRESSBAR
        delete(w);
    end
    
elseif strcmp('.txt', ext)
    
    if PROGRESSBAR
        w = waitbar(PROGRESSBAR, 'Loading the functional attributes...');
    else
        fprintf('\nLoading the functional attributes from %s...\n', PATHTOANNOTATION);
    end
        
    % Text file format: see |sample_annotation_file.txt| in the |/data/|
    % folder.
    D = read_matrix_file(PATHTOANNOTATION,1,1);
    
    % Validity check: if the rows (gene labels) in the annotation file are not unique, take the average of the corresponding values
    if length(unique(D.labels_row)) < length(D.labels_row)
        warning('\nThe row labels (node identifiers) are not unique in the annotation file. Calculating the average of their annotation values...\n');
        [t,t2] = grpstats(D.data, D.labels_row, {'gname','mean'});
        D.labels_row = t;
        D.data = t2;
    end
        
    layout.group_ids = [1:length(D.labels_col)]';
    layout.group_names = D.labels_col;
    layout.label2group = zeros(length(layout.label), length(layout.group_ids));
    
    [~,ix] = ismember(layout.label_orf, D.labels_row);
    layout.label2group(ix>0,:) = D.data(ix(ix>0),:);
    
    layout.group_global_size = sum(D.data,1)';
    layout.group_total_gene_number = length(D.labels_row);
    
    if PROGRESSBAR
        delete(w);
    end
    
else
    error('Annotation file: format unknown. Check the README file for a list of accepted formats.');
end

layout.label2group(isnan(layout.label2group)) = 0;
layout.groupIsTop = ones(size(layout.group_ids));


% Final report
fprintf('\nFunctional standard loaded.\n');
fprintf('Number of functional groups: %d\n', length(layout.group_ids));
fprintf('Number of annotations: \n');
fprintf('\t\t%d NaNs\n', length(find(isnan(layout.label2group))));
fprintf('\t\t%d not NaNs\n', length(find(~isnan(layout.label2group))));
fprintf('\t\t%d zeros\n', length(find(layout.label2group == 0)));
fprintf('\t\t%d positives\n', length(find(layout.label2group > 0)));
fprintf('\t\t%d negatives\n', length(find(layout.label2group < 0)));
fprintf('\n');

