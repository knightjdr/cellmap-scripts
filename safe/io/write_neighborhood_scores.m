function write_neighborhood_scores(layout, varargin)

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
annotationsigns = {'highest'; 'lowest'};

for sgn = 1 : SGN

    % Attribute file with node info. Neighborhood score = opacity = -log10 enrichment.
    fid = fopen([OUTPUTDIR 'neighborhood_scores_annotation-' annotationsigns{sgn} '.txt'],'w');
    
    fprintf(fid, '## \n');
    fprintf(fid, '## This file lists the -log10 scores that represent the enrichment of each node''s neighborhood for each of the tested attributes.\n');
    fprintf(fid, '## High values (close to 1) indicate high enrichment. Low values (close to 0) indicate low enrichment.\n');
    fprintf(fid, '## Any value higher than %.3f is considered significant (-log10 of p-value = 0.05, corrected for multiple testing and scaled to [0,1] range).\n', layout.thresholdOpacity);
    fprintf(fid, '## \n\n');
    
    write_matrix_file(fid, layout.label, layout.group_names, layout.opacity(:,:,sgn), ...
        'ProgressBar', PROGRESSBAR, 'ProgressTitle', ['Writing neighborhood scores for the ' annotationsigns{sgn} ' values...']);
    
    fclose(fid);
    
end