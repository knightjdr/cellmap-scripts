function layout = initialize(path_to_results_folder, networkFile, annotationsMatrix, outputFolder, radius)

%% Check inputs

% This is mostly for Windows where filesep ('\') needs to be escaped in
% regular expressions.
filesep_escaped = regexptranslate('escape', filesep);

% If path to results folder is not provided, use the location of the /safe/ folder.
if isempty(path_to_results_folder)
    p = mfilename('fullpath');
    ind = regexp(p, [filesep_escaped 'io' filesep_escaped 'initialize']);
    path_to_results_folder = p(1:ind);
end

%% Output folder

if ~strcmp(path_to_results_folder(end), filesep)
    path_to_results_folder = [path_to_results_folder, filesep];
end
    
[pathstr, ~, ~] = fileparts(path_to_results_folder);
layout.outputdir = fullfile(pathstr, 'Results', outputFolder);
mkdir(layout.outputdir);

%% Load settings
[keys,~,~] = inifile([pathstr filesep 'safe.ini'],'readall');

for i = 1 : size(keys,1)
    if ~isempty(keys{i,4}) && ~isempty(str2num(keys{i,4}))
        keys{i,4} = str2num(keys{i,4});
    end
    layout.(keys{i,3}) = keys{i,4};
end

if isfield(layout, 'randomSeed') && ~isempty(layout.randomSeed)
    rng(layout.randomSeed);
end

%% Set network file and annotations matrix

layout.annotationfile = annotationsMatrix;
layout.neighborhoodRadius = radius;
layout.networkfile = networkFile;

%% Backward compatibility

if ~isfield(layout, 'layoutAlgorithm')
    layout.layoutAlgorithm = '';
end

if ~isfield(layout, 'neighborhoodRadiusType')
    layout.neighborhoodRadiusType = 'percentile';
end

if ~isfield(layout, 'unimodality')
    if isempty(layout.unimodalityType)
        layout.unimodality = 0;
    else
        layout.unimodality = 1;
    end
end

if ~isfield(layout, 'groupDistance')
    if isempty(layout.groupDistanceType)
        layout.groupDistance = 0;
    else
        layout.groupDistance = 1;
    end
end

if ~isfield(layout, 'plotNetwork')
    layout.plotNetwork = 1;
end

%% Cross-checks

if isempty(layout.annotationsign) || isempty(ismember(layout.annotationsign, {'both','highest','lowest'}))
    layout.annotationsign = 'highest';
    fprintf('\nSetting "annotationsign" unknown. Reverting to default "highest".\n');
end

if ~ismember(layout.neighborhoodRadiusType, {'percentile','absolute','diameter'})
    warning('Unknown value for parameter "neighborhoodRadiusType". Setting the parameter to default value: "percentile". See safe.ini for available options.');
    layout.neighborhoodRadiusType = 'percentile';
end
    
    
