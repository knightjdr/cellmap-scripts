%% Spatial Analysis of Functional Enrichment (SAFE)
%
% Spatial Analysis of Functional Enrichment (SAFE) is an automated network annotation algorithm. It performs local
% enrichment analysis on a biological network to determine which regions of
% the network are over-represented for functional attributes (binary annotations or quantitative phenotypes).
%
% Author: Anastasia Baryshnikova (2016), abarysh@princeton.edu,
% http://www.baryshnikova-lab.org
%
% Modified for interating over a folder of networks and annotation matrices
% by James Knight jknight@lunenfeld.ca
%% 1a. Setup iterator

% parent folder containing files to analyze and where results will be
% stored. Also contains .ini file.
analysisFolder = '/Volumes/Samsung external/Dropbox/Projects/Cell-map/analysis/2018_11_13/safe';
% subfolder containing network files
networkFolder = 'networks/';
% subfolder containg annotation matrices
annotationMatrices = 'matrices/';
% radii to test
radii = [1 1.5 2 2.5 3 3.5 4 4.5 5];

% get network files
networkFullPath = fullfile(analysisFolder, networkFolder);
networkFullPathCys = fullfile(networkFullPath, '*.cys');
networkFiles = dir(networkFullPathCys);

% get matrix files
matrixFullPath = fullfile(analysisFolder, annotationMatrices);
matrixFullPathCys = fullfile(matrixFullPath, '*.txt');
matrixFiles = dir(matrixFullPathCys);

%% 1b. Iterate

for i = 1:length(networkFiles)
    netFile = networkFiles(i).name
    netfullFilePath = fullfile(networkFullPath, netFile);
    corrCutoff = strtok(netFile, '_');
    for j = 1:length(matrixFiles)
        matFile = matrixFiles(j).name
        matfullFilePath = fullfile(matrixFullPath, matFile);
        namespace = strtok(matFile, '_');
        for k = 1:length(radii)
            outputFolder = strcat(corrCutoff, 'cc_', namespace, 'ns_', num2str(radii(k)), 'r', '/');
            fprintf('\n\nLoading a mat-file %s...\n', outputFolder);
            %% 1. Run SAFE
            layout = initialize(analysisFolder, netfullFilePath, matfullFilePath, outputFolder, radii(k));
            safe(layout)
        end
    end
end
