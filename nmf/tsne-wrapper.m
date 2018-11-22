% Uses basis matrix for tSNE

% Load data
dataFolder = '/Volumes/Samsung external/Dropbox/Projects/Cell-map/analysis/2018_5_10/nmf/nmf_rank20/'
inputFile = strcat(dataFolder, 'basis.csv')
outputFile = strcat(dataFolder, 'tsne_nmf.txt')
basisTable = readtable(inputFile)

% subset table to array
noRanks = width(basisTable) - 1
rankNames = table2array(basisTable(1, 2:width(basisTable)))
rowNames = table2array(basisTable(2:height(basisTable), 1))
basisMatrix = table2array(basisTable(2:height(basisTable), 2:width(basisTable)))

% Set parameters
no_dims = 2;
initial_dims = noRanks;
perplexity = 5;

% Run t?SNE
mappedX = tsne(basisMatrix, [], no_dims, initial_dims, perplexity);

% format matrix for output
outputTable = table(rowNames, mappedX(:,1), mappedX(:,2))
outputTable.Properties.VariableNames = {'gene' 'x' 'y'}

% output table
writetable(outputTable,'tsne_nmf.txt','Delimiter','\t') 