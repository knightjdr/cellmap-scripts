function nodePositions = cyRestKK(adjacency, varargin)

% Show GUI-style progress bar
PROGRESSBAR = 0;
if ~isempty(find(strcmpi('ProgressBar', varargin)))
    PROGRESSBAR = varargin{find(strcmpi('ProgressBar', varargin))+1};   
end

% Import node positions
nodePositions = [];
if ~isempty(find(strcmpi('NodePositions', varargin)))
    nodePositions = varargin{find(strcmpi('NodePositions', varargin))+1};   
end

%%

url = 'http://localhost:1234/v1/';

% First, check if the cyREST server is running
try
    data = webread(url);
catch
    
    errortxt = 'To apply the Cytoscape implementation of the Kamada-Kawai algorith, make sure that Cytoscape is running in the background.\n\n';
    errortxt = [errortxt 'Cytoscape must be 3.2 or more recent. In Cytoscape 3.2, make sure the cyREST app is installed and running.\n\n'];
    errortxt = [errortxt 'For more information, visit https://github.com/idekerlab/cyREST/wiki/Installation-Guide.'];
    
    if PROGRESSBAR
        h = errordlg(sprintf(errortxt),'Cytoscape must be running');
        return;
    else
        error(sprintf(errortxt));
        return;
    end
end

stepname = 'Preparing the data to be passed to Cytoscape...';
if PROGRESSBAR
    w = waitbar(0.2, stepname);
else
    fprintf(['\n' stepname '\n']);
end

network_json = network2cytoscape(adjacency, nodePositions);

options = weboptions('MediaType','application/json');

stepname = 'Creating the network in Cytoscape...';
if PROGRESSBAR
    waitbar(0.4, w, stepname);
else
    fprintf(['\n' stepname '\n']);
end

data = webwrite([url 'networks/'], network_json, options);
networkSUID = data.networkSUID;

stepname = 'Creating the view in Cytoscape...';
if PROGRESSBAR
    waitbar(0.6, w, stepname);
else
    fprintf(['\n' stepname '\n']);
end

try
    
    % For small networks, this is sufficient to create a view of the network.
    viewSUID = webread([url 'networks/' num2str(networkSUID) '/views/']);
    
catch
    
    % For large networks, Cytoscape (for some reason) discarts all node
    % position information while creating the view. So, a little more work
    % is required to assign the positions back to the nodes.
    
    % Create the view
    data = urlread2([url 'networks/' num2str(networkSUID) '/views/'],'POST');
    data = loadjson(data);
    viewSUID = data.networkViewSUID;
    
    % Get the SUIDs and labels of all nodes in the network
    data = urlread2([url 'networks/' num2str(networkSUID) '/tables/defaultnode/columns/SUID/']);
    data = loadjson(data);
    nodes_SUIDs = data.values';

    data2 = urlread2([url 'networks/' num2str(networkSUID) '/tables/defaultnode/columns/SHARED_NAME/']);
    data2 = loadjson(data2);
    nodes_names = data2.values';
    nodes_ids = cell2mat(cellfun(@(x) textscan(x,'n%d'), nodes_names));

    [~,ix] = sort(nodes_ids);
    nodes_ids = nodes_ids(ix);
    nodes_names = nodes_names(ix);
    nodes_SUIDs = nodes_SUIDs(ix);

    % Update their positions one-by-one
    % Note: there might be away to do this for all nodes at once, but the
    % json syntax is a bit complex at the moment, so we're going to do it
    % one node at a time
    header = http_createHeader('Content-Type','application/json');
    for i = 1 : length(nodes_SUIDs)
        data = urlread2([url 'networks/' num2str(networkSUID) '/views/' num2str(viewSUID) '/nodes/' num2str(nodes_SUIDs(i)) '/'],'GET');
        data = regexprep(data, ...
            '([\s\n]*)"visualProperty"(\s*):(\s*)"NODE_X_LOCATION",([\s\n]*)"value"(\s*):(\s*)([\w.]+)([\s\n]*)', ...
            ['"visualProperty" : "NODE_X_LOCATION", "value" : ' sprintf('%.3f', nodePositions(i,1))]);
            data = regexprep(data, ...
            '([\s\n]*)"visualProperty"(\s*):(\s*)"NODE_Y_LOCATION",([\s\n]*)"value"(\s*):(\s*)([\w.]+)([\s\n]*)', ...
            ['"visualProperty" : "NODE_Y_LOCATION", "value" : ' sprintf('%.3f', nodePositions(i,2))]);
        [output, extras] = urlread2([url 'networks/' num2str(networkSUID) '/views/' num2str(viewSUID) '/nodes/' num2str(nodes_SUIDs(i)) '/'], ...
            'PUT', data, header);
    end

end

% Disable graph randomization because it is taken care by Matlab and controlled
% by its random number generator
% Must use urlread2 because webwrite does not support PUT
paramjson = urlread2([url 'apply/layouts/kamada-kawai/parameters/']);
param = loadjson(paramjson);
for i = 1 : length(param)
    if strcmp(param{i}.name, 'randomize')
        param{i}.value = false;
    end
end
paramjson = savejson('',param, 'ParseLogical', 1);

header = http_createHeader('Content-Type','application/json');
[output, extras] = urlread2([url 'apply/layouts/kamada-kawai/parameters/'],'PUT', paramjson, header); 

stepname = 'Applying the network layout (may take a few minutes)...';
if PROGRESSBAR
    waitbar(0.8, w, stepname);
else
    fprintf(['\n' stepname '\n']);
end

options = weboptions('Timeout', 12000);
data = webread([url 'apply/layouts/kamada-kawai/' num2str(networkSUID) '/'],options);

stepname = 'Saving the node positions (may take a few minutes)...';
if PROGRESSBAR
    waitbar(1, w, stepname);
else
    fprintf(['\n' stepname '\n']);
end

data = webread([url 'networks/' num2str(networkSUID) '/views/first/']);

% There seems to be difference in how the data comes in on Windows & Mac
if isfield(data, 'elements')
    data = data.elements;   % windows
else
    data = data.data.x__Annotations{8}; % mac
end

nodeNames = zeros(length(data.nodes),1);
nodePositions = zeros(length(data.nodes),2);
for i = 1 : length(data.nodes)
    
    nodeNames(i) = sscanf(data.nodes(i).data.name, 'n%d');
    nodePositions(i,:) = [data.nodes(i).position.x data.nodes(i).position.y];
    
end

[~,ix] = sort(nodeNames);
nodePositions = nodePositions(ix,:);

if PROGRESSBAR
    delete(w);
end


