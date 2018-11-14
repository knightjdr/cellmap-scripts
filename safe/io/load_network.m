function [layout, varargout] = load_network(layout, varargin)

%% Check inputs

% This is mostly for Windows where filesep ('\') needs to be escaped in regular expressions.
filesep_escaped = regexptranslate('escape', filesep);

IMPORTEDGES = 1;
if ~isempty(find(strcmpi('ImportEdges', varargin),1))
    IMPORTEDGES = varargin{find(strcmpi('ImportEdges', varargin))+1};
    if ismember(IMPORTEDGES, [0,1])
        IMPORTEDGES = (IMPORTEDGES == 1);
    else
        error('Check your inputs. The ImportEdges property expects a binary (0 or 1) or a logical (true or false) value.');
    end
end

IMPORTNODEATTRS = 1;
if ~isempty(find(strcmpi('ImportNodeAttributes', varargin),1))
    IMPORTNODEATTRS = varargin{find(strcmpi('ImportNodeAttributes', varargin))+1};
    if ismember(IMPORTNODEATTRS, [0,1])
        IMPORTNODEATTRS = (IMPORTNODEATTRS == 1);
    else
        error('Check your inputs. The ImportNodeAttributes property expects a binary (0 or 1) or a logical (true or false) value.');
    end
end

ISDIRECTED = 0;
if ~isempty(find(strcmpi('IsDirected', varargin),1))
    ISDIRECTED = varargin{find(strcmpi('IsDirected', varargin))+1};
end

if ~isfield(layout, 'outputdir')
    p = mfilename('fullpath');
    ind = regexp(p, [filesep_escaped 'io' filesep_escaped 'load_network']);
    path_to_results_folder = p(1:ind);

    [pathstr, ~, ~] = fileparts(path_to_results_folder);
    layout.outputdir = [pathstr filesep 'safe-' datestr(now,'yyyy-mm-dd-HH-MM-SS') filesep];
end

OUTPUTDIR = layout.outputdir;

PROGRESSBAR = 0;
if ~isempty(find(strcmpi('ProgressBar', varargin),1))
    PROGRESSBAR = varargin{find(strcmpi('ProgressBar', varargin))+1};
end

datadir = ['data' filesep];

%%

% If the input file is not specified, use the default Costanzo et al., 2010, network.
if isempty(layout.networkfile) || strcmp(layout.networkfile, 'Genetic Interaction Similarity Network (Costanzo~Boone, 2010)')

    fprintf('\nLoading the genetic interaction similarity network (Costanzo~Boone, 2010)...\n');
    temp = load([datadir 'layout_Costanzo2010_150831.mat']);
    fieldNames = fieldnames(temp.layout);
    for f = 1 : length(fieldNames)
        layout.(fieldNames{f}) = temp.layout.(fieldNames{f});
    end

    varargout{1} = true;    % Successflag

    return;

end

%% Check format of the file

[~, ~, ext] = fileparts(layout.networkfile);

if strcmpi(ext, '.mat')

    fprintf('\n\nLoading a mat-file %s...\n', layout.networkfile);

    temp = load(layout.networkfile);

    fieldNames = fieldnames(temp.layout);
    for f = 1 : length(fieldNames)
        layout.(fieldNames{f}) = temp.layout.(fieldNames{f});
    end

    return;

end

if ~strcmpi(ext, '.cys')
    % Get the first line to figure out how many columns are there in the file
    [fid, message ]= fopen(layout.networkfile, 'r');
    tLines = fgets(fid);
    numCols = numel(regexp(tLines,'\t'))+1;
    fclose(fid);

    if numCols == 3

        fprintf('\n\nLoading a 3-column network file %s...\n', layout.networkfile);

        fid = fopen(layout.networkfile, 'r');
        C = textscan(fid, '%s %s %f','delimiter','\t');
        fclose(fid);

        layout.label = unique([C{1}; C{2}]);
        layout.label_orf = layout.label;
        layout.edges_weights = nan(length(layout.label));
        layout.edges = zeros(size(layout.edges_weights));

        [~, Locb1] = ismember(C{1}, layout.label);
        [~, Locb2] = ismember(C{2}, layout.label);

        Locb = sub2ind(size(layout.edges_weights), Locb1, Locb2);
        layout.edges_weights(Locb) = C{3};
        layout.edges(layout.edges_weights > 0) = 1;

    elseif numCols == 5

        fprintf('\n\nLoading a 5-column network file %s...\n', layout.networkfile);

        fid = fopen(layout.networkfile, 'r');
        C = textscan(fid, '%s %s %s %s %f','delimiter','\t');
        fclose(fid);

        layout.label = unique([C{1}; C{3}]);
        layout.label_orf = cell(size(layout.label));
        layout.edges_weights = nan(length(layout.label));
        layout.edges = zeros(size(layout.edges_weights));

        [~, Locb1] = ismember(C{1}, layout.label);
        [~, Locb2] = ismember(C{3}, layout.label);

        layout.label_orf(Locb1) = C{2};
        layout.label_orf(Locb2) = C{4};

        Locb = sub2ind(size(layout.edges_weights), Locb1, Locb2);
        layout.edges_weights(Locb) = C{5};
        layout.edges(layout.edges_weights > 0) = 1;

    else

        fprintf('\n\nLoading a network file as an adjacency matrix %s...\n', layout.networkfile);

        D = read_matrix_file(layout.networkfile, 1,1);

        layout.label = unique([D.labels_row; D.labels_col]);
        layout.label_orf = layout.label;
        layout.edges = nan(length(layout.label));

        [~,ind1,ind2] = intersect(layout.label, D.labels_row);
        [~,ind3,ind4] = intersect(layout.label, D.labels_col);

        layout.edges(ind1,ind3) = D.data(ind2,ind4);

    end

    if ~ISDIRECTED
        t = cat(3, layout.edges, layout.edges');
        layout.edges = nanmean(t,3);
    end

    layout.edges(isnan(layout.edges)) = 0;

    varargout{1} = true;    % Successflag

    return;

end

%% CYTOSCAPE SESSION

%% Unpack the CYS file

if ~exist(OUTPUTDIR, 'dir')
    mkdir(OUTPUTDIR);
end

filenames = unzip(layout.networkfile, OUTPUTDIR);

%% Load the network data (nodes and edges)

if IMPORTEDGES

    network_found = false;

    % If the networkname was provided, use that one.
    if ~isempty(layout.networkname)
         % Correct the file names, in case they contain unusual characters
        network_name_file = layout.networkname;
        if ~isempty(strfind(layout.networkname, '('))
            network_name_file = regexprep(network_name_file,'(','%28');
        end
        if ~isempty(strfind(layout.networkname, ')'))
            network_name_file = regexprep(network_name_file,')','%29');
        end

        rx = [filesep_escaped 'networks' filesep_escaped '[a-z0-9\-]*' network_name_file '.xgmml'];
        f = find(~cellfun(@isempty, regexp(filenames, rx)));

        if ~isempty(f)
            network_found = true;
        else
            fprintf('\nThe provided network name didn''t match any existing networks. Will check other options.\n');
        end
    end

    % If the provided network name didn't match anything or was not
    % provided at all, check what's available in the session file.
    if isempty(layout.networkname) || ~network_found

        rx = [filesep_escaped 'networks' filesep_escaped];
        inds = find(~cellfun(@isempty, regexp(filenames, rx)));

        switch length(inds)
            case 0
                error('No networks in this CYS file.');
            case 1
                f = inds(1);
                fprintf('\nThere is only one network in this CYS file, so taking that one:\n\t%s\n', filenames{f});
            otherwise
                network_list = [];
                for i = 1 : length(inds)
                    [~, network_list{i}, ~] = fileparts(filenames{inds(i)});
                end

                if PROGRESSBAR
                    i = cytoscape_network_question('Title', 'Cytoscape Session File', ...
                        'String', 'Several networks exist in this CYS file. Please, choose one:', ...
                        'List', network_list);
                    if i < 0
                        varargout{1} = false;   % Success flag
                        return;
                    end
                else
                    fprintf('\nSeveral networks exist in this CYS file:\n');
                    for i = 1 : length(inds)
                        fprintf('\t%d. %s\n', i, network_list{i});
                    end
                    str = sprintf('\nWhich network would you like to use? 1-%d [1]: ', length(inds));
                    i = input(str);
                    if isempty(i) || ~isnumeric(i)
                        i = 1;
                    end
                end

                f = inds(i);

        end

    end

    [~, cytoscapeNetworkFile, ~] = fileparts(filenames{f});
    c = textscan(cytoscapeNetworkFile, '%d-%s');
    cytoscapeNetworkId = c{1};
    cytoscapeNetworkName = c{2};

    if PROGRESSBAR
        w = waitbar(0.25, 'Loading the nodes and edges (may take a few minutes)...');
        set(get(findobj(w,'Type','axes'),'Title'),'FontSize',12);
    else
        fprintf('\nLoading network edges... (may take a few minutes)\n');
    end

    fprintf('\t\tParsing the network file...\n');
    ntw_conn = parseXML(filenames{f});

    chld = {ntw_conn.Children.Name};
    inds = find(strcmp('att', chld));

    for j = 1 : length(inds)
        chld2 = {ntw_conn.Children(inds(j)).Children.Name};
        inds2 = find(strcmp('graph', chld2));

        for k = 1 : length(inds2)
            intw = (j-1) * length(inds2) + k;

            attributes = {ntw_conn.Children(inds(j)).Children(inds2(k)).Attributes.Name};
            attributes_values = {ntw_conn.Children(inds(j)).Children(inds2(k)).Attributes.Value};

            % Record the network ID
            inds3 = find(strcmp('id', attributes));
            if ~isempty(inds3)
                cys(intw).graphId = str2num(attributes_values{inds3});
            end

            % Record the network name
            inds3 = find(strcmp('label', attributes));
            if ~isempty(inds3)
                cys(intw).graphName = attributes_values{inds3};
            end

            chld3 = {ntw_conn.Children(inds(j)).Children(inds2(k)).Children.Name};

            % Get the list of nodes
            inds4 = find(strcmp('node', chld3));
            cys(intw).nodeId = nan(length(inds4),1);
            cys(intw).nodeLabel = cell(length(inds4),1);

            ninds4 = length(inds4);

            if ninds4 > 0

                % Check the first node
                attributes = {ntw_conn.Children(inds(j)).Children(inds2(k)).Children(inds4(1)).Attributes.Name};
                attributes_values = {ntw_conn.Children(inds(j)).Children(inds2(k)).Children(inds4(1)).Attributes.Value};

                i1 = strcmp('id', attributes);
                i2 = strcmp('label', attributes);
                i3 = strcmp('xlink:href', attributes);

                if ~isempty(find(i1,1))
                    intw_ref = intw;    % this is the reference network that contains full info about the nodes and edges
                elseif ~isempty(find(i3,1))
                    if ~strcmp(cys(intw).graphName, layout.viewname)    % If this is not the network chosen by the user, don't load anything and skip.
                        continue;
                    end
                end

                fprintf('\n\t\tReading node info...\n');
                fprintf(['\t\t|', blanks(100), '|\n']);
                fprintf('\t\t|');
                y = 0;

                for h = 1 : ninds4

                    % -- Progress bar ---
                    if PROGRESSBAR
                        % nothing for now
                    else
                        x = round(h * 100 / ninds4);
                        if x > y
                            fprintf(repmat('*',1,x-y)); y = x;
                        end
                    end
                    % -- Progress bar ---

                    attributes = {ntw_conn.Children(inds(j)).Children(inds2(k)).Children(inds4(h)).Attributes.Name};
                    attributes_values = {ntw_conn.Children(inds(j)).Children(inds2(k)).Children(inds4(h)).Attributes.Value};

                    i1 = strcmp('id', attributes);
                    i2 = strcmp('label', attributes);
                    i3 = strcmp('xlink:href', attributes);  % extracted networks link to the the original network they were extracted from

                    if ~isempty(find(i1,1))
                        cys(intw).nodeId(h) = str2num(attributes_values{i1});
                        cys(intw).nodeLabel(h) = attributes_values(i2);
                    elseif ~isempty(find(i3,1))
                        nodeId = textscan(attributes_values{i3}, '#%d');
                        nodeId = nodeId{1};
                        cys(intw).nodeId(h) = nodeId;
                        cys(intw).nodeLabel(h) = cys(intw_ref).nodeLabel(cys(intw_ref).nodeId == nodeId);
                    end
                end

                if PROGRESSBAR
                    % nothing for now
                else
                    fprintf('|\n');
                end

            end

            % Get the list of edges
            inds4 = find(strcmp('edge', chld3));
            cys(intw).edgeId = zeros(length(cys(intw).nodeId));

            ninds4 = length(inds4);

            if ninds4 > 0

                fprintf('\n\t\tReading edge info...\n');
                fprintf(['\t\t|', blanks(100), '|\n']);
                fprintf('\t\t|');
                y = 0;

                for h = 1 : ninds4

                    % -- Progress bar ---
                    if PROGRESSBAR
                        % nothing for now
                    else
                        x = round(h * 100 / ninds4);
                        if x > y
                            fprintf(repmat('*',1,x-y)); y = x;
                        end
                    end
                    % -- Progress bar ---

                    attribute_names = {ntw_conn.Children(inds(j)).Children(inds2(k)).Children(inds4(h)).Attributes.Name};
                    attribute_values = {ntw_conn.Children(inds(j)).Children(inds2(k)).Children(inds4(h)).Attributes.Value};

                    i1 = strcmp('id', attribute_names);
                    i2 = strcmp('source', attribute_names);
                    i3 = strcmp('target', attribute_names);
                    i4 = strcmp('xlink:href', attribute_names);  % extracted networks link to the the original network they were extracted from

                    if ~isempty(find(i1,1))
                        edgeId = str2num(attribute_values{i1});
                        sourceId = str2num(attribute_values{i2});
                        targetId = str2num(attribute_values{i3});
                    elseif ~isempty(find(i4,1))
                        edgeId = textscan(attribute_values{i4}, '#%d');
                        edgeId = edgeId{1};
                        [r,c] = find(cys(intw_ref).edgeId == edgeId);
                        sourceId = cys(intw_ref).nodeId(r);
                        targetId = cys(intw_ref).nodeId(c);
                    end

                    cys(intw).edgeId(cys(intw).nodeId == sourceId,cys(intw).nodeId == targetId) = edgeId;

                end

                if PROGRESSBAR
                    % nothing for now
                else
                    fprintf('|\n');
                end

            end

        end
    end

    if PROGRESSBAR
        delete(w);
    end
end

%% Load the view (node coordinates)

view_found = false;

% If the view name was provided, use that one.
if ~isempty(layout.viewname)
    % Correct the file names, in case they contain unusual characters
    view_name_file = layout.viewname;
    if ~isempty(strfind(layout.viewname, '('))
        view_name_file = regexprep(view_name_file,'(','%28');
    end
    if ~isempty(strfind(layout.viewname, ')'))
        view_name_file = regexprep(view_name_file,')','%29');
    end

    rx = [filesep_escaped 'views' filesep_escaped view_name_file '.xgmml'];
    f = find(~cellfun(@isempty, regexp(filenames, rx)));     % Cytoscape 3.1.0
    if isempty(f)
        rx = [filesep_escaped 'views' filesep_escaped '[a-z0-9\-]*' view_name_file '(\+View)*.xgmml'];
        f = find(~cellfun(@isempty, regexp(filenames, rx)));    % Cytoscape 3.2.1
    end

    if ~isempty(f)
        view_found = true;
    else
        fprintf('\nThe provided view name didn''t match any existing views. Will check other options.\n');
    end
end

% If the provided view name didn't match anything or was not
% provided at all, check what's available in the session file.
if isempty(layout.viewname) || ~view_found

    % Generate the regular expression to find the views that correspond to the selected network
    cysGraphIds = cellfun(@num2str, {cys(:).graphId},'UniformOutput',0);
    rx = [filesep_escaped 'views' filesep_escaped '(' strjoin(strcat(cysGraphIds,{'-'}),'|') ')'];

    inds = find(~cellfun(@isempty, regexp(filenames, rx)));

    switch length(inds)
        case 0
            error('No views in this CYS file.');
        case 1
            f = inds(1);
            fprintf('\nThere is only one view in this CYS file, so taking that one:\n\t%s\n', filenames{f});
        otherwise

            view_list = [];
            for i = 1 : length(inds)
                [~, view_list{i}, ~] = fileparts(filenames{inds(i)});
            end

            if PROGRESSBAR
                i = cytoscape_network_question('Title', 'Cytoscape Session File', ...
                        'String', 'Several views exist for the network you selected. Please, choose one:', ...
                        'List', view_list);
                if i < 0
                    varargout{1} = false;   % Success flag
                    return;
                end
            else
                fprintf('\nSeveral views exist in this CYS file:\n');
                for i = 1 : length(inds)
                    fprintf('\t%d. %s\n', i, view_list{i});
                end
                str = sprintf('\nWhich view would you like to use? 1-%d [1]: ', length(inds));
                i = input(str);
                if isempty(i) || ~isnumeric(i)
                    i = 1;
                end
            end
            f = inds(i);
    end
end

if PROGRESSBAR
    w = waitbar(0.5, 'Loading the network view...');
    set(get(findobj(w,'Type','axes'),'Title'),'FontSize',12);
else
    fprintf('\nLoading the network view... (may take a few minutes)\n');
end

[~, s, ~] = fileparts(filenames{f});
c = textscan(s, '%d-%d-%s');
cysGraphId = c{1};

intw = find([cys(:).graphId] == cysGraphId);

cys(intw).viewFile = filenames{f};
cys(intw).viewId = c{2};
cys(intw).viewName = c{3}{1};

cys(intw).x = nan(length(cys(intw).nodeId),1);
cys(intw).y = nan(length(cys(intw).nodeId),1);

ntw = parseXML(filenames{f});

chld = {ntw.Children.Name};
inds = find(strcmp('node', chld));

for i = 1 : length(inds)
    attribute_names = {ntw.Children(inds(i)).Attributes.Name};
    attribute_values = {ntw.Children(inds(i)).Attributes.Value};

    nodeId = str2num(attribute_values{strcmp('cy:nodeId', attribute_names)});

    children_names = {ntw.Children(inds(i)).Children.Name};
    igraphics = find(strcmp('graphics', children_names));

    graphics_attribute_names = {ntw.Children(inds(i)).Children(igraphics).Attributes.Name};
    graphics_attribute_values = {ntw.Children(inds(i)).Children(igraphics).Attributes.Value};

    x = str2num(graphics_attribute_values{strcmp('x', graphics_attribute_names)});
    y = str2num(graphics_attribute_values{strcmp('y', graphics_attribute_names)});

    cys(intw).x(cys(intw).nodeId == nodeId) = x;
    cys(intw).y(cys(intw).nodeId == nodeId) = y;
end

if PROGRESSBAR
    delete(w);
end

%% Get label ORFs

if IMPORTNODEATTRS

    % Get the exact network name
    exact_network_name = [num2str(cytoscapeNetworkId) '-' cytoscapeNetworkName{1}];
    rx = [filesep 'tables' filesep exact_network_name filesep 'SHARED_ATTRS-org.cytoscape.model.CyNode'];
    f = find(~cellfun(@isempty, regexp(filenames, regexptranslate('escape', rx))));

    % Get the second line to figure out how many columns are there in the
    % file
    fid = fopen(filenames{f},'r');
    tLines = fgets(fid);
    tLines = fgets(fid);
    numCols = numel(strfind(tLines,','))+1;
    fclose(fid);

    fid = fopen(filenames{f},'r');
    C = textscan(fid, repmat('%q ', 1, numCols),'delimiter',',');
    fclose(fid);

    for i = 1 : numCols
        tmp_struct.(matlab.lang.makeValidName(C{i}{2})) = C{i};
    end

    tmp_struct.SUID = cellfun(@str2num, tmp_struct.SUID,'UniformOutput',0);
    inds = find(cellfun(@isempty, tmp_struct.SUID));
    flds = fieldnames(tmp_struct);
    for i = 1 : length(flds)
        tmp_struct.(flds{i})(inds) = [];
    end
    tmp_struct.SUID = cell2mat(tmp_struct.SUID);

    cys(intw).nodeLabelSystematic = cell(size(cys(intw).nodeId));

    if (~isempty(find(ismember(fieldnames(tmp_struct), 'sharedName'))) | ~isempty(find(ismember(fieldnames(tmp_struct), 'canonicalName')))) && ~isempty(find(ismember(fieldnames(tmp_struct),'ORF')))
        [~,ind1,ind2] = intersect(cys(intw).nodeId, tmp_struct.SUID);
        cys(intw).nodeLabelSystematic(ind1) = tmp_struct.ORF(ind2);
    end

    inds = find(cellfun(@isempty, cys(intw).nodeLabelSystematic));

    if ~isempty(inds)
        warning('%d out of %d nodes don''t have an ORF attribute.', length(inds), length(cys(intw).nodeId));
    end

end

%% Save the relevant information to the layout structure

layout.networkname = cys(intw).graphName;
layout.viewname = cys(intw).viewName;

layout.label = cys(intw).nodeLabel;
layout.label_orf = cys(intw).nodeLabelSystematic;

layout.x = cys(intw).x;
layout.y = cys(intw).y;

layout.edges = zeros(size(cys(intw).edgeId));
layout.edges(cys(intw).edgeId>0) = 1;

if ~ISDIRECTED
    layout.edges = layout.edges + layout.edges';
    layout.edges(layout.edges > 1) = 1;
end

layout.edges = uint8(layout.edges);

%%

if ~PROGRESSBAR
    % Final report
    fprintf('\n\nNetwork loaded.\n');
    fprintf('Number of nodes: %d\n', length(layout.label));
    fprintf('Number of edges: %d\n', length(find(layout.edges>0)));
    if issymmetric(double(layout.edges))
        fprintf('The edges are not directed.\n');
    else
        fprintf('The edges are directed.\n');
    end
end

varargout{1} = true;
