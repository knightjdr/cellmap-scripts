function layout = load_network(layout, varargin)

%% Check inputs

IMPORTEDGES = 1;
if ~isempty(find(strcmpi('ImportEdges', varargin)))
    IMPORTEDGES = varargin{find(strcmpi('ImportEdges', varargin))+1};
    if ismember(IMPORTEDGES, [0,1])
        IMPORTEDGES = (IMPORTEDGES == 1);
    else
        error('Check your inputs. The ImportEdges property expects a binary (0 or 1) or a logical (true or false) value.');
    end   
end

IMPORTNODEATTRS = 0;
if ~isempty(find(strcmpi('ImportNodeAttributes', varargin)))
    IMPORTNODEATTRS = varargin{find(strcmpi('ImportNodeAttributes', varargin))+1};
    if ismember(IMPORTNODEATTRS, [0,1])
        IMPORTNODEATTRS = (IMPORTNODEATTRS == 1);
    else
        error('Check your inputs. The ImportNodeAttributes property expects a binary (0 or 1) or a logical (true or false) value.');
    end   
end

ISDIRECTED = 0;
if ~isempty(find(strcmpi('IsDirected', varargin)))
    ISDIRECTED = varargin{find(strcmpi('IsDirected', varargin))+1};
end

OUTPUTDIR = pwd;
if ~isempty(find(strcmpi('OutputDir', varargin)))
    OUTPUTDIR = varargin{find(strcmpi('OutputDir', varargin))+1};
end
  
%%

% If a Cytoscape session is not specified, use the default Costanzo et al., 2010, network.
if isempty(layout.filename)
    
    load layout_Costanzo2010.mat
    return;
    
end

%% Get the X and Y coordinates of the nodes using the View file in the Cytoscape CYS session
filenames = unzip(layout.filename, OUTPUTDIR);

if isempty(layout.viewname)
    f = find(~cellfun(@isempty, regexp(filenames, '/views/')),1);
    fprintf('\nView name not provided. Taking the first one in the list: %s\n', filenames{f});
else
    % Correct the file names, in case they contain unusual characters
    view_name_file = layout.viewname;
    if ~isempty(strfind(layout.viewname, '('))
        view_name_file = regexprep(view_name_file,'(','%28');
    end
    if ~isempty(strfind(layout.viewname, ')'))
        view_name_file = regexprep(view_name_file,')','%29');
    end
    
    f = find(~cellfun(@isempty, regexp(filenames, ['/views/' view_name_file '.xgmml'])));     % Cytoscape 3.1.0
    if isempty(f)
        f = find(~cellfun(@isempty, regexp(filenames, ['/views/' view_name_file '\+View.xgmml'])));    % Cytoscape 3.2.1
    end
    
end

fprintf('\nLoading network layout... (may take a few minutes)\n');

layout.xySource = filenames{f};
ntw = parseXML(filenames{f});


chld = {ntw.Children.Name};
inds = find(strcmp('node', chld));

for i = 1 : length(inds)
    attribute_names = {ntw.Children(inds(i)).Attributes.Name};

    ilabel = find(strcmp('label', attribute_names));
    iid = find(strcmp('cy:nodeId', attribute_names));

    layout.id(i,1) = str2num(ntw.Children(inds(i)).Attributes(iid).Value);
    layout.label{i,1} = ntw.Children(inds(i)).Attributes(ilabel).Value;

    children_names = {ntw.Children(inds(i)).Children.Name};
    igraphics = find(strcmp('graphics', children_names));

    graphics_attribute_names = {ntw.Children(inds(i)).Children(igraphics).Attributes.Name};

    ix = find(strcmp('x', graphics_attribute_names));
    iy = find(strcmp('y', graphics_attribute_names));

    layout.x(i,1) = str2num(ntw.Children(inds(i)).Children(igraphics).Attributes(ix).Value);
    layout.y(i,1) = str2num(ntw.Children(inds(i)).Children(igraphics).Attributes(iy).Value);
end

%% Get edges

layout.edges = zeros(length(layout.label));

if IMPORTEDGES

    if isempty(layout.networkname)
        f = find(~cellfun(@isempty, regexp(filenames, '/networks/')),1);
        fprintf('\nNetwork name not provided. Taking the first one in the list: %s\n', filenames{f});
    else
        f = find(~cellfun(@isempty, regexp(filenames, ['/networks/' layout.networkname '.xgmml'])));
    end

    fprintf('\nLoading network edges... (may take a few minutes)\n');
    
    layout.edgesSource = filenames{f};
    ntw_conn = parseXML(filenames{f});

    chld = {ntw_conn.Children.Name};
    inds = find(strcmp('att', chld));

    for j = 1 : length(inds)
        chld2 = {ntw_conn.Children(inds(j)).Children.Name};
        inds2 = find(strcmp('graph', chld2));

        for k = 1 : length(inds2)
            attributes = {ntw_conn.Children(inds(j)).Children(inds2(k)).Attributes.Name};
            inds3 = find(strcmp('label', attributes));
            if isempty(layout.viewname) | strcmp(layout.viewname, ntw_conn.Children(inds(j)).Children(inds2(k)).Attributes(inds3).Value)
                chld3 = {ntw_conn.Children(inds(j)).Children(inds2(k)).Children.Name};

                inds4 = find(strcmp('edge', chld3));

                for h = 1 : length(inds4)
                    attributes = {ntw_conn.Children(inds(j)).Children(inds2(k)).Children(inds4(h)).Attributes.Name};
                    inds5 = find(strcmp('label', attributes));
                    tmp = regexp(ntw_conn.Children(inds(j)).Children(inds2(k)).Children(inds4(h)).Attributes(inds5).Value, ' ', 'split');

                    i1 = strcmp(tmp{1}, layout.label);
                    i2 = strcmp(tmp{3}, layout.label);

                    layout.edges(i1,i2) = 1;
                end
            end
        end
    end
    
    if ~ISDIRECTED
        layout.edges = layout.edges + layout.edges';
        layout.edges(layout.edges > 1) = 1;
    end

    layout.edges = uint8(layout.edges);
end
    
%% Get label ORFs

if IMPORTNODEATTRS

    % Get the exact network name
    [pathstr, exact_network_name, ext] = fileparts(filenames{f});

    f = find(~cellfun(@isempty, regexp(filenames, ['/tables/' exact_network_name '/SHARED_ATTRS-org.cytoscape.model.CyNode'])));

    % Get the second line to figure out how many columns are there in the
    % file
    fid = fopen(filenames{f},'r');
    tLines = fgets(fid);
    tLines = fgets(fid);
    numCols = numel(strfind(tLines,','))+1;
    fclose(fid);

    fid = fopen(filenames{f},'r');
    C = textscan(fid, repmat('%q ', 1,numCols),'delimiter',',');
    fclose(fid);

    for i = 1 : numCols
        tmp_struct.(matlab.lang.makeValidName(C{i}{2})) = C{i};
    end

    layout.label_orf = cell(size(layout.label));

    if (~isempty(find(ismember(fieldnames(tmp_struct), 'sharedName'))) | ~isempty(find(ismember(fieldnames(tmp_struct), 'canonicalName')))) && ~isempty(find(ismember(fieldnames(tmp_struct),'ORF')))
        [~,ind1,ind2] = intersect(layout.label, tmp_struct.sharedName);
        layout.label_orf(ind1) = tmp_struct.ORF(ind2);
    end

    inds = find(cellfun(@isempty, layout.label_orf));
    layout.label_orf(inds) = layout.label(inds);

end

% Final report
fprintf('\n\nNetwork loaded.\n');
fprintf('Number of nodes: %d\n', length(layout.label));
fprintf('Number of edges: %d\n', length(find(layout.edges>0)));
if issymmetric(double(layout.edges))
    fprintf('The edges are not directed.\n');
else
    fprintf('The edges are directed.\n');
end
    
   
    
