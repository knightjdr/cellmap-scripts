function layout = calculate_node_distances(layout, varargin)

%% Check inputs

if mod(nargin,2) ~= 1 
    error('Check your inputs. This function requires a set of property name/property value pairs.');  
end

% Show GUI-style progress bar
PROGRESSBAR = 0;
if ~isempty(find(strcmpi('ProgressBar', varargin)))
    PROGRESSBAR = varargin{find(strcmpi('ProgressBar', varargin))+1};   
end

%% Calculate distances

if PROGRESSBAR
    w = waitbar(0.5, 'Calculating node distances...');
end
    
switch layout.nodeDistanceType
    case 'shortpath'
        layout.nodeDistance = graphallshortestpaths(sparse(double(layout.edges)));
    case 'shortpath_weighted_layout'   
        [x1,x2] = meshgrid(layout.x);
        [y1,y2] = meshgrid(layout.y);
        layout.nodeDistanceEuclid = sqrt( (x1-x2).^2 + (y1-y2).^2 ); % Physical distances between all pairs of nodes
        clear x* y*;
        layout.nodeDistance = graphallshortestpaths(sparse(double(layout.edges)), 'Weights', reshape(layout.nodeDistanceEuclid(layout.edges>0),[],1));
    case 'shortpath_weighted_edge'
        layout.nodeDistance = graphallshortestpaths(sparse(double(layout.edges)), 'Weights', double(reshape(1-layout.edges_weights(layout.edges>0),[],1)));
    otherwise
        error(['Unknown type of node distance: ' layout.nodeDistanceType '. Check the list of available options in the settings.m file.']);
end

if PROGRESSBAR
    delete(w);
end