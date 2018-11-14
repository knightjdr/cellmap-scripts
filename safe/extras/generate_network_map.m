function layout = generate_network_map(layout, varargin)

%% Check inputs
if mod(nargin,2) ~= 1 
    error('Check your inputs. This function requires a set of property name/property value pairs.');  
end

PROGRESSBAR = 0;
if ~isempty(find(strcmpi('ProgressBar', varargin)))
    PROGRESSBAR = varargin{find(strcmpi('ProgressBar', varargin))+1};
end

%%

if ~isfield(layout,'x') || ~isempty(layout.randomSeed)
    % Assign positions at random
    X = rand([size(layout.label,1) 2]);
else    
    X = [layout.x layout.y];
end
    
if strcmp(layout.layoutAlgorithm, 'Fruchterman-Reingold (beta)')
    
    if PROGRESSBAR
        w = waitbar(0.5, 'Generating the network map... (this make take a few minutes)');
    else
        fprintf('\nGenerating the network map... (this make take a few minutes)\n');
    end
    
    X = fruchterman_reingold_force_directed_layout(sparse(double(layout.edges)), 'edge_weight', 'matrix', 'progressive', X);
    
    if PROGRESSBAR
        delete(w);
    end
    
elseif strcmp(layout.layoutAlgorithm, 'Kamada-Kawai (beta)')
    
    if PROGRESSBAR
        w = waitbar(0.5, 'Generating the network map... (this make take a few minutes)');
    else
        fprintf('\nGenerating the network map... (this make take a few minutes)\n');
    end
    
    X = kamada_kawai_spring_layout(sparse(double(layout.edges)), 'edge_weight', 'matrix', 'progressive', X);
    
    if PROGRESSBAR
        delete(w);
    end
    
elseif strcmp(layout.layoutAlgorithm, 'Kamada-Kawai (Cytoscape)')
    
    X = cyRestKK(double(layout.edges), 'NodePositions', X, 'ProgressBar', PROGRESSBAR);
    
else
    
    error('Layout:unknown','Unknown layout algorithm %s.\nCheck the safe.ini file for the list of available options.', layout.layoutAlgorithm);
    
end

layout.x = X(:,1);
layout.y = X(:,2);


    