function layout = cluster_groups(layout, varargin)

%% Check inputs

if mod(nargin,2) ~= 1 
    error('Check your inputs. This function requires as inputs the variable |layout| and a set of property name/property value pairs.');  
end

% Show GUI-style progress bar
PROGRESSBAR = 0;
if ~isempty(find(strcmpi('ProgressBar', varargin)))
    PROGRESSBAR = varargin{find(strcmpi('ProgressBar', varargin))+1};   
end

%%

SGN = size(layout.pval,3);
GRP = length(layout.group_ids);

if strcmp(layout.annotationsign,'both')
    SGN_NAME = {'highest','lowest'};
else
    SGN_NAME = {layout.annotationsign};
end

layout.groupColor = zeros(GRP, SGN);

if PROGRESSBAR
    w = waitbar(PROGRESSBAR, 'Computing attribute similarity...');
else
    fprintf('Computing attribute similarity...\n');
end

for sgn = 1 : SGN

    GRPT = numel(find(layout.groupIsTop(:,sgn)));
    
    if ~layout.groupDistance || GRPT < 2
        
        layout.groupColor(layout.groupIsTop(:,sgn),sgn) = [1:GRPT]'; 

        if layout.groupDistance
            
            if GRP < 2
                fprintf('\nValues %s: You have less than 2 attributes. Skipping attribute grouping.\n', SGN_NAME{sgn});
            elseif GRPT < GRP
                fprintf('\nValues %s: Of your %d attributes, less than 2 make the criteria to be considered for grouping. Skipping attribute grouping.\n', SGN_NAME{sgn}, GRP);
                fprintf('Next time, you may consider changing your safe.ini settings to:\n');
                if layout.unimodality
                    fprintf('\t- set the unimodalityType option to empty\n');
                end
                fprintf('\t- decrease the groupEnrichmentMinSize threshold\n');
                fprintf('\t- increasing the neighborhoodRadius threshold\n');
            end
        end
        
    else

        switch layout.groupDistanceType
            case 'jaccard'
                dist = pdist(double(layout.opacity_01(:,layout.groupIsTop(:,sgn),sgn))','jaccard');
            case 'correlation'
                dist = pdist(layout.opacity(:,layout.groupIsTop(:,sgn), sgn)','correlation');
            case 'imed'
                %     dist = imed2(layout.opacity(:,layout.groupIsTop) .* layout.opacity_01(:,layout.groupIsTop),layout.shortpath_weighted);
                %
                %     dist = dist';
                %     dist = dist(:);
                %     dist(isnan(dist)) = [];
                %     dist = dist';
                %
                %     dist_total = dist./(1-dist_jaccard);
                %     dist_total(isinf(dist_total)) = max(dist_total(~isinf(dist_total)));
                %
                %     % dist_square = squareform(dist);
                %     % dist_jaccard_square = squareform(dist_jaccard);
                %     % dist_total_square = squareform(dist_total);
                %     dist = dist_total;
            otherwise
                error('Unknown |layout.groupDistanceType|. Check the list of available options in the settings file.');
        end

        Z = linkage(dist,'average');

        % groupOrder = optimalleaforder(Z, dist);
        % 
        % [~,ix] = sort(groupOrder);
        % layout.groupOrder = ones(size(layout.group_ids));
        % layout.groupOrder(layout.groupIsTop) = 1+ix;

        % Group colors = group clusters
        layout.groupColor(layout.groupIsTop(:,sgn),sgn) = cluster(Z,'cutoff',layout.groupDistanceThreshold*max(Z(:,3)),'criterion','distance');
        
    end

    if PROGRESSBAR
        delete(w);
    end 
end

