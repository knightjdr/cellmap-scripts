function layout = compute_unimodality(layout, varargin)

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

if ~layout.unimodality
    
    % If the unimodality setting is OFF, all attributes are to be considered.
    layout.groupIsTop = true(size(layout.pval,2), size(layout.pval,3));

else
    
    % If the unimodality setting in ON, identify the region-specific attributes.

    GRP = numel(layout.group_ids);
    SGN = size(layout.opacity,3);

    layout.groupUnimode = nan(GRP,SGN);

    inds_ref = find(~isinf(layout.nodeDistance) & ~isnan(layout.nodeDistance));

    if PROGRESSBAR
        w = waitbar(0, 'Calculating a unimode coefficient for every attribute...', 'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
        setappdata(w, 'canceling', 0);
    else    
        fprintf('Calculating a unimode coefficient for every attribute...\n');
        fprintf(['|', blanks(100), '|\n']);
        fprintf('|');
    end

    y = 0;

    for grp = 1 : GRP
        
         % Check for Cancel button press
        if PROGRESSBAR && getappdata(w,'canceling')
            break;
        end

        % -- Progress bar ---
        if PROGRESSBAR
            waitbar(grp / GRP, w);
        else
            x = round(grp * 100 / GRP);
            if x > y
                fprintf(repmat('*',1,x-y)); y = x;
            end
        end
        % -- Progress bar ---

        for sgn = 1 : SGN
            inds = find(layout.opacity_01(:,grp,sgn)>0);

            if length(inds)>=5

                d = sort(reshape(layout.nodeDistance(inds,inds),[],1));

                switch layout.unimodalityType
                    case 'hartigans'    
                        [dip, p_value, xlow, xup] = HartigansDipSignifTest(d, 500); 
                        layout.groupUnimode(grp,sgn) = p_value;
                    case 'radius' 
                        r = prctile(d,[65]);
                        layout.groupUnimode(grp,sgn) = (r <= 2*layout.R);
                    case 'subclust'
                        bounds = [min(layout.x) min(layout.y); max(layout.x) max(layout.y)];
                        [C,~] = subclust([layout.x(inds), layout.y(inds)],0.15, bounds);
                        layout.groupUnimode(grp,sgn) = (size(C,1)==1);
                    otherwise
                        warning('Unknown unimodality type setting. Check the list of available options in the settings.m file.');
                end

            end
        end

    end
    
    if PROGRESSBAR
        delete(w);
    else
        fprintf('|\n');
    end

    layout.groupIsTop = layout.numLabelsEnrichedGroup > layout.groupEnrichmentMinSize & layout.groupUnimode > 0;
    
end
