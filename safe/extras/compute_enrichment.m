function layout = compute_enrichment(layout, varargin)

if mod(nargin,2) ~= 1 
    error('Check your inputs. This function requires as inputs the variable |layout| and a set of property name/property value pairs.');  
end

% Show GUI-style progress bar
PROGRESSBAR = 0;
if ~isempty(find(strcmpi('ProgressBar', varargin)))
    PROGRESSBAR = varargin{find(strcmpi('ProgressBar', varargin))+1};   
end

BACKGROUND = 'map';
if isfield(layout, 'background') && ~isempty(layout.background)
    BACKGROUND = layout.background;
end

%% Defining the neighborhood radius

inds = find(~isinf(layout.nodeDistance) & ~isnan(layout.nodeDistance));

if ~isfield(layout, 'neighborhoodRadiusType') || isempty(layout.neighborhoodRadiusType) || strcmp(layout.neighborhoodRadiusType, 'percentile')
    layout.R = prctile(layout.nodeDistance(inds), [layout.neighborhoodRadius]);
elseif strcmp(layout.neighborhoodRadiusType, 'diameter')
    layout.R = layout.neighborhoodRadius * (max(layout.x) - min(layout.x)) / 100;
else
    layout.R = layout.neighborhoodRadius;
end

%% Defining neighborhoods

% tau = 0.015;
% diam = max(layout.nodeDistance(inds));
% layout.neighborhoods = exp(-layout.nodeDistance/(tau*diam));
layout.neighborhoods = layout.nodeDistance <= layout.R;

%% Compute enrichments

NLBL = length(layout.label);
NGRP = length(layout.group_ids);

if isempty(BACKGROUND) || strcmp('map', BACKGROUND)
    % Total number of nodes on the map (in matrix format):
    N = zeros(NLBL, NGRP) + NLBL;  

    % Number of nodes (on the map) annotated to this attribute
    Ng = repmat(nansum(layout.label2group,1), NLBL,1);
else
    % Total number of nodes in the standard (in matrix format):
    N = zeros(NLBL, NGRP) + layout.group_total_gene_number;

    % Number of nodes (in the standard) annotated to this attribute
    Ng = repmat(layout.group_global_size', length(layout.label),1);

end

% Number of nodes in each node's neighborhood
Ni = repmat(sum(layout.neighborhoods,2), 1, length(layout.group_ids));

% Number of nodes in each node's neighborhood annotated to each attribute
Nig = layout.neighborhoods * layout.label2group;

layout.fld = (Nig./Ni) ./ (Ng./N);

ixNotNaN = find(~isnan(layout.label2group(:)));

if isempty(find(~ismember(layout.label2group(ixNotNaN), [0 1])))
        
    fprintf('\nThe annotation standard is binary. Using the standard hypergeometric enrichment schema.\n');

    % Some networks are too big to run hypergeometric analysis on the
    % entire matrix (Matlab runs out of memory).
    % Simple solution: run it on a attribute by attribute basis.

    layout.pval = nan(size(layout.fld));
    
    if PROGRESSBAR
        w = waitbar(0, 'Calculating the hypergeometric enrichment p-value...','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
        setappdata(w, 'canceling', 0);
    else
        fprintf('\nCalculating the hypergeometric enrichment p-value...\n');
        fprintf(['|', blanks(100), '|\n']);
        fprintf('|');
    end

    y = 0;

    for gr = 1 : NGRP
        
        % Check for Cancel button press
        if PROGRESSBAR && getappdata(w,'canceling')
            delete(w);
            error('Canceled.');
        end

        % -- Progress bar ---
        if PROGRESSBAR
            waitbar(gr / NGRP,w);
        else
            x = round(gr * 100 / NGRP);
            if x > y
                fprintf(repmat('*',1,x-y)); y = x;
            end
        end
        % -- Progress bar ---

        layout.pval(:,gr) = 1-hygecdf(Nig(:,gr)-1, N(:,gr), Ng(:,gr), Ni(:,gr));

    end

    if PROGRESSBAR
        delete(w);
    else
        fprintf('|\n');
    end
    
else
    
    fprintf('\nThe annotation standard is quantitative. Using a more general permutation approach to calculate enrichment.\n');
    
    nPermutations = 1000;

    Sr = nan(NLBL,NGRP,nPermutations);
    
    if PROGRESSBAR
        w = waitbar(0, 'Running random permutations...','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
        setappdata(w, 'canceling', 0);
    else
        fprintf('\nRunning random permutations to calculate significance...\n');
        fprintf(['|', blanks(100), '|\n']);
        fprintf('|');
    end

    y = 0;
    
    for r = 1 : nPermutations
        
        % Check for Cancel button press
        if PROGRESSBAR && getappdata(w,'canceling')
            delete(w);
            error('Canceled.');
        end
        
        % -- Progress bar ---
        if PROGRESSBAR
            waitbar(r / nPermutations, w);
        else
            x = round(r * 100 / nPermutations);
            if x > y
                fprintf(repmat('*',1,x-y)); y = x;
            end
        end
        % -- Progress bar ---
        
        ixPerm = randperm(NLBL);
        Wr = layout.label2group(ixPerm,:);

        Sr(:,:,r) = layout.neighborhoods * Wr;
    end
    
    if PROGRESSBAR
        delete(w);
    else
        fprintf('|\n');
    end
    
    
    layout.pval = nan(NLBL,NGRP);

    if PROGRESSBAR
        w = waitbar(0, 'Calculating the empirical enrichment p-value...','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
        setappdata(w, 'canceling', 0);
    else
        fprintf('\nCalculating the empirical enrichment p-value...\n');
        fprintf(['|', blanks(100), '|\n']);
        fprintf('|');
    end

    y = 0;

    for grp = 1 : NGRP
        
        % Check for Cancel button press
        if PROGRESSBAR && getappdata(w,'canceling')
            delete(w);
            error('Canceled.');
        end

        % -- Progress bar ---
        if PROGRESSBAR
            waitbar(grp / NGRP, w);
        else
            x = round(grp * 100 / NGRP);
            if x > y
                fprintf(repmat('*',1,x-y)); y = x;
            end
        end
        % -- Progress bar ---
        
        Sm = nanmean(Sr(:,grp,:),3);
        Ss = nanstd(Sr(:,grp,:),[],3);

        Z = (Nig(:,grp) - Sm)./Ss;
            
        switch layout.annotationsign
            case 'highest'
                layout.pval(:,grp) = 1-normcdf(Z);
%                 layout.pval(:,grp) = sum(Sr(:,grp,:) >= repmat(Nig(:,grp),1,1,1000),3)/nPermutations;
            case 'lowest'
                layout.pval(:,grp) = normcdf(Z);
%                 layout.pval(:,grp) = sum(Sr(:,grp,:) <= repmat(Nig(:,grp),1,1,1000),3)/nPermutations;
            otherwise
                layout.pval(:,grp,1) = 1-normcdf(Z);
                layout.pval(:,grp,2) = normcdf(Z);                
%                 layout.pval(:,grp,1) = sum(Sr(:,grp,:) >= repmat(Nig(:,grp),1,1,1000),3)/nPermutations;
%                 layout.pval(:,grp,2) = sum(Sr(:,grp,:) <= repmat(Nig(:,grp),1,1,1000),3)/nPermutations;

        end
                

    end

    if PROGRESSBAR
        delete(w);
    else
        fprintf('|\n');
    end
    
end

% Transforming enrichment p-values into scores

layout.opacity = -log10(layout.pval);

% This version was used in older versions of SAFE:
% m = max(max(layout.opacity(~isinf(layout.opacity))));

m = layout.MAX_LOG10_PVAL;

layout.opacity(layout.opacity > m) = m;
layout.opacity = layout.opacity ./ m;

% Calculate the minimum opacity corresponding to significant enrichment (after Bonferroni multiple testing correction)
layout.thresholdOpacity = -log10(layout.THRESHOLD_ENRICHMENT/length(layout.group_ids))/layout.MAX_LOG10_PVAL;    

layout.opacity_01 = layout.opacity > layout.thresholdOpacity;

layout.numLabelsEnrichedGroup = permute(sum(layout.opacity_01,1),[2 3 1]);
layout.numLabelsGroup = permute(nansum(layout.label2group,1), [2 3 1]);

layout.numGroupsEnrichedLabel = permute(sum(layout.opacity_01,2),[1 3 2]);

