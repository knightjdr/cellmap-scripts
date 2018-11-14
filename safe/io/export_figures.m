function export_figures(layout, varargin)

%% Check inputs

axs = gca;
if ~isempty(find(strcmpi('Axes', varargin)))
    axs = varargin{find(strcmpi('Axes', varargin))+1};
end

if ~isfield(layout, 'outputdir')
    p = mfilename('fullpath');
    ind = regexp(p, '/io/load_network');
    path_to_results_folder = p(1:ind);
    
    [pathstr, ~, ~] = fileparts(path_to_results_folder);
    layout.outputdir = [pathstr '/safe-' datestr(now,'yyyy-mm-dd-HH-MM-SS') '/'];  
end

OUTPUTDIR = layout.outputdir;

if ~exist(OUTPUTDIR, 'dir')
    mkdir(OUTPUTDIR);
end

PROGRESSBAR = 0;
if ~isempty(find(strcmpi('ProgressBar', varargin)))
    PROGRESSBAR = varargin{find(strcmpi('ProgressBar', varargin))+1};
end

%% Saving figures to PDF

% axs_names = {'original_network','sample_attribute','composite_map'};

if PROGRESSBAR
    w = waitbar(0, 'Exporting figures to PDF...');
else
    fprintf('\nExporting figures to PDF...\n');
end

for i = 1 : length(axs)
    
    ax_id = textscan(get(axs(i),'Tag'),'axes%d');
    ax_id = ax_id{1};
    
    if isempty(ax_id)
        ax_id = 10+i;
    end
    

    fig2 = figure('visible','off');
    newax = copyobj(axs(i),fig2);
    set(newax, 'units', 'normalized', 'position', [0.13 0.11 0.775 0.815]);

    % rgb_colors = uint8(get(scatterPoints, 'Cdata')*255);
    % rgb_colors = reshape(rgb_colors, size(rgb_colors,1),1,size(rgb_colors,2));
    % [i,map] = rgb2ind(rgb_colors, size(map_colors,1)*2);
    % set(scatterPoints, 'Cdata', double(i));
    % colormap(map);
    
    filename = [layout.outputdir 'figure' num2str(ax_id) '.pdf'];
    
    n = 1;
    while exist(filename)
        filename = [layout.outputdir 'figure' num2str(ax_id) '-' num2str(n) '.pdf'];
        n = n+1;
    end
        
    print(fig2, filename, '-painters', '-dpdf', '-r600', '-noui');

    close(fig2);
    
    if PROGRESSBAR
        waitbar(i/length(axs), w);
    end
    
end

if PROGRESSBAR
    delete(w);
end

    
% NOTE: as of Feb. 4, 2015 saving to PDF in Matlab 2014b doesn't work very well.
% The nodes appear as hexagons, instead of circles, and it just doesn't
% look good. To save a publication quality figure, use Matlab 2014a or
% earlier. Alternatively, open the saved PDF in Adobe Illustrator, select all nodes 
% and choose Effect > Convert to shape > Ellipse. Set the absolute width and height of the
% nodes to 4 px or similar.