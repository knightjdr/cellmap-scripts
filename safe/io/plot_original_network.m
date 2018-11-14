function layout = plot_original_network(layout, varargin)

%% Check inputs

if mod(nargin,2) ~= 1 
    error('Check your inputs. This function requires as inputs the variable |layout| and a set of property name/property value pairs.');  
end

if ~isempty(find(strcmpi('Axes', varargin)))
    ax = varargin{find(strcmpi('Axes', varargin))+1};
end

PROGRESSBAR = 0;
if ~isempty(find(strcmpi('ProgressBar', varargin)))
    PROGRESSBAR = varargin{find(strcmpi('ProgressBar', varargin))+1};
end

if ~isempty(find(strcmpi('Axes', varargin)))
    ax = varargin{find(strcmpi('Axes', varargin))+1};
end

%%

if PROGRESSBAR
    w = waitbar(0, 'Plotting the network...','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    setappdata(w, 'canceling', 0);
else
    fprintf('\nPlotting the network...\n');
end
    
if ~exist('ax')
    figure();
    set(gcf,'Color','k');
    ax = gca;
    ax_tag = 'axes1';
else
    ax_tag = get(ax,'Tag');
    set(0,'CurrentFigure', get(ax,'Parent'));
    set(ax,'NextPlot','replace');
end

if ~isfield(layout,'x')
    % Assign positions at random
    X = rand([size(layout.label,1) 2]);
    layout.x = X(:,1);
    layout.y = X(:,2);
end    
    
plot(layout.x, layout.y, 'w.', 'Parent', ax);
set(ax,'DataAspectRatio',[1 1 1]);
set(ax,'Color','none');
set(ax,'YDir','reverse');
set(ax,'Tag', ax_tag);

axis(ax,'off');
hold on;

[r,c] = find(triu(layout.edges) > 0);
for i = 1 : length(r)
    
     % Check for Cancel button press
    if PROGRESSBAR && getappdata(w,'canceling')
        break;
    end
    
    if PROGRESSBAR && mod(i, 1000) == 0
        waitbar(i/length(r), w);
    end
    
    xs = [layout.x(r(i)), layout.x(c(i))];
    ys = [layout.y(r(i)) layout.y(c(i))];
    l = patchline(xs, ys, 'EdgeColor', 'w', 'LineWidth', 0.5, 'EdgeAlpha', 0.2, 'Parent', ax);
end

if PROGRESSBAR
    delete(w);
end
