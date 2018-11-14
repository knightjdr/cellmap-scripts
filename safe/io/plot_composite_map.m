function [layout, varargout] = plot_composite_map(layout, varargin)

%% Check inputs

if mod(nargin,2) ~= 1 
    error('Check your inputs. This function requires as inputs the variable |layout| and a set of property name/property value pairs.');  
end

if ~isempty(find(strcmpi('Axes', varargin)))
    ax = varargin{find(strcmpi('Axes', varargin))+1};
end

annotationsign = 'highest';
if ~isempty(find(strcmpi('AnnotationSign', varargin)))
    annotationsign = varargin{find(strcmpi('AnnotationSign', varargin))+1};
end

%%

if strcmp('both', layout.annotationsign)
    switch annotationsign
        case 'highest'
            sgn = 1;
        case 'lowest'
            sgn = 2;
        otherwise
            error('Annotation sign unknown.');
    end
else
    annotationsign = layout.annotationsign;
    sgn = 1;
end
    
if ~exist('ax')
    figure();
    set(gcf,'Color','k');
    ax = gca;
    set(ax,'Units','pixels');
    ax_tag = 'axes3';
else
    ax_tag = get(ax, 'Tag');
    set(0,'CurrentFigure', get(ax,'Parent'));
    set(ax,'NextPlot','replace');
end

% Final list of colors, each corresponding to a functional region.
ucolors = unique(layout.groupColor(:,sgn));
if ucolors(1) ~= 1
    ucolors = [1; ucolors];
end

% Colors. 
map_colors = colormap(hsv(length(ucolors)-1)); 
map_colors = map_colors(randperm(size(map_colors,1)),:);

% Region #1 (containing nodes with no enrichment or nodes that didn't make it into one of the final regions) should be black.
map_colors = [0 0 0; map_colors]; 

layout.mapColors = map_colors;

% The final color of every node
layout.geneRGB = zeros(length(layout.label),3);
for i = 1 : length(layout.label)

    inds = find(layout.opacity_01(i,:,sgn) .* layout.groupIsTop(:,sgn)');

    if ~isempty(inds)
        layout.geneRGB(i,:) = sum(map_colors(layout.groupColor(inds,sgn),:) .* repmat(layout.opacity(i,inds,sgn)'.^2,1,3),1)/length(inds);
    end

end

%     layout.geneRGB = map_colors(layout.labelColor,:);


% Adjust the brightness
layout.geneRGB = layout.geneRGB .* 1.5;
layout.geneRGB(layout.geneRGB > 1) = 1;

[r, ~] = find(layout.geneRGB > 1);
r = unique(r);
layout.geneRGB(r,:) = layout.geneRGB(r,:) ./ repmat(max(layout.geneRGB(r,:),[],2),1,3);

% Sort nodes by brightness. 
tmp = sum(layout.geneRGB,2);
[~,inds] = sort(tmp,'ascend');  

%     f(sgn) = figure();
%     set(f(sgn),'Color','k');

scatterPoints = scatter(ax, layout.x(inds), layout.y(inds), 120, layout.geneRGB(inds,:),'filled');

set(ax,'DataAspectRatio',[1 1 1]);
set(ax,'YDir','reverse');
set(ax,'Color','none');
set(ax,'TitleFontSizeMultiplier',1,'TitleFontWeight','normal');
set(ax,'Tag', ax_tag);

axis(ax, 'off');

hold(ax, 'on');

%     % Plot legend for radius
%     minx = min(layout.x(inds));
%     maxy = max(layout.y(inds));
% 
%     theta = linspace(0,2*pi);
%     xc = layout.R*cos(theta)+minx;
%     yc = layout.R*sin(theta)+maxy;
% 
%     hold on;
%     plot(xc,yc,'w--');


% Contour
k = convhull(layout.x, layout.y);
plot(layout.x(k), layout.y(k), 'w--', 'Parent', ax);

% Title
ttl = sprintf('Composite map -- %s values', annotationsign);
t = title(ax, ttl,'FontSize',12,'Color','w','Interpreter','none', 'Units','pixels');
Pax = get(ax,'Position');
Pt = get(t,'Position');
set(t,'Position', [Pt(1) Pax(4) Pt(3)],'VerticalAlignment','middle');

varargout{1} = scatterPoints; 

