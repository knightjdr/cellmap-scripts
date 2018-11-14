function hText = plot_labels(layout, varargin)

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
    ax = gca;
else
    set(0,'CurrentFigure', get(ax,'Parent'));
    set(ax,'NextPlot','replace');
end

ucolors = unique(layout.groupColor(:,sgn));
ucolors(ucolors==1) = [];

hText = cell(size(ucolors));

for i = 1 : length(ucolors)
%             inds = find(layout.labelColor(:,sgn) == ucolors(i));
    indsg = find(layout.groupColor(:,sgn) == ucolors(i));
    [r,~,~] = find(layout.opacity_01(:,indsg,sgn)>0);
    inds = unique(r);
    centroid_x = sum((layout.labelOpacity(inds,sgn).^2) .* layout.x(inds)) ./ sum((layout.labelOpacity(inds,sgn).^2));
    centroid_y = sum((layout.labelOpacity(inds,sgn).^2) .* layout.y(inds)) ./ sum((layout.labelOpacity(inds,sgn).^2));
    hold on;
    hText{i} = text(centroid_x, centroid_y, sprintf('%d',ucolors(i)),'FontSize',20,'Color','w','HorizontalAlignment','center', 'Parent', ax);
end
    