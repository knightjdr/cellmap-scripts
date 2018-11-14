function plot_sample_groups(layout, varargin)

%% Check inputs

if mod(nargin,2) ~= 1 
    error('Check your inputs. This function requires as inputs the variable |layout| and a set of property name/property value pairs.');  
end

if nargin == 1
    varargin(end+1) = {'Random'};
    varargin(end+1) = {10};
end

if ~isempty(find(strcmpi('Random', varargin))) & ~isempty(find(strcmpi('Range', varargin)))
    error('Check your inputs. Specify either a number of random groups (''Random'') or a specific range of groups (''Range'') to plot.');  
end

if ~isempty(find(strcmpi('Random', varargin)))
    NUMRANDOM = varargin{find(strcmpi('Random', varargin))+1}; 
    
    % Only consider "top" attributes (i.e., "region-specific" if the unimodality setting is ON or "all" if the unimodality setting is OFF)
    indsTop = find(layout.groupIsTop);
    NUMPLOTS = min(length(indsTop), NUMRANDOM);
    GRPS = randsample(indsTop, NUMPLOTS);
    
    fprintf('\nPlotting a set of %d random attributes...\n', NUMPLOTS);
end

if ~isempty(find(strcmpi('Range', varargin)))
    GRPS = varargin{find(strcmpi('Range', varargin))+1};
    fprintf('\nPlotting the specified set of attributes...\n');
end

if ~isempty(find(strcmpi('Axes', varargin)))
    ax = varargin{find(strcmpi('Axes', varargin))+1};
end

%%

range = [0:0.25:1]; % The opacity steps for the figure legend.

for i = 1 : length(GRPS)
    
    if ~exist('ax')
        figure();
        set(gcf,'Color','k');
        ax = gca;
        set(ax,'Units','pixels');
        ax_tag = 'axes2';
    else
        ax_tag = get(ax,'Tag');
        set(0,'CurrentFigure', get(ax,'Parent'));
        set(ax,'NextPlot','replace');
    end
    
    if strcmp('both', layout.annotationsign)
        
        % Bicolor
        
        o = cat(3, layout.opacity(:,:,1), layout.opacity(:,:,2));
        o = max(o,[],3);
        
        inds = find(~isnan(o(:,GRPS(i))));
        [~,ix] = sort(o(inds,GRPS(i)),'ascend');
        
        map_colors = [255 204 0; 0 204 255]/255;    % yellow/blue
%         map_colors = [0 255 0; 255 0 0]/255;      % green/red
        
        c1 = repmat(map_colors(1,:),length(ix),1) .* repmat(layout.opacity(inds(ix),GRPS(i),1),1,3);
        c2 = repmat(map_colors(2,:),length(ix),1) .* repmat(layout.opacity(inds(ix),GRPS(i),2),1,3);
        
        c = cat(3,c1,c2);
        c = nansum(c,3);
        
        % Color legend
        l1 = repmat(map_colors(1,:), length(range),1) .* repmat(range',1,3);
        l2 = repmat(map_colors(2,:), length(range),1) .* repmat(range',1,3);
        
        l = [flipud(l1); l2(2:end,:)];
        
    else
    
        inds = find(~isnan(layout.opacity(:,GRPS(i))));
        [~,ix] = sort(layout.opacity(inds,GRPS(i)),'ascend');

        c = repmat(layout.opacity(inds(ix),GRPS(i)),1,3);
        
        l = repmat(range',1,3);
    
    end
    
    scatterPoints = scatter(ax, layout.x(inds(ix)), layout.y(inds(ix)), 120, c,'filled');
    
    set(ax,'DataAspectRatio',[1 1 1]);
    set(ax,'YDir','reverse');
    set(ax,'Color','none');
    if isprop(ax,'TitleFontSizeMultiplier')
        set(ax,'TitleFontSizeMultiplier',1);
    end
    if isprop(ax, 'TitleFontWeight')
        set(ax, 'TitleFontWeight','normal');
    end
    set(ax,'Tag', ax_tag);
    
    axis(ax, 'off');
    
    hold(ax, 'on');

    % Plot contour
    k = convhull(layout.x, layout.y, 'simplify', true);
    plot(layout.x(k), layout.y(k), 'w--', 'Parent', ax);
    
    
    % Plot legend for radius
    if strcmp(layout.nodeDistanceType, 'shortpath_weighted_layout')
        maxx = max(layout.x);
        maxy = max(layout.y);

        theta = linspace(0,2*pi);
        xc = layout.R*cos(theta)+maxx;
        yc = layout.R*sin(theta)+maxy-layout.R;

        plot(xc,yc,'w:','Parent', ax);
        text(mean(xc),mean(yc),{'Max Distance';'Threshold';'(approx.)'},'Color', 'w', 'FontSize', 10, 'HorizontalAlignment','center', 'Parent', ax);
    end
    
    % Plot legend on color  
    minx = min(layout.x)*0.8;
    maxx = max(layout.x)*0.8;
    miny = min(layout.y);
    maxy = max(layout.y)*1.1;
    
    lx = linspace(minx, maxx, size(l,1));
    ly = repmat(maxy, length(lx), 1);
    s = scatter(ax, lx, ly, 120, l, 'filled');
    
    TeXString = {};
    for r = 1 : length(range)
        TeXString{r} = texlabel(['10^(-' num2str(range(r)*16) ')']);
    end
    
    shft = 0.25*(maxx-minx)/size(l,1);

    if strcmp('both', layout.annotationsign)
        TeXString = [flipud(TeXString'); TeXString(2:end)'];
        shft = shft * 2;
    end
    
    text(lx,ly+shft, TeXString, 'Color', 'w', 'FontSize', 10, 'HorizontalAlignment', 'center', 'Parent', ax);
    text(mean(lx), mean(ly)+2*shft, 'Enrichment p-value','Color','w','FontSize', 10, 'HorizontalAlignment', 'center','Parent', ax);
    
    if strcmp('both', layout.annotationsign)
        text(min(lx),min(ly)-shft, 'Highest Values', 'Color','w', 'FontSize', 10, 'HorizontalAlignment', 'center', 'Parent', ax);
        text(max(lx),max(ly)-shft, 'Lowest Values', 'Color','w', 'FontSize', 10, 'HorizontalAlignment', 'center', 'Parent', ax);
    end

    % Title
    ttl = sprintf('Attribute %d -- %s', layout.group_ids(GRPS(i)), layout.group_names{GRPS(i)});
    t = title(ax, ttl,'FontSize',12,'Color','w','Interpreter','none');
%     Pax = get(ax,'Position');
%     Pt = get(t,'Position');
%     set(t,'Position', [Pax(3)/2 Pax(4) Pt(3)],'VerticalAlignment','middle');
%     set(t,'Position', [Pt(1) miny Pt(3)],'VerticalAlignment','middle');
    
%     % Plot cluster markers
%     Bounds = [min(layout.x) min(layout.y); max(layout.x) max(layout.y)];
%     Center = [nanmean(layout.x) nanmean(layout.y)];
%     Radius = max(sqrt((Center(1)-layout.x).^2 + (Center(2)-layout.y).^2))*1.05;
% 
%     inds = find(layout.opacity_01(:,GRPS(i))>0);
%     
%     if ~isempty(inds)
%         [C,S] = subclust([layout.x(inds), layout.y(inds)],0.1,Bounds);
% 
%         for j = 1 : size(C,1)
%     %         plot(C(:,1),C(:,2),'g*','MarkerSize',20);
% 
%             n = length(find(sqrt((layout.x(inds)-C(j,1)).^2 + (layout.y(inds)-C(j,2)).^2) < nanmean(S)));
% 
%             if n > 3
% 
%                 slope = (C(j,2)-Center(2))/(C(j,1)-Center(1));
%                 intercpt = -slope*Center(1)+Center(2);
% 
%                 [xout,yout] = linecirc(slope, intercpt, Center(1), Center(2), Radius);
% 
%                 [~,imin] = min(sqrt((C(j,1)-xout).^2 + (C(j,2)-yout).^2));
%                 plot([C(j,1) xout(imin)],[C(j,2) yout(imin)],'w-');
%             end
%         end
%     end
    
    
% % --- Uncomment this to save the plots to individual files ----
%     % Save as PNG
%     set(gcf, 'InvertHardCopy', 'off');
%     print('-painters','-dpng', '-r100', [ttl '.png']);
% 
%     rgb_colors = uint8(get(scatterPoints, 'Cdata')*255);
%     rgb_colors = reshape(rgb_colors, size(rgb_colors,1),1,size(rgb_colors,2));
%     [k,map] = rgb2ind(rgb_colors, 100);
%     set(scatterPoints, 'Cdata', double(k));
%     colormap(map);
% % ---

    clear ax;
end