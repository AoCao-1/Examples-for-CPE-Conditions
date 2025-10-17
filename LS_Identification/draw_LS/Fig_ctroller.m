clear; clc; close all;
load('u_MCPE_all.mat')
load('u_CCPE_all.mat')
load('u_HCPE_all.mat')

controllers = {u_MCPE_all, u_CCPE_all, u_HCPE_all};
names = {'MCPE','CCPE','HCPE'};
numMethods = numel(controllers);

% Known Hankel matrix column lengths
HankelCols = [61, 21, 43];

% Color and line style configuration (deep blue + gray)
baseColor = [0 0.2 0.6];
altColor = [0.5 0.5 0.5];
styles = {'-','--','-.'};

numSegments = max(cellfun(@numel, controllers));
% Alternate colors so that segments with the same index have consistent colors across subplots
colors = repmat(baseColor, numSegments, 1);
for i = 1:numSegments
    if mod(i,2)==0
        colors(i,:) = altColor;
    end
end

% Figure size and font settings
figure('Units','centimeters','Color','w');
fontMain = 'Times New Roman';
set(groot,'defaultAxesFontName',fontMain,'defaultTextFontName',fontMain);

for c = 1:numMethods
    u_all = controllers{c};
    ax = subplot(numMethods,1,c);
    hold(ax,'on');
    % Title (using \textemdash to avoid LaTeX errors)
    title(ax, ['\textbf{' names{c} ' \textemdash\ Control Input (1st Dim)}'], ...
        'Interpreter','latex','FontSize',13);
    xlabel(ax, '\textbf{Time step (concatenated)}', 'Interpreter','latex','FontSize',11);
    ylabel(ax, '\textbf{Control input}', 'Interpreter','latex','FontSize',11);

    offset = 0;
    y_min = inf; y_max = -inf;
    mids = []; Ts = []; colorIdx = [];
    styleIdxs = zeros(1,numel(u_all));
    
    % Plot the curve and segment boundaries; record midpoints and lengths for labeling
    for s = 1:numel(u_all)
        seg = u_all{s};
        if isempty(seg)
            continue;
        end
        [r,cseg] = size(seg);
        if r == m
            u1 = seg(1,:); Tseg = cseg;
        elseif cseg == m
            u1 = seg(:,1)'; Tseg = r;
        else
            error('Invalid data dimension in method %d, segment %d', c, s);
        end

        tseg = (1:Tseg) + offset;
        style = styles{mod(s-1,numel(styles))+1};
        plot(ax, tseg, u1, style, 'Color', colors(s,:), 'LineWidth', 1.5);

        % Segment boundary line
        xline(ax, offset + Tseg, ':', 'Color', [0.78 0.78 0.78], 'LineWidth', 0.9);

        % Record information for labeling
        mids(end+1) = offset + floor(Tseg/2); %#ok<SAGROW>
        Ts(end+1) = Tseg; %#ok<SAGROW>
        colorIdx(end+1) = s; %#ok<SAGROW>
        styleIdxs(s) = mod(s-1,numel(styles))+1;

        y_min = min(y_min, min(u1));
        y_max = max(y_max, max(u1));

        offset = offset + Tseg;
    end

    % Adjust x-axis limits to keep the last segment visible with some margin
    xlim(ax, [0, offset + 1]);

    % Adjust y-axis limits with padding for label placement
    pad = 0.20*(abs(y_max - y_min) + eps);
    ylim(ax, [-2.15, y_max + pad]);

    % Add LaTeX annotations for each segment (ensuring T_10 also displays)
    for k = 1:length(mids)
        y_text = y_max + 0.06*(abs(y_max - y_min) + eps);
        text(ax, mids(k), y_text, sprintf('$T_{%d} = %d$', k, Ts(k)), ...
            'Interpreter','latex', 'HorizontalAlignment','center', ...
            'FontSize',12, 'Color', colors(colorIdx(k),:), 'FontWeight','bold');
    end

    % Add Length/Hankel information box in the top-right corner
    xl = xlim(ax); yl = ylim(ax);
    total_len = offset;
    box_str = sprintf('$\\mathbf{Length} = %d,\\ \\mathbf{Hankel\\ cols} = %d$', total_len, HankelCols(c));
    tb = text(ax, xl(2), yl(2), box_str, 'Interpreter','latex', ...
        'HorizontalAlignment','right','VerticalAlignment','top', ...
        'FontSize',11, 'FontWeight','bold', 'BackgroundColor','white', 'EdgeColor','none', 'Margin',8);
    % If needed for older MATLAB versions, tb.Position can be slightly adjusted manually.

    % Beautify axes: remove top/right border, outward ticks, light grid
    ax.Box = 'off';
    ax.TickDir = 'out';
    ax.XGrid = 'on'; ax.YGrid = 'on';
    ax.GridColor = [0.92 0.92 0.92];
    ax.GridAlpha = 0.4;

    hold(ax,'off');
end

% Global font settings
set(findall(gcf,'-property','FontName'),'FontName',fontMain);
set(findall(gcf,'-property','FontSize'),'FontSize',12);

% Optional: export as high-quality vector graphic
% exportgraphics(gcf,'control_inputs_comparison_boxed.pdf','ContentType','vector','BackgroundColor','none');
