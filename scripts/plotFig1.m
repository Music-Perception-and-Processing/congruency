clear
close all

set(0,'defaultTextInterpreter','latex')
set(0,'defaultAxesTickLabelInterpreter','latex')
set(0,'defaultLegendInterpreter','latex')


load('data/expData_1.mat')
load('data/analysis_v2.mat')
load('data/RangeIdx.mat')
load('data/NoteData.mat')
load('data/pca_data_v5.mat')
run plot_properties.m

figh1 = figure('Visible',pp.visible,...
    'DefaultTextFontName',pp.fname,...
    'DefaultAxesFontName',pp.fname);
axh1 = axes('Parent',figh1);
hold(axh1,'on');
box(axh1,'on');

colors = turbo(length(NoteData.F0));

classes = unique(pcaData.classLabel, 'stable');
nClasses = numel(classes)-1;
instruments = {'VocalAlto', ...
               'Violin', ...
               'Flute', ...
               'RecorderAlto', ...
               'ClarinetBb', ...
               'SaxophoneAlto', ...
               'OboeFrench', ...
               'Bassoon', ...
               'TrumpetC', ...
               'TromboneTenor', ...
               'Horn', ...
               'Tuba'};
instrNames = {'Alto voice',...
    'Violin',...
    'Flute',...
    'Alto recorder',...
    'Clarinet',...
    'Alto saxophone',...
    'Oboe',...
    'Bassoon',...
    'Trumpet',...
    'Trombone',...
    'Horn',...
    'Tuba'};
nInstr = numel(instruments);

t = tiledlayout(3, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
xlh = xlabel(t, ['PC1 (',num2str(round(pcaData.explained(1),2)),'\%)'],...
    'FontSize',pp.fsize+1,...
    'Color','none','Interpreter','latex');
% xlh.Position(2) = xlh.Position(2) - 0.1;
annotation('textbox',[0.369 0.025 0.25 0.05],...
    'String',['PC1 (',num2str(round(pcaData.explained(1),2)),'\%)'],...
    'EdgeColor','none','FontSize',pp.fsize+1,'Interpreter','latex')
ylabel(t, ['PC2 (',num2str(round(pcaData.explained(2),2)),'\%)'],...
    'FontSize', pp.fsize+1,'Interpreter','latex')
for iInstr = 1:nInstr
    instrName = instruments{iInstr};
    h(iInstr) = nexttile(t);
    hold on
    % plot black dots
    scatter(pcaData.score(:,1), pcaData.score(:,2), pp.scatterMsize-5, [0.7 0.7 0.7], 'Marker' ,'.');
    line([-9 8], [0 0], 'Color', 'k', 'LineWidth', pp.linewidth-0.7, 'LineStyle', '-');
    line([0 0], [-6 7], 'Color', 'k', 'LineWidth', pp.linewidth-0.7, 'LineStyle', '-');
    % plot instrument
    instrIdx = find(strcmp(pcaData.instrLabel, instrName));
    data = pcaData.score(instrIdx,:);
    colorIdx = rangeIdx.(instrName);
    % Low register
    lowRegInd = an.regis_ind_adj.(instrName).low;
    scatter(data(lowRegInd,1), data(lowRegInd,2), pp.scatterMsize-3,...
        colors(colorIdx(lowRegInd),:), 'filled',...
        'Marker', 'v', 'MarkerEdgeColor', 'k', 'LineWidth', pp.linewidth-0.6)
    pLow = plot(NaN,NaN,'vk','MarkerSize',pp.markersize-2);
    % Middle register
    midRegInd = an.regis_ind_adj.(instrName).mid;
    scatter(data(midRegInd,1), data(midRegInd,2), pp.scatterMsize-3,...
        colors(colorIdx(midRegInd),:), 'filled',...
        'Marker', 'o', 'MarkerEdgeColor', 'k', 'LineWidth', pp.linewidth-0.6)
    pMid = plot(NaN,NaN,'ok','MarkerSize',pp.markersize-2);
    % High register
    higRegInd = an.regis_ind_adj.(instrName).hig;
    scatter(data(higRegInd,1), data(higRegInd,2), pp.scatterMsize-3,...
        colors(colorIdx(higRegInd),:), 'filled',...
        'Marker', '^', 'MarkerEdgeColor', 'k', 'LineWidth', pp.linewidth-0.6)
    pHig = plot(NaN,NaN,'^k','MarkerSize',pp.markersize-2);

%     text(-8.5, 6.5, instrName, 'FontSize', pp.fsize-5.5, 'FontWeight', 'bold')
    hold off
    title(['\textbf{',instrNames{iInstr},'}'])
%     xlabel('1st principal component')
%     ylabel('2nd principal component')
    xlim([-9 8])
    ylim([-6 7])
%     if iInstr ~= 1 || 5 || 9
%         set(gca, 'YTickLabel', []);
%     else
%         set(gca, 'YTickLabel', [-5 0 5])
%     end
%     if iInstr ~= 9 || 10 || 11 || 12
%         set(gca, 'XTickLabel', []);
%     else
%         set(gca, 'XTickLabel', [-10 -5 0 5 10])
%     end
    set(gca, 'FontSize', pp.fsize-2, 'LineWidth', pp.linewidth,'Box','on')
    pbaspect([sum(abs(xlim))/sum(abs(ylim)), 1, 1])
%     xlm = xlim;
%     ylm = ylim;
%     title(instrName, 'FontSize', pp.fsize-2,...
%         'Position', [xlm(1) ylm(2)], 'HorizontalAlignment', 'left')
    if iInstr == 1
        lgd = legend([pLow,pMid,pHig],{'Low register', 'Middle register', 'High register'},...
            'Position',[0.596 0.022 0.32 0.045],...
            'FontSize', pp.fsize-4,...
            'LineWidth',pp.linewidth);
        lgd.NumColumns = 3;
        lgd.ItemTokenSize(1) = 2;
    end
end
set(h, 'Colormap', turbo)
cbNoteNames = {'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8'};
mp2cb = linspace(0,1,numel(NoteData.MP));
for iNoteName = 1:numel(cbNoteNames)
    noteIdx = find(strcmp(NoteData.MP, cbNoteNames{iNoteName}));
    cbIdx(iNoteName) = mp2cb(noteIdx);
end
cbh = colorbar(h(end), 'Ticks', cbIdx, 'TickLabels', cbNoteNames,...
    'LineWidth', pp.linewidth, 'FontSize', pp.fsize-2);
cbh.Layout.Tile = 'east';
cbh.Label.String = 'F0 (pitch)';
cbh.Label.FontSize = pp.fsize+1;
cbh.Label.Interpreter = 'latex';
set(cbh,'TickLabelInterpreter','latex')
% cbh.Label.FontWeight = 'bold';

% cbh.Position(3) = cbh.Position(3) - 0.01;


if isfield(pp, 'figwidth')
    if ~isempty(pp.figwidth)
        set(figh1, 'PaperPositionMode', 'manual');
        set(figh1, 'PaperUnits', 'centimeters');
        set(figh1, 'PaperPosition', [0 0 pp.figwidth pp.figheight]);
        set(figh1, 'Toolbar', 'none')
        set(figh1, 'Menubar', 'none')
        set(figh1, 'PaperSize', [pp.figwidth pp.figheight])
    end
end

print(figh1, [fig_folder filesep 'cluster.pdf'], '-dpdf');
