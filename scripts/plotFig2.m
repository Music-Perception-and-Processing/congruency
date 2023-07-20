clear
close all

run plot_properties.m

set(0,'defaultAxesFontName',pp.fname)
set(0,'defaultTextInterpreter','tex')
set(0,'defaultAxesTickLabelInterpreter','tex')
set(0,'defaultLegendInterpreter','tex')

load('data/pca_data_v5.mat')
load('data/F0data.mat')
load('data/gridPoints.mat')
load('data/NoteData.mat')

addpath('functions/')

figh2 = figure('Visible',pp.visible,...
    'DefaultTextFontName',pp.fname,...
    'DefaultAxesFontName',pp.fname);
axh2 = axes('Parent',figh2);
hold(axh2,'on');
box(axh2,'on');

% plot settings
tl = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
xlabel(tl,'PC1','FontSize',pp.fsize+1,'Interpreter','tex','FontName',pp.fname)
ylabel(tl,'PC2','FontSize',pp.fsize+1,'Interpreter','tex','FontName',pp.fname)
c = turbo(numel(NoteData.F0));

% variables
F0 = F0data.quality.congruent;
notes = freq2muspitch(F0);
dim1 = gridPoints.A(:,1);
dim2 = gridPoints.A(:,2);
% plot space
for iValue = 1:numel(F0)

    cIdx(iValue) = find(strcmp(NoteData.MP,notes(iValue)));

end

h = nexttile;
scatter(dim1,dim2,pp.scatterMsize+30,c(cIdx,:),'filled',...
    'MarkerEdgeColor','k',...
    'Marker','square',...
    'LineWidth',pp.linewidth-0.8)

xlim([-25 25])
xticks([-20 -10 0 10 20])
% xlabel('PC1','FontWeight','bold')
ylim([-20 20])
yticks([-15 -7.5 0 7.5 15])
% ylabel('PC2','FontWeight','bold')


% pbaspect([sum(abs(xlim))/sum(abs(ylim)), 1, 1])

set(h, 'Colormap', c)
% cbNoteNames = {'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8'};
% mp2cb = linspace(0,1,numel(NoteData.MP));
% for iNoteName = 1:numel(cbNoteNames)
%     noteIdx = find(strcmp(NoteData.MP, cbNoteNames{iNoteName}));
%     cbIdx(iNoteName) = mp2cb(noteIdx);
% end
% cbh = colorbar(h(end),...
%     'Ticks', cbIdx,...
%     'TickLabels', cbNoteNames,...
%     'LineWidth',pp.linewidth,...
%     'FontSize', pp.fsize-4);

% cbh.Layout.Tile = 'east';
% cbh.Label.String = 'F0';
% cbh.Label.FontSize = pp.fsize;
% cbh.Label.HorizontalAlignment = 'center';

set(gca,...
    'LineWidth',pp.linewidth,...
    'FontSize',pp.fsize-2,...
    'Layer','top',...
    'Box','on','FontName',pp.fname)

title('A','Position',[-30 17 0],'FontSize',pp.fsize)


% incongruent
clear cIdx

% variables
F0 = F0data.quality.incongruent;
notes = freq2muspitch(F0);
dim1 = gridPoints.A(:,1);
dim2 = gridPoints.A(:,2);
% plot space
for iValue = 1:numel(F0)

    cIdx(iValue) = find(strcmp(NoteData.MP,notes(iValue)));

end

h = nexttile;
scatter(dim1,dim2,pp.scatterMsize+30,c(cIdx,:),'filled',...
    'MarkerEdgeColor','k',...
    'Marker','square',...
    'LineWidth',pp.linewidth-0.8)

xlim([-25 25])
xticks([-20 -10 0 10 20])
% xlabel('PC1','FontWeight','bold')
ylim([-20 20])
yticks([-15 -7.5 0 7.5 15])
% ylabel('PC2','FontWeight','bold')


% pbaspect([sum(abs(xlim))/sum(abs(ylim)), 1, 1])

set(h, 'Colormap', c)
cbNoteNames = {'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8'};
mp2cb = linspace(0,1,numel(NoteData.MP));
for iNoteName = 1:numel(cbNoteNames)
    noteIdx = find(strcmp(NoteData.MP, cbNoteNames{iNoteName}));
    cbIdx(iNoteName) = mp2cb(noteIdx);
end
cbh = colorbar(h(end),...
    'Ticks', cbIdx,...
    'TickLabels', cbNoteNames,...
    'LineWidth',pp.linewidth,...
    'FontSize', pp.fsize-4);

cbh.Layout.Tile = 'east';
cbh.Label.String = 'F0 (pitch)';
cbh.Label.FontSize = pp.fsize;
cbh.Label.HorizontalAlignment = 'center';
cbh.Label.Interpreter = 'tex';
set(cbh,'TickLabelInterpreter','tex')

set(gca,...
    'LineWidth',pp.linewidth,...
    'FontSize',pp.fsize-2,...
    'Layer','top',...
    'Box','on','FontName',pp.fname)
title('B','Position',[-30 17 0],'FontSize',pp.fsize)


if isfield(pp, 'figwidth')
    if ~isempty(pp.figwidth)
        set(figh2, 'PaperPositionMode', 'manual');
        set(figh2, 'PaperUnits', 'centimeters');
        set(figh2, 'PaperPosition', [0 0 pp.figwidth/2 8]);
        set(figh2, 'Toolbar', 'none')
        set(figh2, 'Menubar', 'none')
        set(figh2, 'PaperSize', [pp.figwidth/2 8])
    end
end

if pp.print
    print(figh2, [fig_folder filesep 'componentSpaceGrid.pdf'], '-dpdf');
end
