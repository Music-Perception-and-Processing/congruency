clear
close all

load('data/expData_norm.mat');
expData = expData_norm;
clear expData_norm

load('data/old_expData_norm.mat');
expData_old = expData_norm;
clear expData_norm

load('data/levels.mat')
newLvl = 0.1531;

run plot_properties.m

instrShorts = {'vn','alt','cl','tb'};

%% plotting

figh6 = figure('visible', pp.visible,...
    'DefaultTextFontName', pp.fname,...
    'DefaultAxesFontName', pp.fname);
axh1 = axes('Parent', figh6);
hold(axh1, 'on');
box(axh1, 'on');

tl = tiledlayout(1,4,'TileSpacing','loose','Padding','compact');
c = lines(3);
c(4,:) = [0 0 0];
m = {'v','o','^','s'};

for iInstr = 1:nInstr
    instrName = instruments{iInstr};

    nexttile
    hold on

    for iCond = 1:nCond
        condName = conditions{iCond};
        p1(iCond) = plot(1:19,mean(expData.quality.(instrName)(:,:,iCond))-mean(expData_old.quality.(instrName)(:,:,iCond)),...
            'Color',c(iCond,:),...
            'LineWidth',pp.linewidth);

        plot(22,mean(mean(expData.quality.(instrName)(:,:,iCond))-mean(expData_old.quality.(instrName)(:,:,iCond)),'omitnan'),...
            'Color',c(iCond,:),...
            'Marker',m{iCond},...
            'MarkerSize',pp.markersize-1,...
            'LineWidth',pp.linewidth);

    end

    hold off

    if iInstr == 4
        legh = legend(p1,{'LRA','MRA','HRA','congr.'},...
            'FontSize',pp.fsize-2,...
            'NumColumns',4,...
            'Location','best');
        legh.ItemTokenSize = [12 1 0];
        legh.Position(1) = legh.Position(1)-0.4;
        legh.Position(2) = legh.Position(2)+0.5;

    end

    xlim([-1 24])
    xticks([1,4,7,10,13,16,19,22])
    xticklabels({'1','2','3','4','5','6','7','\mu'})
    xtickangle(0)

    ylim([-1 1])
    yticks([-1 0 1])

    title(instrumentNames{iInstr})

    set(gca,'FontSize',pp.fsize-2,'LineWidth',pp.linewidth,'Layer','top')
    box on

end

xlabel(tl,'Pitch re F#','FontSize',pp.fsize+1,'FontWeight','normal')
ylabel(tl,'Diff. response','FontSize',pp.fsize+1,'FontWeight','normal')

if isfield(pp, 'figwidth')
    if ~isempty(pp.figwidth)
        set(figh6, 'PaperPositionMode', 'manual');
        set(figh6, 'PaperUnits', 'centimeters');
        set(figh6, 'PaperPosition', [0 0 pp.figwidth 4.5]);
        set(figh6, 'Toolbar', 'none')
        set(figh6, 'Menubar', 'none')
        set(figh6, 'PaperSize', [pp.figwidth 4.5])
    end
end

print(figh6, [fig_folder filesep 'response_diff.pdf'], '-dpdf');

%%
figh7 = figure('visible', pp.visible,...
    'DefaultTextFontName', pp.fname,...
    'DefaultAxesFontName', pp.fname);
axh1 = axes('Parent', figh7);
hold(axh1, 'on');
box(axh1, 'on');

hold on
for iInstr = 1:nInstr
    instrName = instruments{iInstr};

    for iCond = 1:nCond
        condName = conditions{iCond};
        lvl(iCond,iInstr) = 20*log10(level(1,iCond,iInstr)/max(max(level)));
        lvl(iCond,iInstr) = 20*log10(level(1,iCond,iInstr)/newLvl);
        diff(iCond,iInstr) = mean(mean(expData.quality.(instrName)(:,:,iCond))-mean(expData_old.quality.(instrName)(:,:,iCond)),'omitnan');
        p(iCond) = scatter(lvl(iCond,iInstr),diff(iCond,iInstr),pp.scatterMsize+20,c(iCond,:),...
            'Marker',m{iCond},...
            'LineWidth',pp.linewidth);
        text(lvl(iCond,iInstr)-0.1,diff(iCond,iInstr)+0.02,instrShorts{iInstr},...
            'Color',c(iCond,:),...
            'FontSize',pp.fsize-2);

    end

end

x = reshape(lvl,numel(lvl),1);
X = [ones(size(x)) x];
y = reshape(diff,numel(diff),1);

[b,~,~,~,stats] = regress(y,X);
xCalc = -2:5;
yCalc = b(1) + b(2)*xCalc;

p2 = plot(xCalc,yCalc,...
    'Color',[0.5 0.5 0.5],...
    'LineWidth',pp.linewidth);
text(0,-0.36,['R^2 = ',num2str(round(stats(1),2))],...
    'Color',[0.5 0.5 0.5],...
    'FontSize',pp.fsize)

uistack(p2,'bottom')

hold off

legend(p,conditionNames,...
    'Location','southeast')

set(gca,'FontSize',pp.fsize-2,'LineWidth',pp.linewidth,'Layer','top')
box on

xlim([-2 5])
ylim([-0.6 0])

xlabel('Relative level / dB','FontSize',pp.fsize+1,'FontWeight','normal')
ylabel('Diff. response','FontSize',pp.fsize+1,'FontWeight','normal')

if isfield(pp, 'figwidth')
    if ~isempty(pp.figwidth)
        set(figh7, 'PaperPositionMode', 'manual');
        set(figh7, 'PaperUnits', 'centimeters');
        set(figh7, 'PaperPosition', [0 0 pp.figwidth pp.figheight]);
        set(figh7, 'Toolbar', 'none')
        set(figh7, 'Menubar', 'none')
        set(figh7, 'PaperSize', [pp.figwidth pp.figheight])
    end
end

print(figh7, [fig_folder filesep 'diff_vs_level.pdf'], '-dpdf');