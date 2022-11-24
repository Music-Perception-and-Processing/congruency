clear
close all

load('/home/simon/ma/code/soundQuality/data/expData_norm.mat')
expData = expData_norm;
run plot_properties.m

ee.Violin = 19;
ee.VocalAlto = 13;
ee.ClarinetBb = 16;
ee.Tuba = 10;

for iInstr = 1:nInstr
    instrName = instruments{iInstr};
    meanCDiff.(instrName) = nan(3,ee.(instrName));
end

meanDiff = [-0.2188,...
    -0.2314,...
    -0.3963,...
    -0.2462,...
    -0.1364,...
    -0.2348,...
    0.0856,...
    -0.1932,...
    -0.1487,...
    0.0000,...
    0.0455,...
    -0.0455];

figh4 = figure('visible', pp.visible,...
    'DefaultTextFontName', pp.fname,...
    'DefaultAxesFontName', pp.fname);
axh4 = axes('Parent', figh4);
hold(axh4, 'on');
box(axh4, 'off');

tl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
c = lines(3);
c(4,:) = [0 0 0];
m = {'v','o','^','s'};
t = {'(a)','(b)','(c)','(d)'};

diff_matrix = zeros(3,19,4);
for iInstr = 1:nInstr
    instrName = instruments{iInstr};

    nexttile
    hold on

    plot([0 12],[0 0],...
        'Color',[.7 .7 .7],...
        'LineStyle','-',...
        'LineWidth',pp.linewidth-0.3)

    pma = plot(0:11,meanDiff,...
        'Color',[0.3 0.3 0.3],...
        'LineStyle',':',...
        'LineWidth',pp.linewidth);

    mean_all = mean(expData.quality.(instrName)(:,:,4));

    for iCond = 1:nCond-1
        condName = conditions{iCond};

        mean_cond = mean(expData.quality.(instrName)(:,:,iCond));
        diff_matrix(iCond,:,iInstr) = mean_cond - mean_all;
        meanCDiff.(instrName)(iCond,1:ee.(instrName)-r.(instrName)(iCond)+1) = diff_matrix(iCond,r.(instrName)(iCond):ee.(instrName),iInstr);

        p(iCond) = plot(0:ee.(instrName)-r.(instrName)(iCond),diff_matrix(iCond,r.(instrName)(iCond):ee.(instrName),iInstr),...
            'Color',c(iCond,:),...
            'LineWidth',pp.linewidth,...
            'Marker',m{iCond},...
            'MarkerSize',pp.markersize-1.5);

    end

    pmc = plot(0:ee.(instrName)-1,mean(meanCDiff.(instrName),'omitnan'),...
        'Color','k',...
        'LineStyle','--',...
        'LineWidth',pp.linewidth);
    uistack(pmc,'down',3)

    legh = legend([p(1),p(2),p(3),pmc,pma],{'low','mid','high','instr. mean','overall mean'},...
        'Location','southeast',...
        'NumColumns',2);
    legh.ItemTokenSize = [15 0 0];

%     text(0.2,0.85,instrName,'FontSize',pp.fsize-2,'FontWeight','bold')

    xlim([0 12])
    xticks([0 3 6 9 12])
    xticklabels({'0','1','2','3','4'})
    ylim([-2 1])
    yticks([-2 -1 0 1])
    %     title(instrName)
    set(gca,'FontSize',pp.fsize-2,'LineWidth',pp.linewidth,'Layer','top')
    title(instrumentNames{iInstr})

end

% title(tl,['N = ',num2str(nPart)],'FontSize',24,'FontWeight','bold')
xlabel(tl,'Incongruency level / octaves','FontSize',pp.fsize+1,'FontWeight','normal')
ylabel(tl,'Pleasantness rating difference','FontSize',pp.fsize+1,'FontWeight','normal')
set(gca,'FontSize',pp.fsize-2,'LineWidth',pp.linewidth,'Layer','top')

if isfield(pp, 'figwidth')
    if ~isempty(pp.figwidth)
        set(figh4, 'PaperPositionMode', 'manual');
        set(figh4, 'PaperUnits', 'centimeters');
        set(figh4, 'PaperPosition', [0 0 pp.figwidth pp.figheight]);
        set(figh4, 'Toolbar', 'none')
        set(figh4, 'Menubar', 'none')
        set(figh4, 'PaperSize', [pp.figwidth pp.figheight])
    end
end

print(figh4, [fig_folder filesep 'iclvl.pdf'], '-dpdf');