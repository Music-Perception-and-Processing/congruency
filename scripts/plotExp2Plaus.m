clear
close all

run plot_properties.m

set(0,'defaultAxesFontName',pp.fname)
set(0,'defaultTextInterpreter','tex')
set(0,'defaultAxesTickLabelInterpreter','tex')
set(0,'defaultLegendInterpreter','tex')

load('data/exp2_data.mat');
expData = expData_norm;

figh2 = figure('visible', pp.visible,...
    'DefaultTextFontName', pp.fname,...
    'DefaultAxesFontName', pp.fname);
axh2 = axes('Parent', figh2);
hold(axh2, 'on');
box(axh2, 'on');

regIdx.Violin = [2 3 4];
regIdx.VocalAlto = [2 3 4];
regIdx.ClarinetBb = [2 3 4];
regIdx.Tuba = [1 3 4];

tl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
c = lines(3);
c(4,:) = [0 0 0];
m = {'v','o','^','s'};
% t = {'A','(b)','(c)','(d)'};
% os1 = [0 0 0 0];
os1 = [-0.1 -0.033333 0.033333 0.1];
os2 = [0 0 0];

for iInstr = 1:nInstr
    instrName = instruments{iInstr};

    nexttile
    hold on

    %     % find ~nan idx of 'all'
    %     idx1 = find(~isnan(expData.quality.(instrName)(1,:,4)));
    %
    %     % plot instrument range region
    %     patch([idx1 fliplr(idx1)],[ones(1,numel(idx1))*6 ones(1,numel(idx1))],...
    %         'k',...
    %         'FaceAlpha',0.04,...
    %         'EdgeColor','none')
    %
    %      % plot register locations
    %     for ii = 1:3
    %         plot([r.(instrName)(ii) r.(instrName)(ii)],[1 6],...
    %             'Color',c(ii,:),...
    %             'LineWidth',pp.linewidth,...
    %             'LineStyle',':')
    %     end

    for iCond = 1:nCond
        condName = conditions{iCond};

        if iCond == 1
            %             patch([0.7 3.3 3.3 0.7],[6 6 1 1],...
            %                 'k',...
            %                 'FaceAlpha',0.04,...
            %                 'EdgeColor','none')

            %             patch([21.5 22.5 22.5 21.5],[6 6 1 1],...
            %                 'k',...
            %                 'FaceAlpha',0.04,...
            %                 'EdgeColor','none')

            pl = plot([0 5],[3.5 3.5],...
                'Color',[0.5 0.5 0.5],...
                'LineStyle',':',...
                'LineWidth',pp.linewidth-0.3);
            uistack(pl,'bottom')

        end

        data = expData.plausib.(instrName)(:,~isnan(expData.plausib.(instrName)(1,:,iCond)),iCond);

        if iCond < 4
            data = data(:,regIdx.(instrName));
        end

        nNote = numel(data(1,:));
        for iNote = 1:nNote
            CI(:,iNote) = bootci(1000,@mean,data(:,iNote));

        end

        instrMean = mean(data);

        meanCI = bootci(1000,@mean,reshape(data(:,:),numel(data),1));
        meanNegE = mean(instrMean - meanCI(2));
        meanPosE = meanCI(1) - mean(instrMean);

        negE = instrMean-CI(2,:);
        posE = CI(1,:)-instrMean;

        instrIdx = find(~isnan(expData.plausib.(instrName)(1,:,iCond)));

        p(iCond) = errorbar((1:3)+os1(iCond),instrMean,negE,posE,...
            'Color',c(iCond,:),...
            'Marker',m{iCond},...
            'MarkerSize',pp.markersize-1,...
            'LineStyle','-',...
            'LineWidth',pp.linewidth,...
            'CapSize',0);

        errorbar(4+os1(iCond),mean(instrMean),meanNegE,meanPosE,...
            'Color',c(iCond,:),...
            'Marker',m{iCond},...
            'MarkerSize',pp.markersize-1,...
            'LineWidth',pp.linewidth,...
            'CapSize',0)

        clear CI

    end

    hold off

    xlim([0.5 4.5])
    ylim([1 6])

    xticks([1 2 3 4])
    xticklabels({'low','mid','high','\mu'})
    xtickangle(0)
    yticks([1 3.5 6])

    if iInstr == 4
        legh = legend(p,{'lrSE','mrSE','hrSE','congr.'},...
            'FontSize',pp.fsize-2,...
            'NumColumns',2,...
            'Location','southwest');
            legh.ItemTokenSize = [15 15 0];
    end

    %     text(-0.4,5.75,instrName,'FontSize',pp.fsize-2,'FontWeight','bold')

    title(instrumentNames{iInstr})

    %     title(instrName)

    set(gca,'FontSize',pp.fsize-2,'LineWidth',pp.linewidth,'Layer','top','FontName',pp.fname)
    box on

end
% title(tl,['N = ',num2str(nPart)],'FontSize',24,'FontWeight','bold')
xlabel(tl,'Register','FontSize',pp.fsize+1,'Interpreter','tex','FontName',pp.fname)
ylabel(tl,'Plausibility rating','FontSize',pp.fsize+1,'Interpreter','tex','FontName',pp.fname)
set(gca,'FontSize',pp.fsize-2,'LineWidth',pp.linewidth,'Layer','top')

legh.Position(1) = legh.Position(1)-0.24;
legh.Position(2) = legh.Position(2)+0.08;

if isfield(pp, 'figwidth')
    if ~isempty(pp.figwidth)
        set(figh2, 'PaperPositionMode', 'manual');
        set(figh2, 'PaperUnits', 'centimeters');
        set(figh2, 'PaperPosition', [0 0 7 7]);
        set(figh2, 'Toolbar', 'none')
        set(figh2, 'Menubar', 'none')
        set(figh2, 'PaperSize', [7 7])
    end
end

if pp.print
    print(figh2, [fig_folder filesep 'plausibility.pdf'], '-dpdf');
end
