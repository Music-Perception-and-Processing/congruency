clear
close all

load('/home/simon/ma/code/soundQuality/data/expData_norm.mat')
expData = expData_norm;
run plot_properties.m

figh2 = figure('visible', pp.visible,...
    'DefaultTextFontName', pp.fname,...
    'DefaultAxesFontName', pp.fname);
axh2 = axes('Parent', figh2);
hold(axh2, 'on');
box(axh2, 'off');

regIdx.Violin = [2 3 4];
regIdx.VocalAlto = [2 3 4];
regIdx.ClarinetBb = [2 3 4];
regIdx.Tuba = [1 3 4];

tl = tiledlayout(1,4,'TileSpacing','loose','Padding','compact');
c = lines(3);
c(4,:) = [0 0 0];
m = {'v','o','^','s'};
% t = {'A','(b)','(c)','(d)'};
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

        nNote = numel(data(1,:));
        for iNote = 1:nNote
            CI(:,iNote) = bootci(1000,@mean,data(:,iNote));

        end

        instrMean = mean(data);

        if iCond < 4
            meanCI = bootci(1000,@mean,reshape(data(:,regIdx.(instrName)),numel(data(:,regIdx.(instrName))),1));
            meanNegE = mean(instrMean(regIdx.(instrName))) - meanCI(2);
            meanPosE = meanCI(1) - mean(instrMean(regIdx.(instrName)));
        else
            meanCI = bootci(1000,@mean,reshape(data,numel(data),1));
            meanNegE = mean(instrMean) - meanCI(2);
            meanPosE = meanCI(1) - mean(instrMean);
        end

        negE = instrMean-CI(2,:);
        posE = CI(1,:)-instrMean;

        instrIdx = find(~isnan(expData.plausib.(instrName)(1,:,iCond)));

        if iInstr == 1
            if iCond < 4
                p(iCond) = errorbar([1 2 3]+os1(iCond),instrMean(2:end),negE(2:end),posE(2:end),...
                    'Color',c(iCond,:),...
                    'Marker',m{iCond},...
                    'MarkerSize',pp.markersize-1,...
                    'LineStyle','-',...
                    'LineWidth',pp.linewidth,...
                    'CapSize',0);
%                 errorbar(0+os2(iCond),instrMean(1),negE(1),posE(1),...
%                     'Color',c(iCond,:),...
%                     'Marker',m{iCond},...
%                     'MarkerSize',pp.markersize-1,...
%                     'LineWidth',pp.linewidth,...
%                     'CapSize',0);

                errorbar(4+os1(iCond),mean(instrMean(2:end)),meanNegE,meanPosE,...
                    'Color',c(iCond,:),...
                    'Marker',m{iCond},...
                    'MarkerSize',pp.markersize-1,...
                    'LineWidth',pp.linewidth,...
                    'CapSize',0)
            end
            clear CI
        end

        if iInstr == 2
            if iCond < 4
                p(iCond) = errorbar([1 2 3]+os1(iCond),instrMean(2:end-1),negE(2:end-1),posE(2:end-1),...
                    'Color',c(iCond,:),...
                    'Marker',m{iCond},...
                    'MarkerSize',pp.markersize-1,...
                    'LineStyle','-',...
                    'LineWidth',pp.linewidth,...
                    'CapSize',0);

%                 errorbar([0 4]+os2(iCond),instrMean([1 end]),negE([1 end]),posE([1 end]),...
%                     'Color',c(iCond,:),...
%                     'Marker',m{iCond},...
%                     'MarkerSize',pp.markersize-1,...
%                     'LineStyle','none',...
%                     'LineWidth',pp.linewidth,...
%                     'CapSize',2);

                errorbar(4+os1(iCond),mean(instrMean(2:end-1)),meanNegE,meanPosE,...
                    'Color',c(iCond,:),...
                    'Marker',m{iCond},...
                    'MarkerSize',pp.markersize-1,...
                    'LineWidth',pp.linewidth,...
                    'CapSize',0)
            end
            clear CI
        end

        if iInstr == 3
            if iCond < 4
                p(iCond) = errorbar([1 2 3]+os1(iCond),instrMean(2:end-1),negE(2:end-1),posE(2:end-1),...
                    'Color',c(iCond,:),...
                    'Marker',m{iCond},...
                    'MarkerSize',pp.markersize-1,...
                    'LineStyle','-',...
                    'LineWidth',pp.linewidth,...
                    'CapSize',0);
% 
%                 errorbar([0 4]+os2(iCond),instrMean([1 end]),negE([1 end]),posE([1 end]),...
%                     'Color',c(iCond,:),...
%                     'Marker',m{iCond},...
%                     'MarkerSize',pp.markersize-1,...
%                     'LineStyle','none',...
%                     'LineWidth',pp.linewidth,...
%                     'CapSize',2);

                errorbar(4+os1(iCond),mean(instrMean(2:end-1)),meanNegE,meanPosE,...
                    'Color',c(iCond,:),...
                    'Marker',m{iCond},...
                    'MarkerSize',pp.markersize-1,...
                    'LineWidth',pp.linewidth,...
                    'CapSize',0)
            end
            clear CI
        end

        if iInstr == 4
            if iCond < 4
                p(iCond) = errorbar([1 2 3]+os1(iCond),instrMean([1 3 4]),negE([1 3 4]),posE([1 3 4]),...
                    'Color',c(iCond,:),...
                    'Marker',m{iCond},...
                    'MarkerSize',pp.markersize-1,...
                    'LineStyle','-',...
                    'LineWidth',pp.linewidth,...
                    'CapSize',0);

                %                 errorbar(4+os2(iCond),instrMean(2),negE(2),posE(2),...
                %                     'Color',c(iCond,:),...
                %                     'Marker',m{iCond},...
                %                     'MarkerSize',pp.markersize-2,...
                %                     'LineWidth',pp.linewidth-0.5);

                errorbar(4+os1(iCond),mean(instrMean([1 3 4])),meanNegE,meanPosE,...
                    'Color',c(iCond,:),...
                    'Marker',m{iCond},...
                    'MarkerSize',pp.markersize-1,...
                    'LineWidth',pp.linewidth,...
                    'CapSize',0)
            end

            clear CI
        end

        % cond all
        if iCond == 4
            p(iCond) = errorbar([1 2 3]+os1(4),instrMean,negE,posE,...
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
        end

        clear CI

    end

    hold off

    xlim([0.5 4.5])
    ylim([1 6])

    xticks([1 2 3 4])
    xticklabels({'low','mid','high','\mu'})
    xtickangle(0)
    yticks([1 6])

    if iInstr == 4
        legh = legend(p,{'LRA','MRA','HRA','congr.'},...
            'FontSize',pp.fsize-2,...
            'NumColumns',2,...
            'Location','southwest');
            legh.ItemTokenSize = [15 15 0];
    end

    %     text(-0.4,5.75,instrName,'FontSize',pp.fsize-2,'FontWeight','bold')

    title(instrumentNames{iInstr})

    %     title(instrName)

    set(gca,'FontSize',pp.fsize-2,'LineWidth',pp.linewidth,'Layer','top')

end
% title(tl,['N = ',num2str(nPart)],'FontSize',24,'FontWeight','bold')
xlabel(tl,'Register','FontSize',pp.fsize+1,'FontWeight','normal')
ylabel(tl,'Plausibility rating','FontSize',pp.fsize+1,'FontWeight','normal')
set(gca,'FontSize',pp.fsize-2,'LineWidth',pp.linewidth,'Layer','top')

legh.Position(1) = legh.Position(1)-0.16;
legh.Position(2) = legh.Position(2)+0.17;

if isfield(pp, 'figwidth')
    if ~isempty(pp.figwidth)
        set(figh2, 'PaperPositionMode', 'manual');
        set(figh2, 'PaperUnits', 'centimeters');
        set(figh2, 'PaperPosition', [0 0 pp.figwidth 4.5]);
        set(figh2, 'Toolbar', 'none')
        set(figh2, 'Menubar', 'none')
        set(figh2, 'PaperSize', [pp.figwidth 4.5])
    end
end

print(figh2, [fig_folder filesep 'plausibility.pdf'], '-dpdf');


