clear
close all

run plot_properties.m

set(0,'defaultAxesFontName',pp.fname)
set(0,'defaultTextInterpreter','tex')
set(0,'defaultAxesTickLabelInterpreter','tex')
set(0,'defaultLegendInterpreter','tex')

load('data/expData_norm.mat')
expData = expData_norm;

v1.Violin = [3 8 11 16];
v1.VocalAlto = [3 8 10 13 18];
v1.ClarinetBb = [3 7 10 13 18];
v1.Tuba = [2 5 8];

figh3 = figure('visible', pp.visible,...
    'DefaultTextFontName', pp.fname,...
    'DefaultAxesFontName', pp.fname);
axh3 = axes('Parent', figh3);
hold(axh3, 'on');
box(axh3, 'on');

tl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
c = lines(3);
c(4,:) = [0 0 0];
m = {'v','o','^','s'};
t = {'(a)','(b)','(c)','(d)'};
legPos = [0.05,0.05,0.025,0.025];

for iInstr = 1:nInstr
    instrName = instruments{iInstr}

    nexttile
    hold on

    plot([1 6],[1 6],...
        'LineStyle',':',...
        'Color',[0.5 0.5 0.5],...
        'LineWidth',pp.linewidth)

    x_all = [];
    y_all = [];
    X_all = [];

    for iCond = 1:nCond-1
        condName = conditions{iCond};

        x = mean(expData.quality.(instrName)(:,v1.(instrName),iCond))';
        y = mean(expData.plausib.(instrName)(:,~isnan(expData.plausib.(instrName)(1,:,iCond)),iCond))';

        if iInstr == 4 && iCond < 4
            y = [y(1); y(3); y(4)];
        end

        x_all = [x_all;x];
        y_all = [y_all;y];

        X = [ones(size(x)) x];

        X_all = [X_all;X];

        [b,~,~,~,stats(iCond,:)] = regress(y,X);

        yhat = b(1) + b(2)*x;

        s(iCond) = scatter(x,y,pp.scatterMsize,'blue',...
            'Marker',m{iCond},...
            'LineWidth',pp.linewidth,...
            'MarkerEdgeAlpha',0.75);
        %         p(iCond) = plot(x,yhat,'Color',c(iCond,:),'LineWidth',pp.linewidth);

    end

    x = mean(expData.quality.(instrName)(:,r.(instrName),4))';
    y = mean(expData.plausib.(instrName)(:,~isnan(expData.plausib.(instrName)(1,:,4)),4))';

    x_all = [x_all;x];
    y_all = [y_all;y];

    X = [ones(size(x)) x];

    X_all = [X_all;X];

    [b,~,~,~,stats(4,:)] = regress(y,X);

    yhat = b(1) + b(2)*x;


    %     text(4.3,3,['R2 = ',num2str(round(stats(1,1),2)),', p = ',num2str(round(stats(1,3),3))],...
    %         'FontSize',pp.fsize-4,...
    %         'FontWeight','bold',...
    %         'Color',c(1,:))
    %     text(4.3,2.5,['R2 = ',num2str(round(stats(2,1),2)),', p = ',num2str(round(stats(2,3),3))],...
    %         'FontSize',pp.fsize-4,...
    %         'FontWeight','bold',...
    %         'Color',c(2,:))
    %     text(4.3,2,['R2 = ',num2str(round(stats(3,1),2)),', p = ',num2str(round(stats(3,3),3))],...
    %         'FontSize',pp.fsize-4,...
    %         'FontWeight','bold',...
    %         'Color',c(3,:))
    %     text(4.3,1.5,['R2 = ',num2str(round(stats(4,1),2)),', p = ',num2str(round(stats(4,3),3))],...
    %         'FontSize',pp.fsize-4,...
    %         'FontWeight','bold',...
    %         'Color','k')

    s(4) = scatter(x,y,pp.scatterMsize,'blue',...
        'Marker','s',...
        'LineWidth',pp.linewidth,...
        'MarkerEdgeAlpha',0.75);
    %     p(4) = plot(x,yhat,'Color','k','LineWidth',pp.linewidth);


    [b_all,~,~,~,stats_all] = regress(y_all,X_all);

    b_all
    stats_all

    yhat_all = b_all(1) + b_all(2)*x_all;

    p = plot(x_all,yhat_all,'Color','red','LineWidth',pp.linewidth);
    text(3.4,1.7,['R^2 = .',num2str(round(stats_all(1),2)*100)],...
        'Color','red',...
        'FontSize',pp.fsize-2,'Interpreter','tex','FontName',pp.fname,'FontWeight','normal')

    hold off

    xlim([1 6])
    xticks([1 2 3 4 5 6])
    ylim([1 6])
    yticks([1 2 3 4 5 6])

%     if iInstr == 1
%         legh = legend([s p],{'low data','mid data','high data','all data','linear fit'},...
%             'FontSize',pp.fsize-2,...
%             'Location','northwest',...
%             'NumColumns',5);
%         legh.ItemTokenSize = [5 0 0];
%     end

    %     text(1.1,5.7,instrName,'FontSize',pp.fsize-2,'FontWeight','bold')

    set(gca,'FontSize',pp.fsize-2,'LineWidth',pp.linewidth,'Layer','top','FontName',pp.fname)
    title(instrumentNames{iInstr})

    pbaspect([1 1 1])

    box on

end

xlabel(tl,'Pleasantness rating','FontSize',pp.fsize+1,'Interpreter','tex','FontName',pp.fname)
ylabel(tl,'Plausibility rating','FontSize',pp.fsize+1,'Interpreter','tex','FontName',pp.fname)

% legh.Position(1) = legh.Position(1)+0.12;

if isfield(pp, 'figwidth')
    if ~isempty(pp.figwidth)
        set(figh3, 'PaperPositionMode', 'manual');
        set(figh3, 'PaperUnits', 'centimeters');
        set(figh3, 'PaperPosition', [0 0 7 7]);
        set(figh3, 'Toolbar', 'none')
        set(figh3, 'Menubar', 'none')
        set(figh3, 'PaperSize', [7 7])
    end
end

if pp.print
    print(figh3, [fig_folder filesep 'regression.pdf'], '-dpdf');
end
