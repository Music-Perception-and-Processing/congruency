clear
close all

run plot_properties.m

set(0,'defaultAxesFontName',pp.fname)
set(0,'defaultTextInterpreter','tex')
set(0,'defaultAxesTickLabelInterpreter','tex')
set(0,'defaultLegendInterpreter','tex')

load('data/expData_1.mat')
load('data/F0data.mat')
load('data/pca_data_v5.mat')

figh4 = figure('Visible',pp.visible,...
    'DefaultTextFontName',pp.fname,...
    'DefaultAxesFontName',pp.fname);
axh4 = axes('Parent',figh4);
hold(axh4,'on');
box(axh4,'on');

% define variables
d = {'dim1','dim2'};
pcName = {'PC1','PC2'};
dP.dim1 = 1:9;
dP.dim2 = 10:16;
dC.dim1 = linspace(-20,20,9);
dC.dim2 = linspace(-15,15,7);
c = {'congruentF0','incongruentF0','fixedF0'};
cs = {'cF0','iF0','fF0'};
os.dim1 = [-0.4 0 0.4];
os.dim2 = [-0.3 0 0.3];

coob = {3:7,3:5};

% define tiled layout
t = tiledlayout(1,3,'TileSpacing','loose','Padding','compact');
ylabel(t,'Brightness rating','FontSize',pp.fsize,'Interpreter','tex','FontName',pp.fname)
color = lines(3);
marker = {'o','s','^'};
tileNum = [1,3];

for id = 1:numel(d)

    nexttile,hold on
    dN = d{id};

    line([min(dC.(dN))*1.1 max(dC.(dN))*1.1],[3.5 3.5],'Color',[0.5 0.5 0.5],'LineStyle',':','LineWidth',pp.linewidth)

    
    for ic = 1:numel(c)

        cN = c{ic};
        data = mean(expData.brightness.(cN).normStimResp(dP.(dN),:),2);
%         subjData = expData.brightness.(cN).normStimResp(dP.(dN),:);

        for idp = 1:length(data)
            
            allData = expData.brightness.(cN).normStimResp(dP.(dN)(idp),:);
%             SEM = std(allData)/sqrt(length(allData));         % Standard Error
%             ts = tinv([0.025  0.975],length(allData)-1);      % T-Score
%             CI(idp,:) = mean(allData) + ts*SEM;               % Confidence Intervals
            ci(idp,:) = bootci(1000,@mean,allData);



        end

        disp(cN)
        [rho,pval] = corr(dC.(dN)',data,'type','Spearman')

        negE = data - ci(:,1);
        posE = ci(:,2) - data;
        e(ic) = errorbar(dC.(dN)+os.(dN)(ic),data,negE,posE,...
            'Color',color(ic,:),...
            'Marker',marker(ic),...
            'MarkerSize',pp.markersize-1,...
            'LineWidth',pp.linewidth,...
            'CapSize',0);

%         plot(dC.(dN)(coob{id})+os.(dN)(ic),data(coob{id}),...
%             'Color',color(ic,:),...
%             'Marker',marker(ic),...
%             'MarkerFaceColor',color(ic,:),...
%             'LineStyle','none')

        if ismember(ic,[1 2])

            set([e(ic).Bar e(ic).Line],...
                'ColorType','truecoloralpha',...
                'ColorData',[e(ic).Line.ColorData(1:3);255*0.4])
            set([e(ic).Cap e(ic).MarkerHandle],...
                'EdgeColorType','truecoloralpha',...
                'EdgeColorData',[e(ic).Line.ColorData(1:3);255*0.4])

        end

    end

    hold off

    if id == 1
        title('A','Position',[-28 5.9 0],'HorizontalAlignment','center')
    else
        title('B','Position',[-21 5.9 0],'HorizontalAlignment','center')
    end
    if id == 2
        leg = legend(e,{'congruent','incongruent','fixed'},...
            'Location','south',...
            'Orientation','horizontal',...
            'NumColumns',3,...
            'FontSize',pp.fsize-2);
        leg.ItemTokenSize = [10 10 0];
    end
    xlim([-max(dC.(dN))*1.1, max(dC.(dN))*1.1])
    ylim([1 6])
    yticks([1 3.5 6])

    set(gca,'LineWidth',pp.linewidth,'FontSize',pp.fsize-2,'Layer','top','Box','on','FontName',pp.fname)
    xlabel([pcName{id},' (',num2str(round(pcaData.explained(id),2)),'%)'],...
        'FontSize',pp.fsize)
    clear ci

end

% next one
conditions = {'congruentF0', 'incongruentF0'};
altConditions = {'congruent', 'incongruent'};
nCond = numel(conditions);
space = {'inner','outer'};

nexttile, hold on
colors = lines(2);
loc = [1.65,4.2;2.6,5.5];

for iCond = 1:nCond

    condName = conditions{iCond};
    altCondName = altConditions{iCond};

    x = log10(F0data.brightness.(altCondName));
    X = [ones(length(x),1) x];
    y = expData.brightness.(condName).meanNormStimResp;

    b = X\y;
    yCalc = X*b;
    R2(iCond) = 1 - sum((y - yCalc).^2)/sum((y - mean(y)).^2);
    R = corrcoef(x,y);
    [~,~,~,~,stats.(condName)] = regress(y,X);
    stats.(condName)

    line([1.5 3.5],[3.5 3.5],...
        'Color',[0.5 0.5 0.5],...
        'LineWidth',pp.linewidth,...
        'LineStyle',':')

    sc(iCond) = scatter(x,y,...
        pp.scatterMsize-1.5,...
        colors(iCond,:),'filled',...
        'Marker',marker{iCond},...
        'MarkerEdgeAlpha',0.5,...
        'MarkerFaceAlpha',0.5);

    fit(iCond) = plot(x,yCalc,'-','LineWidth',pp.linewidth,...
        'Color',color(iCond,:));

    text(loc(iCond,1),loc(iCond,2),['R^2 = .',num2str(100*round(R2(iCond),2))],...
        'Color',colors(iCond,:),'Interpreter','tex','FontNAme',pp.fname)

end

leg1 = legend([sc fit],{'congr. data',...
    'incongr. data',...
    'linear fit',...
    'linear fit'},...
    'Location','southeast','NumColumns',1,'FontSize',pp.fsize-2);
leg1.ItemTokenSize = [8 8 0];


hold off

title('C','Position',[1.23 5.9 0],'HorizontalAlignment','center')


xlim([1.5 3.5])
xticks([1.5 2 2.5 3 3.5])

ylim([1 6])
yticks([1 3.5 6])

set(gca,'LineWidth',pp.linewidth,'FontSize',pp.fsize-2,'Layer','top','Box','on','FontName',pp.fname)
xlabel('Log frequency','FontWeight','normal','FontSize',pp.fsize)
% ylabel('Sound brightness rating','FontWeight','normal','FontSize',pp.fsize)

leg.Position(1) = leg.Position(1)-0.17;
leg.Position(2) = leg.Position(2)+0.155;
leg.Position(3) = leg.Position(3)+0.05;
leg.Position(4) = leg.Position(4)+0.06;
% leg.Position(4) = leg.Position(4)+0.02;

leg1.Position(1) = leg1.Position(1) + 0.018;
leg1.Position(2) = leg1.Position(2) + 0.155;
leg1.Position(3) = leg1.Position(3) + 0.03;
leg1.Position(4) = leg1.Position(4) + 0.15;

if isfield(pp, 'figwidth')
    if ~isempty(pp.figwidth)
        set(figh4, 'PaperPositionMode', 'manual');
        set(figh4, 'PaperUnits', 'centimeters');
        set(figh4, 'PaperPosition', [0 0 pp.figwidth 4.5]);
        set(figh4, 'Toolbar', 'none')
        set(figh4, 'Menubar', 'none')
        set(figh4, 'PaperSize', [pp.figwidth 4.5])
    end
end

print(figh4, [fig_folder filesep 'SBRcombined.pdf'], '-dpdf');
