clear
close all

set(0,'defaultTextInterpreter','latex')
set(0,'defaultAxesTickLabelInterpreter','latex')
set(0,'defaultLegendInterpreter','latex')

load('data/expData_norm.mat');
expData = expData_norm;
run plot_properties.m

figh5 = figure('visible', pp.visible,...
    'DefaultTextFontName', pp.fname,...
    'DefaultAxesFontName', pp.fname);
axh5 = axes('Parent', figh5);
hold(axh5, 'on');
box(axh5, 'on');

% concatenate data
data = [];
c = [];
i = [];
varCond = [2 3 4 1];
hold on
for iInstr = 1:nInstr
    instrName = instruments{iInstr};
    for iCond = 1:nCond
        condName = conditions{iCond};
        data = [data;mean(expData.quality.(instrName)(:,:,iCond),'omitnan')];
%         s = scatter(1:19,mean(expData.quality.(instrName)(:,:,iCond),'omitnan'),pp.scatterMsize-6,[.3 .3 .3],'filled');
        c = [c;iCond*ones(1,19)];
        i = [i;iInstr*ones(1,19)];
    end
end

% % calculate CI
% nNote = 19;
% for iNote = 1:nNote
%     CI(:,iNote) = bootci(1000,{@(v)mean(v,'omitnan'),data(:,iNote)});
% end
% 
% confInt = patch([1:19 fliplr(1:19)],[CI(1,:) fliplr(CI(2,:))],...
%     'k',...
%     'FaceAlpha',0.1,...
%     'EdgeColor','none');

% % regression
% x1 = reshape(repmat(1:19,16,1),numel(data),1);
% x2 = reshape(c,numel(c),1);
% x3 = reshape(i,numel(i),1);
% y = reshape(data,numel(data),1);
% 
% % remove nans
% nanidx = ~isnan(y);
% x1 = x1(nanidx);
% x2 = x2(nanidx);
% x3 = x3(nanidx);
% y = y(nanidx);
% 
% tbl = table(x1,categorical(x2),categorical(x3),y,'VariableNames',{'pitch','cond','instr','resp'});
% lm = fitlm(tbl,'resp ~ 1 + pitch*pitch','CategoricalVars',[2 3])

nPart = 20;
partNum = [];
colors = lines(nPart);
for iPart = 1:nPart
    tempData = [];
    for iInstr = 1:nInstr
        instrName = instruments{iInstr};
        for iCond = 1:nCond
            tempData = [tempData;expData.quality.(instrName)(iPart,:,iCond)];
        end
    end
    partData(iPart,:) = mean(tempData,'omitnan');
    partNum = [partNum;iPart*ones(19,1)];

%     s = scatter(1:19,partData(iPart,:),pp.scatterMsize-6,[.3 .3 .3],'filled');


    s(iPart) = plot(1:19,partData(iPart,:),...
        'Color',[colors(iPart,:),0.3],...
        'LineWidth',pp.linewidth-0.99);
end

nNote = 19;
for iNote = 1:nNote
    CI(:,iNote) = bootci(1000,@mean,partData(:,iNote));
end

confInt = patch([1:19 fliplr(1:19)],[CI(1,:) fliplr(CI(2,:))],...
    [.5 .5 .5],...
    'FaceAlpha',0.5,...
    'EdgeColor','none');

%% lme
x = reshape(repmat(1:19,20,1),numel(partData),1);
y = reshape(partData,numel(partData),1);

% remove nans
nanidx = ~isnan(y);
x = x(nanidx);
y = y(nanidx);

tbl = table(x,categorical(partNum),y,'VariableNames',{'pitch','part','resp'});
lme = fitlme(tbl,'resp ~ 1 + pitch*pitch + (1|part)',...
        'FitMethod','REML','CheckHessian',1,'DummyVarCoding','effects')

%% participant clustering
Y = pdist(reshape(lme.residuals>0,nPart,nNote));
Z = linkage(Y);
T = cluster(Z,'maxclust',3);

%% plotting
fit = plot(x,lme.Fitted,...
    'Color','red',...
    'LineWidth',pp.linewidth+2);

p = plot(1:19,mean(partData),'-k','LineWidth',pp.linewidth);

% text(1,5.5,['R^2 = ',num2str(round(lme.Rsquared.Adjusted,2))],...
%     'Color','red')

hold off

legend([s(1),p,confInt,fit],{'participant data','mean','$95\%$ CI',['quadratic fit ($R^2 =$ .',num2str(100*round(lme.Rsquared.Adjusted,2)),')']},...
    'Location','best')

xlim([0 20])
xticks([1,4,7,10,13,16,19])
xticklabels({'F\#1','F\#2','F\#3','F\#4','F\#5','F\#6','F\#7'})
xtickangle(0)

ylim([1 6])
yticks([1 2 3 4 5 6])

set(gca,'FontSize',pp.fsize-2,'LineWidth',pp.linewidth,'Layer','top')

xlabel('Pitch','FontSize',pp.fsize)
ylabel('Sound pleasantness rating','FontSize',pp.fsize)

if isfield(pp, 'figwidth')
    if ~isempty(pp.figwidth)
        set(figh5, 'PaperPositionMode', 'manual');
        set(figh5, 'PaperUnits', 'centimeters');
        set(figh5, 'PaperPosition', [0 0 pp.figwidth pp.figheight]);
        set(figh5, 'Toolbar', 'none')
        set(figh5, 'Menubar', 'none')
        set(figh5, 'PaperSize', [pp.figwidth pp.figheight])
    end
end

print(figh5, [fig_folder filesep 'pleasantness_distribution.pdf'], '-dpdf');