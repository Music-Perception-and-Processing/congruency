clear
close all

set(0,'defaultTextInterpreter','latex')
set(0,'defaultAxesTickLabelInterpreter','latex')
set(0,'defaultLegendInterpreter','latex')

load('data/expData_1.mat')
load('data/F0data.mat')
load('data/NoteData.mat')
run plot_properties.m

figh3 = figure('Visible',pp.visible,...
    'DefaultTextFontName',pp.fname,...
    'DefaultAxesFontName',pp.fname);
axh3 = axes('Parent',figh3);
hold(axh3,'on');
box(axh3,'on');

% plot settings
tl = tiledlayout(1,5,'TileSpacing','compact','Padding','compact');

nexttile(1,[1 3]), hold on

% define variables
conds = {'congruentF0','incongruentF0','fixedF0'};
cShort = {'cF0','iF0','fF0'};

space = {'inner','outer'};
sIdx.inner = 1:45;
sIdx.outer = 46:81;
sOS = [-0.15 0.15];
sOSName = [-0.15 0.2];
sAlpha = [0.2 0.2];
sNameShort = {'IN','OUT'};

colors = lines(3);
cW = [1,0.75];
markers = {'o','s','^'};

line([0 4],[3.5 3.5],'Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',pp.linewidth)

for iCond = 1:numel(conds)

    cName = conds{iCond};

    for iSpace = 1:numel(space)

        sName = space{iSpace};

        data.(cName).(sName) = mean(expData.quality.(cName).normStimResp(sIdx.(sName),:),2);
        subjData.(cName).(sName) = mean(expData.quality.(cName).normStimResp(sIdx.(sName),:))';


%         SEM = std(subjData.(cName).(sName))/sqrt(length(subjData.(cName).(sName)));
%         ts = tinv([0.025 0.975],length(subjData.(cName).(sName))-1);
%         err = ts*SEM;
%         CI = mean(subjData.(cName).(sName)) + err;

        CI = bootci(1000,@mean,subjData.(cName).(sName));
        negE = mean(subjData.(cName).(sName)) - CI(1);
        posE = CI(2) - mean(subjData.(cName).(sName));

%         y = sort(subjData.(cName).(sName));
%         loQ = y(floor(0.25*length(y))) - median(y);
%         hiQ = y(floor(0.75*length(y))) - median(y);

        sL = length(subjData.(cName).(sName));
        scatter(iCond+sOS(iSpace)+randn(sL,1).*0.01,subjData.(cName).(sName),...
            pp.scatterMsize,...
            colors(iCond,:).*cW(iSpace),'filled',...
            'MarkerEdgeAlpha',sAlpha(iSpace),...
            'MarkerFaceAlpha',sAlpha(iSpace),...
            'Marker',markers{iCond})

        e = errorbar(iCond+sOS(iSpace),mean(subjData.(cName).(sName)),negE,posE,...
            'Color',colors(iCond,:).*cW(iSpace),...
            'Marker',markers{iCond},...
            'LineWidth',pp.linewidth,...
            'MarkerSize',pp.markersize+2);

        mean(subjData.(cName).(sName))
        std(subjData.(cName).(sName))

%         if iSpace == 2
% 
%             set([e.Bar, e.Line],...
%                 'ColorType', 'truecoloralpha',...
%                 'ColorData', [e.Line.ColorData(1:3); 255*0.5])
%             set(e.Cap, 'EdgeColorType', 'truecoloralpha',...
%                 'EdgeColorData', [e.Cap.EdgeColorData(1:3); 255*0.5])
%             eMarkers = e.MarkerHandle;
%             set(eMarkers, 'EdgeColorData', [e.Cap.EdgeColorData(1:3); 255*0.5])
% 
% 
%         end

        text(iCond+sOSName(iSpace),1.2,['\textbf{',sNameShort{iSpace},'}',],...
            'HorizontalAlignment','center',...
            'FontWeight','bold',...
            'FontSize',pp.fsize-3)

    end

end

hold off

xlim([0.5 3.5])
xticks([1 2 3])
xticklabels({'congruent','incongruent','fixed'})
xtickangle(0)


ylim([1 6])
% yticks([1 2.25 4.75 6])
% yticklabels({'','unpleasant','pleasant',''})
yticks([1 2 3 4 5 6])

% ytickangle(90)

title('\textbf{A}','Position',[0.25 5.8 0],'HorizontalAlignment','center')

set(gca,'FontSize',pp.fsize-2,'LineWidth',pp.linewidth,'Layer','top','Box','on')
xlabel('Conditions','FontWeight','normal','FontSize',pp.fsize)
ylabel('Sound pleasantness rating','FontWeight','normal','FontSize',pp.fsize)

nexttile(4,[1 2]), hold on
% variables
nStim = 45;

cF0 = freq2muspitch(F0data.quality.congruent(1:nStim));
iF0 = freq2muspitch(F0data.quality.incongruent(1:nStim));
fF0 = freq2muspitch(F0data.quality.fixed(1:nStim));

for iStim = 1:nStim

    cIdx(iStim) = find(strcmp(NoteData.MP,cF0{iStim}));
    iIdx(iStim) = find(strcmp(NoteData.MP,iF0{iStim}));
    fIdx(iStim) = find(strcmp(NoteData.MP,fF0{iStim}));

end

iICLVL = abs((cIdx - iIdx)/2);
fICLVL = abs(cIdx - fIdx);

% regression
xI = iICLVL';
xF = fICLVL';
x = xI;
X = [ones(size(x)) x];
yC = mean(expData.quality.congruentF0.normStimResp(1:nStim,:),2);
yI = mean(expData.quality.incongruentF0.normStimResp(1:nStim,:),2);
yF = mean(expData.quality.fixedF0.normStimResp(1:nStim,:),2);
y = yI;
y1 = yC;

[b,~,~,~,stats] = regress(y1,X);
stats
[rho,pval] = corr(x,y1,'type','Spearman')
yCalc = X*b;

[b,~,~,~,stats] = regress(y,X);
stats
[rho,pval] = corr(x,y,'type','Spearman')
yCalc = X*b;

% quadratic fit
x0 = [1 1 1];
fitfun = fittype( 'squareFit( x, a, b, c )' );
[fitted_curve,gof] = fit(x,y,fitfun,'StartPoint',x0);
curve = fitted_curve(min(x):0.01:max(x));
coeffvals = coeffvalues(fitted_curve);

c = lines(2);

line([0 35],[3.5 3.5], ...
    'Color',[0.5 0.5 0.5], ...
    'LineStyle','--', ...
    'LineWidth',pp.linewidth)
s = scatter(x,y,pp.scatterMsize,c(2,:),'filled',...
    'Marker','square',...
            'MarkerEdgeAlpha',0.5,...
            'MarkerFaceAlpha',0.5);
% scatter(x,y1,pp.scatterMsize,c(1,:),'filled',...
%     'Marker','o',...
%     'MarkerEdgeAlpha',0.5,...
%     'MarkerFaceAlpha',0.5)
p = plot(x,yCalc,...
    'LineWidth',pp.linewidth, ...
    'Color',c(2,:));
% plot(min(x):0.01:max(x),curve)

text(16,2.6,['{\boldmath$R^2 =$}','\textbf{ .',num2str(100*round(stats(1),2)),'}'], ...
    'FontSize',pp.fsize,'Color',c(2,:))

hold off

ylim([1 6])
yticks([1 2 3 4 5 6])

xlim([0 35])
xticks([0 6 12 18 24 30])
xtickangle(0)

legend([s,p],'incongruent data','fit \it{(m x + b)}','Location','best')

title('\textbf{B}','Position',[-4.7 5.8 0],'HorizontalAlignment','center')

set(gca,'FontSize',pp.fsize-2,'LineWidth',pp.linewidth,'Layer','top','Box','on')
xlabel('ICLVL / semitones','FontWeight','normal','FontSize',pp.fsize)
ylabel('Sound pleasantness rating','FontWeight','normal','FontSize',pp.fsize)


if isfield(pp, 'figwidth')
    if ~isempty(pp.figwidth)
        set(figh3, 'PaperPositionMode', 'manual');
        set(figh3, 'PaperUnits', 'centimeters');
        set(figh3, 'PaperPosition', [0 0 pp.figwidth 7]);
        set(figh3, 'Toolbar', 'none')
        set(figh3, 'Menubar', 'none')
        set(figh3, 'PaperSize', [pp.figwidth 7])
    end
end

print(figh3, [fig_folder filesep 'SQRcombined.pdf'], '-dpdf');
