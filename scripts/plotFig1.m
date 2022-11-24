clear
close all

load('/home/simon/ma/code/soundQuality/data/expData_norm.mat');
run plot_properties.m

expData = expData_norm;

%% LME
% including "all" condition, using full pitch range, but no interactions
% with quadratic pitch term

numParts = size(expData.quality.Violin,1);
% convert to long format, so that lme can read data
for iInstr = 1:nInstr
    instrName = instruments{iInstr};
    X = [];
    for nP = 1:numParts
        for nS = 1:19
            for nC = 1:4 %% !!!
                X = [X; nP, nS, nC, expData.quality.(instrName)(nP, nS, nC)];
            end
        end
    end
    expData.quality_long.(instrName) = X;
end

% convert into table format
clear tab
for iInstr = 1:nInstr
    instrName = instruments{iInstr};
    Y = expData.quality_long.(instrName);
    tab.(instrName) = table(categorical(Y(:,1)), ...
        (Y(:,2)-1)/3, categorical(Y(:,3)), Y(:,4)); % scale pitches to be in octave units
    tab.(instrName).Properties.VariableNames = {'part', 'pitch', 'cond', 'resp'};
end

% compute lme models
for iInstr = 1:nInstr
    instrName = instruments{iInstr};
    exp2.quality.(instrName) = fitlme(tab.(instrName), 'resp ~ 1 + pitch*cond + pitch:pitch + (1 | part)', ...
        'FitMethod', 'REML', 'CheckHessian', 1, 'DummyVarCoding', 'effects')
    exp2.quality.(instrName).Rsquared
end
% DD = designMatrix(exp2.quality.(instrName{nInstr})); % design matrix of


% get LME results in proper format
for iInstr = 1:nInstr
    instrName = instruments{iInstr};
    [beta, betname, s] = fixedEffects(exp2.quality.(instrName));
    ss = table(s.Name, num2str(round(s.Estimate,2)), num2str(round(s.Lower,2)), num2str(round(s.Upper,2)), ...
        num2str(round(s.tStat,2)), num2str(round(s.pValue,3)));
    ss.Properties.VariableNames = {'Name', '\beta', 'CI low', 'CI high', 't-value', 'p-value'};
    disp(instrName)
    ss
end

%%

figh1 = figure('visible', pp.visible,...
    'DefaultTextFontName', pp.fname,...
    'DefaultAxesFontName', pp.fname);
axh1 = axes('Parent', figh1);
hold(axh1, 'on');
box(axh1, 'off');

tl = tiledlayout(1,4,'TileSpacing','loose','Padding','compact');
c = lines(3);
c(4,:) = [0 0 0];
m = {'v','o','^','s'};
t = {'A','B','C','D'};
pOrder = [0 3 6 9];
mSize = [2 2 2 2];

regBounds.Violin = [7.25 20];
regBounds.VocalAlto = [6.75 13.75];
regBounds.ClarinetBb = [5.5 16.25];
regBounds.Tuba = [0 10.25];

for iInstr = 1:nInstr
    instrName = instruments{iInstr};
    pit_ind.(instrName) = find(~isnan(expData.quality.(instrName)(1, :, 4)));
end
stepSiz = 19;
f.cond = [1 0 0; 0 1 0; 0 0 1; -1 -1 -1];

% compute lme models
for iInstr = 1:nInstr
    instrName = instruments{iInstr};
    exp2.quality.(instrName) = fitlme(tab.(instrName), 'resp ~ 1 + pitch*cond + pitch:pitch + (1 | part)', ...
        'FitMethod', 'REML', 'CheckHessian', 1, 'DummyVarCoding', 'effects')
    exp2.quality.(instrName).Rsquared
end

for iInstr = 1:nInstr
    instrName = instruments{iInstr}

    nexttile
    hold on

    % find ~nan idx of 'all'
    idx1 = find(~isnan(expData.quality.(instrName)(1,:,4)));

    % plot instrument range region
%     patch([regBounds.(instrName)(1) idx1 regBounds.(instrName)(2) regBounds.(instrName)(2) fliplr(idx1) regBounds.(instrName)(1)],[ones(1,numel(idx1)+2)*6 ones(1,numel(idx1)+2)],...
%         'k',...
%         'FaceAlpha',0.04,...
%         'EdgeColor','none')

    plot([regBounds.(instrName)(1) regBounds.(instrName)(1)],[1 6],...
        'Color',[.5 .5 .5 .4],...
        'LineStyle','-',...
        'LineWidth',pp.linewidth-0.5)
    plot([regBounds.(instrName)(2) regBounds.(instrName)(2)],[1 6],...
        'Color',[.5 .5 .5 .4],...
        'LineStyle','-',...
        'LineWidth',pp.linewidth-0.5)

%     % plot mean region
%     patch([22.5 23.5 23.5 22.5],[6 6 1 1],...
%         'k',...
%         'FaceAlpha',0.04,...
%         'EdgeColor','none')

    % plot middle line
    plot([0 23],[3.5 3.5],...
        'Color',[0.5 0.5 0.5],...
        'LineStyle',':',...
        'LineWidth',pp.linewidth-0.3)

    % plot register locations
    for ii = 1:3
        plot([r.(instrName)(ii) r.(instrName)(ii)],[1 6],...
            'Color',[c(ii,:) 0.4],...
            'LineWidth',pp.linewidth-0.5,...
            'LineStyle','-')
    end

    % plot actual data
    for iCond = 1:nCond
        condName = conditions{iCond};

        % compute LME
        f.pitch = linspace(0, 6, stepSiz)'; % factor pitch
        if iCond == 4
            % rr = (registers.(instrName{nInstr})-1)/3; % a little too small of a range...
            f.pitch = f.pitch(pit_ind.(instrName));
        end
        Xdes = [ones(length(f.pitch),1), f.pitch, ...
            ones(length(f.pitch),1)*f.cond(iCond, :), f.pitch*f.cond(iCond,:),...
            f.pitch.^2];%, (f.pitch.^2)*f.cond(nC,:)]; % design matrix
        [beta,~,stats] = fixedEffects(exp2.quality.(instrName));
        yhat = Xdes*beta;

        % compute CI
        idx2 = find(~isnan(expData.quality.(instrName)(1,:,iCond)));

        Y = mean(expData.quality.(instrName)(:,idx2,iCond))';
        [R,P] = corrcoef(Y,yhat);
        R2(iCond) = round(R(2),2);
        pval(iCond) = round(P(2),3);

        % compute CI for individual data points
        nNote = numel(idx2);
        for iNote = 1:nNote
            CI(:,iNote) = bootci(1000,@mean,expData.quality.(instrName)(:,idx2(iNote),iCond));

        end

        % compute CI for mean
        CImean = bootci(1000,@mean,reshape(expData.quality.(instrName)(:,idx1,iCond),numel(expData.quality.(instrName)(:,idx1,iCond)),1));
        negE = mean(mean(expData.quality.(instrName)(:,idx1,iCond)),'omitnan') - CImean(2);
        posE = CImean(1) - mean(mean(expData.quality.(instrName)(:,idx1,iCond)),'omitnan');

%         confInt(iCond) = patch([idx2 fliplr(idx2)],[CI(1,:) fliplr(CI(2,:))],...
%             c(iCond,:),...
%             'FaceAlpha',0.1,...
%             'EdgeColor','none');

        lme(iCond) = plot(idx2,yhat,...
            'Color',c(iCond,:),...
            'LineWidth',pp.linewidth);

%         uistack(lme(iCond),'down',pOrder(iCond))

%         p(iCond) = plot(1:19,mean(expData.quality.(instrName)(:,:,iCond)),...
%             'Color',c(iCond,:),...
%             'LineWidth',pp.linewidth);
%         pause(1)
%         hMarkers = p(iCond).MarkerHandle;
%         hMarkers.EdgeColorData = uint8(255*[c(iCond,:)';0.5]);
%         hMarkers.EdgeColorType = 'truecoloralpha';

%         set(hMarkers,'EdgeColorData',uint8(255*[c(iCond,:)';0.5]))


%         plot(22,mean(mean(expData.quality.(instrName)(:,idx1,iCond)),'omitnan'),...
%             'Color',c(iCond,:),...
%             'Marker',m{iCond},...
%             'MarkerSize',pp.markersize-1,...
%             'LineWidth',pp.linewidth)

        clear CI

    end

    hold off

    xlim([-1 21])
    xticks([1,4,7,10,13,16,19])
    xticklabels({'1','2','3','4','5','6','7'})
    xtickangle(0)

    ylim([1 6])
    yticks([1 3.5 6])

    if iInstr == 4
        legh = legend(lme,{'LRA','MRA','HRA','congr.'},...
            'FontSize',pp.fsize-2,...
            'NumColumns',4,...
            'Location','best');
        legh.ItemTokenSize = [12 1 0];

    end

%     if iInstr == 3
%         legh = legend([p(4) confInt(4)],{'emp','95% CI'},...
%             'FontSize',pp.fsize-2,...
%             'NumColumns',1,...
%             'Location','northeast');
%         legh.ItemTokenSize = [12 12 0];
%     end

%     text(0.25,5.75,instrName,'FontSize',pp.fsize-2,'FontWeight','bold')
% 
%     title(t{iInstr},...
%         'Position',[-2.3 6 0],...
%         'HorizontalAlignment','left')
    title(instrumentNames{iInstr})

    set(gca,'FontSize',pp.fsize-2,'LineWidth',pp.linewidth,'Layer','top')

    xlabel(tl,'Pitch re F#','FontSize',pp.fsize+1,'FontWeight','normal')
    ylabel(tl,'Pleasantness rating','FontSize',pp.fsize+1,'FontWeight','normal')

    if iInstr == 4
        legh.Position(1) = legh.Position(1)-0.45;
        legh.Position(2) = legh.Position(2)-0.09;
    end
%     if iInstr == 3
%         legh.Position(1) = legh.Position(1)+0.135;
%         legh.Position(2) = legh.Position(2)+0.025;
%     end

    if isfield(pp, 'figwidth')
        if ~isempty(pp.figwidth)
            set(figh1, 'PaperPositionMode', 'manual');
            set(figh1, 'PaperUnits', 'centimeters');
            set(figh1, 'PaperPosition', [0 0 pp.figwidth 4.5]);
            set(figh1, 'Toolbar', 'none')
            set(figh1, 'Menubar', 'none')
            set(figh1, 'PaperSize', [pp.figwidth 4.5])
        end
    end

    R2
    pval
    clear R2 p

end

print(figh1, [fig_folder filesep 'pleasantness_model.pdf'], '-dpdf');