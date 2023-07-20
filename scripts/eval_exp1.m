clear
close all

%% EXPERIMENT 1
load('data/expData1.mat')
load('data/F0.mat')
load('data/ICLVL.mat')
load('data/grid.mat')

% factors:  

% inner space vs. outer space 
% congruency  -- alternatively: ICLVL 
% pitch? 


%% EXPERIMENT 1 - PLEASANTNESS RATINGS 
numParts = 41; 
congruency = {'congruentF0', 'incongruentF0', 'fixedF0'}; 
space = {'inner', 'outer'}; 

% get incongruency level into format 
for nS = 1:2
    IC.incongruentF0.(space{nS}) = ICLVL.(space{nS})(:,:);
    IC.congruentF0.(space{nS}) = 0*ICLVL.(space{nS})(:,:);
    IC.fixedF0.(space{nS}) = nan*ICLVL.(space{nS})(:,:); 
end

X = []; 
for nCongr = 1:length(congruency)
    for nP = 1:numParts 
        for nS = 1:2 % space
            if nS == 1
                xx = expData1.quality.(congruency{nCongr}).(space{nS})(:,:, nP);
                f0 = F0.quality.(congruency{nCongr}).(space{nS})(:,:); 
                iclvl = IC.(congruency{nCongr}).(space{nS})(:,:); 
            elseif nS == 2
                xx = expData1.quality.(congruency{nCongr}).(space{nS})(:,nP);
                f0 = F0.quality.(congruency{nCongr}).(space{nS})(:,:); 
                iclvl = IC.(congruency{nCongr}).(space{nS})(:,:); 
            end
            X = [X; [ones(size(xx(:),1),1)*[nP, nS, nCongr], xx(:), f0(:), abs(iclvl(:))]];  % no structural analysis of pc dimensions 
        end
    end
end


%% convert into table format 
clear tab
Y = X; 
tab.exp1 = table(categorical(Y(:,1)), categorical(Y(:,2)), categorical(Y(:,3)), Y(:,4), log2(Y(:,5)));  % scale pitches to be in octave units 
tab.exp1.Properties.VariableNames = {'part', 'space', 'congr', 'resp', 'pitch'};

% compute lme model 
exp1_quality_2 = fitlme(tab.exp1, 'resp ~ 1 + space*congr + (1 + space + congr | part)', ...
'FitMethod', 'REML', 'CheckHessian', 1, 'DummyVarCoding', 'effects')
exp1_quality_2.Rsquared

% but when pitch is included, effect of congr get very small (only in
% interaction term) 
exp1_quality_3 = fitlme(tab.exp1, 'resp ~ 1 + space*congr + pitch^2 + (1 + space + congr | part)', ...
'FitMethod', 'REML', 'CheckHessian', 1, 'DummyVarCoding', 'effects')
exp1_quality_3.Rsquared

% get LME results in proper format 
[beta, betname, s] = fixedEffects(exp1_quality_3);
ss = table(s.Name, num2str(round(s.Estimate,2)), num2str(round(s.Lower,2)), num2str(round(s.Upper,2)), ...
    num2str(round(s.tStat,2)), num2str(round(s.pValue,3)));
ss.Properties.VariableNames = {'Name', '\beta', 'CI low', 'CI high', 't-value', 'p-value'}; 
ss


%% consider only inner space without fixed F0 condition 

clear tab
indd = ~isnan(X(:,6)) & (X(:,2) == 1); 
Y = X(indd,:); 
tab.exp1 = table(categorical(Y(:,1)), categorical(Y(:,2)), categorical(Y(:,3)), Y(:,4), log2(Y(:,5)), Y(:,6));  % scale pitches to be in octave units 
tab.exp1.Properties.VariableNames = {'part', 'space', 'congr', 'resp', 'pitch', 'iclvl'};

exp1_quality_2 = fitlme(tab.exp1, 'resp ~ 1 + iclvl + congr + pitch^2 + (1 | part)', ...
'FitMethod', 'REML', 'CheckHessian', 1, 'DummyVarCoding', 'effects')
exp1_quality_2.Rsquared


%% EXPERIMENT 1 - BRIGHTNESS RATINGS 
numParts = 41; 
congruency = {'congruentF0', 'incongruentF0', 'fixedF0'}; 
% space = {'inner', 'outer'}; 

pc = {'pc1', 'pc2'}; 
bright_ind.pc1 = [1:9];
bright_ind.pc2 = [10:16]; 
clear X
for nPC = 1:2
    X.(pc{nPC}) = []; 
    for nCongr = 1:length(congruency)
        for nP = 1:numParts 
            xx = expData1.brigntness.(congruency{nCongr})(bright_ind.(pc{nPC}),nP);
            f0 = F0.brightness.(congruency{nCongr})(bright_ind.(pc{nPC}));
            PC_grid = grid.brightness(bright_ind.(pc{nPC}), nPC); 
            X.(pc{nPC}) = [X.(pc{nPC}); [ones(size(xx(:),1),1)*[nP, nCongr, nPC], xx(:), f0(:), PC_grid]];  % no structural analysis of pc dimensions 
        end
    end
end


%% compute models 
clc 
% FIXED F0
for nPC = 1:2
    clear tab
    indd = (X.(pc{nPC})(:, 2) == 3); 
    Y = X.(pc{nPC})(indd, :); 
    tab.exp1 = table(categorical(Y(:,1)), categorical(Y(:,2)), categorical(Y(:,3)), Y(:,4), log2(Y(:,5)), Y(:,6));  % scale pitches to be in octave units 
    tab.exp1.Properties.VariableNames = {'part', 'congr', 'PCdim', 'resp', 'F0', 'PC'};

    % compute lme model 
    exp1_bright = fitlme(tab.exp1, 'resp ~ 1 + PC + (1 | part)', ...
    'FitMethod', 'REML',  'DummyVarCoding', 'effects');
    exp1_bright.Rsquared;

    % get LME results in proper format 
    [beta, betname, s] = fixedEffects(exp1_bright);
    ss = table(s.Name, num2str(round(s.Estimate,2)), num2str(round(s.Lower,2)), num2str(round(s.Upper,2)), ...
        num2str(round(s.tStat,2)), num2str(round(s.pValue,3)));
    ss.Properties.VariableNames = {'Name', '\beta', 'CI low', 'CI high', 't-value', 'p-value'}; 
    ss
end

% congruent / incongruent F0
for nPC = 1:2
    clear tab
    indd = ~(X.(pc{nPC})(:, 2) == 3); 
    Y = X.(pc{nPC})(indd, :); 
    tab.exp1 = table(categorical(Y(:,1)), categorical(Y(:,2)), categorical(Y(:,3)), Y(:,4), log2(Y(:,5)), Y(:,6));  % scale pitches to be in octave units 
    tab.exp1.Properties.VariableNames = {'part', 'congr', 'PCdim', 'resp', 'F0', 'PC'};

    % compute lme model 
    exp1_bright = fitlme(tab.exp1, 'resp ~ 1 + PC*congr + F0 + (1 | part)', ...
    'FitMethod', 'REML',  'DummyVarCoding', 'effects');
    exp1_bright.Rsquared;

    % get LME results in proper format 
    [beta, betname, s] = fixedEffects(exp1_bright);
    ss = table(s.Name, num2str(round(s.Estimate,2)), num2str(round(s.Lower,2)), num2str(round(s.Upper,2)), ...
        num2str(round(s.tStat,2)), num2str(round(s.pValue,3)));
    ss.Properties.VariableNames = {'Name', '\beta', 'CI low', 'CI high', 't-value', 'p-value'}; 
    ss
end



%% EXPERIMENT 2
load('data/expData_norm.mat')
expData = expData_norm;
load('data/pitches.mat')
load('data/registers.mat')

% expData = P x S x C
% P: participants (N=29)
% S: stimuli (N=19)
% C: conditions (N=4) -> low, mid, hig, all

% big picture 
instrName = {'Violin','VocalAlto','ClarinetBb','Tuba'};
%parts_incl = [1:20, 22:26, 28:29]; % if deemed outliers 
for nInstr = 1:length(instrName)
    %figure; plot([1:19], squeeze(nanmean(expData.quality.(instrName{nInstr})(parts_incl,:,:))), 's-')
    figure; plot([1:19], squeeze(nanmean(expData.quality.(instrName{nInstr}))), 's-')
    title(instrName{nInstr})
    set(gca, 'XTick', [1:3:19])
    set(gca, 'XTickLabel', {pitches{[1:3:19]}}) 
end


%% check behavior of individual participants, averaged across instruments
% and conditions 
datP = cat(3, squeeze(nanmean(expData.quality.(instrName{1})(:,:,:),3)), ...
    squeeze(nanmean(expData.quality.(instrName{2})(:,:,:),3)), ...
    squeeze(nanmean(expData.quality.(instrName{3})(:,:,:),3)), ...
    squeeze(nanmean(expData.quality.(instrName{4})(:,:,:),3))); 

X = nanmean(datP, 3)';     
for nP = 1:29
    figure; 
    plot([1:19], X(:,nP), 's-'); hold on;
    set(gca, 'XTick', [1:3:19])
    set(gca, 'XTickLabel', {pitches{[1:3:19]}}) 
    ylim([1 6]) 
end



%% LME
% including "all" condition, using full pitch range, but no interactions
% with quadratic pitch term 

numParts = size(expData.quality.(instrName{nInstr}),1);
% convert to long format, so that lme can read data 
for nInstr = 1:length(instrName)
    X = [];
    for nP = 1:numParts 
        for nS = 1:19
            pitch = muspitch2freq(pitches{nS});
            for nC = 1:4 %% !!!
                 X = [X; nP, nS, nC, expData.quality.(instrName{nInstr})(nP, nS, nC), nInstr, pitch];
            end
        end
    end
    expData.quality_long.(instrName{nInstr}) = X; 
end

% convert into table format 
clear tab
sumf4 = [];
for nInstr = 1:4
    Y = expData.quality_long.(instrName{nInstr}); 
    tab.(instrName{nInstr}) = table(categorical(Y(:,1)), ...
        (Y(:,2)-1)/3, categorical(Y(:,3)), Y(:,4)); % scale pitches to be in octave units 
    tab.(instrName{nInstr}).Properties.VariableNames = {'part', 'pitch', 'cond', 'resp'};

    % make combined table
    % participants(38)
    % pitch(19)
    % condition(4)
    % instrument(4)
    % response
    f1 = categorical(Y(:,1));
    f2 = Y(:,6);
    f3 = categorical(Y(:,3));
    f3 = renamecats(f3,{'LRSE','MRSE','HRSE','congr.'});
    f3 = reordercats(f3,{'congr.','LRSE','MRSE','HRSE'});
    f4 = categorical(Y(:,5));
    sumf4 = [sumf4;f4];
    f5 = Y(:,4);

    if nInstr == 1
        pleasantness = table(f1,f2,f3,f4,f5);
        pleasantness.Properties.VariableNames = {'participants','pitch',...
            'condition','instrument','response'};
    else
        pTemp = table(f1,f2,f3,f4,f5);
        pTemp.Properties.VariableNames = {'participants','pitch',...
            'condition','instrument','response'};
        pleasantness = [pleasantness;pTemp];
    end

    pleasantness.Properties.VariableNames = {'participants','pitch',...
        'condition','instrument','response'};
end
sumf4 = table(renamecats(sumf4,{'violin','vocal','clarinet','tuba'}));
sumf4.Properties.VariableNames = {'instrument'};
pleasantness(:,4) = sumf4;

% compute lme models 
for nInstr = 1:4
    exp2.quality.(instrName{nInstr}) = fitlme(tab.(instrName{nInstr}), 'resp ~ 1 + cond + pitch + pitch:pitch + (1 | part)', ...
    'FitMethod', 'REML', 'CheckHessian', 1, 'DummyVarCoding', 'reference')
    exp2.quality.(instrName{nInstr}).Rsquared
end
% DD = designMatrix(exp2.quality.(instrName{nInstr})); % design matrix of


% get LME results in proper format 
for nInstr = 1:4
    [beta, betname, s] = fixedEffects(exp2.quality.(instrName{nInstr}));
    ss = table(s.Name, num2str(round(s.Estimate,2)), num2str(round(s.Lower,2)), num2str(round(s.Upper,2)), ...
        num2str(round(s.tStat,2)), num2str(round(s.pValue,3)));
    ss.Properties.VariableNames = {'Name', '\beta', 'CI low', 'CI high', 't-value', 'p-value'}; 
    disp(instrName{nInstr})
    ss
end


%% plot fixed effects predictions 
clear pit_ind
for nInstr = 1:length(instrName)
    pit_ind.(instrName{nInstr}) = find(~isnan(expData.quality.(instrName{nInstr})(1, :, 4))); 
end

stepSiz = 19; 
f.cond = [1 0 0; % row: condition level; colums: dummy coding
          0 1 0;
          0 0 1;
        -1 -1 -1]; 

% hit it!
figure
cmap = colormap('lines'); 
cmap(4,:) = [0.3 0.3 0.3];
for nInstr = 1:4
    subplot(1,4,nInstr)
    clear Xdes
    for nC = 1:4
        f.pitch = linspace(0, 6, stepSiz)'; % factor pitch
        if nC == 4
           % rr = (registers.(instrName{nInstr})-1)/3; % a little too small of a range... 
           f.pitch = f.pitch(pit_ind.(instrName{nInstr})); 
        end
        Xdes = [ones(length(f.pitch),1), f.pitch, ...
            ones(length(f.pitch),1)*f.cond(nC, :), f.pitch*f.cond(nC,:),... 
            f.pitch.^2];%, (f.pitch.^2)*f.cond(nC,:)]; % design matrix 
        beta = fixedEffects(exp2.quality.(instrName{nInstr})); 
        yhat = Xdes*beta;
        %yhat = predict(exp2.quality.(instrName{nInstr}), 'Conditional',0);
        plot(f.pitch, yhat, 'color', cmap(nC,:),'linewidth', 4); hold on
        plot(([1:19]-1)/3, squeeze(nanmean(expData.quality.(instrName{nInstr})(:,:,nC))), '-', 'color', cmap(nC,:), 'linewidth', 1); hold on
        plot(([1:19]-1)/3, squeeze(nanmean(expData.quality.(instrName{nInstr})(:,:,nC))), 's', 'color', cmap(nC,:), 'linewidth', 1)
        
    end
    title(instrName{nInstr})
    xlim([-.25, 6.25])
    ylim([1 6])
end


