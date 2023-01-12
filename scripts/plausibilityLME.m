clear
close all

load('data/expData_norm.mat');
expData = expData_norm;
run plot_properties.m

%% LME

nPart = size(expData.quality.Violin,1);
regIdx.Violin = [2 4 3];        % what should be the reference register?
regIdx.VocalAlto = [2 4 3];
regIdx.ClarinetBb = [2 4 3];
regIdx.Tuba = [1 4 3];

for iInstr = 1:nInstr
    instrName = instruments{iInstr};
    X = [];

    for iCond = 1:nCond

        nonnan = ~isnan(expData.plausib.(instrName)(1,:,iCond));
        tempData = expData.plausib.(instrName)(:,nonnan,iCond);

        if iCond < 4
            pData.(instrName)(:,:,iCond) = tempData(:,regIdx.(instrName));
        else
            pData.(instrName)(:,:,iCond) = tempData;
        end

    end

    x = [];
    for iPart = 1:nPart
        for iStim = 1:3
            for iCond = 1:nCond

                X = [X; iPart, iStim, iCond, pData.(instrName)(iPart,iStim,iCond)];

            end
        end
    end
    expData.plausib_long.(instrName) = X;
end

% convert into table format
clear tab
for iInstr = 1:nInstr
    instrName = instruments{iInstr};
    Y = expData.plausib_long.(instrName);
    tab.(instrName) = table(categorical(Y(:,1)), ...
        categorical(Y(:,2)), categorical(Y(:,3)), Y(:,4)); % scale pitches to be in octave units
    tab.(instrName).Properties.VariableNames = {'part', 'reg', 'cond', 'resp'};
end

% compute lme models
for iInstr = 1:nInstr
    instrName = instruments{iInstr};
    exp2.plausib.(instrName) = fitlme(tab.(instrName), 'resp ~ 1 + cond*reg + (1|part)', ...
        'FitMethod', 'REML', 'DummyVarCoding', 'effects')
    exp2.plausib.(instrName).Rsquared
end

% get LME results in proper format
for iInstr = 1:nInstr
    instrName = instruments{iInstr};
    [beta, betname, s] = fixedEffects(exp2.plausib.(instrName));
    ss = table(s.Name, num2str(round(s.Estimate,2)), num2str(round(s.Lower,2)), num2str(round(s.Upper,2)), ...
        num2str(round(s.tStat,2)), num2str(round(s.pValue,3)));
    ss.Properties.VariableNames = {'Name', '\beta', 'CI low', 'CI high', 't-value', 'p-value'};
    disp(instrName)
    ss
    anova(exp2.plausib.(instrName))
end