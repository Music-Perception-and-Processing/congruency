clear
close all

load('data/expData_norm.mat');
expData = expData_norm;
run plot_properties.m

condNames = {'LRSE','MRSE','HRSE','congr.'};
regNames = {'low','mid','high'};

%% LME

nPart = size(expData.plausib.Violin,1);
regIdx.Violin = [2 3 4];        % what should be the reference register?
regIdx.VocalAlto = [2 3 4];
regIdx.ClarinetBb = [2 3 4];
regIdx.Tuba = [1 3 4];

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
        categorical(Y(:,2)), categorical(Y(:,3)), Y(:,4));
    tab.(instrName).Properties.VariableNames = {'part', 'reg', 'cond', 'resp'};

    f1 = categorical(Y(:,1));
    f2 = categorical(Y(:,2));
    f2 = renamecats(f2,{'low';'mid';'high'});
    f3 = categorical(Y(:,3));
    f3 = renamecats(f3,{'LRSE';'MRSE';'HRSE';'congr.'});
    f4 = Y(:,4);
    f5 = categorical(iInstr.*ones(456,1));
    f5 = renamecats(f5,instrumentNames{iInstr});
    plausibility.(instrName) = table(f1,f2,f3,f4);
    plausibility.(instrName).Properties.VariableNames = {'participants',...
        'register','condition','response'};

%     if iInstr == 1
%         plausibility = table(f1,f2,f3,f4,f5);
%         plausibility.Properties.VariableNames = {'participants','register',...
%             'condition','response','instrument'};
%     else
%         pTemp = table(f1,f2,f3,f4,f5);
%         pTemp.Properties.VariableNames = {'participants','register',...
%             'condition','response','instrument'};
%         plausibility = [plausibility;pTemp];
%     end
% 
%     plausibility.Properties.VariableNames = {'participants','register',...
%         'condition','response','instrument'};

end

% compute lme models
for iInstr = 1:nInstr
    instrName = instruments{iInstr}
    exp2.plausib.(instrName) = fitlme(tab.(instrName), 'resp ~ 1 + reg*cond + (1|part)', ...
        'FitMethod', 'REML', 'DummyVarCoding', 'effects');
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