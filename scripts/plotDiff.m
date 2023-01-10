clear
close all

load('data/expData_norm.mat');
expData = expData_norm;
clear expData_norm

load('data/old_expData_norm.mat');
expData_old = expData_norm;
clear expData_norm

run plot_properties.m

%% plotting

tl = tiledlayout(2,2);
c = lines(3);
c(4,:) = [0 0 0];
m = {'v','o','^','s'};

for iInstr = 1:nInstr
    instrName = instruments{iInstr};

    nexttile
    hold on

    for iCond = 1:nCond
        condName = conditions{iCond};
        p(iCond) = plot(1:19,mean(expData.quality.(instrName)(:,:,iCond))-mean(expData_old.quality.(instrName)(:,:,iCond)),...
            'Color',c(iCond,:));

        plot(21,mean(mean(expData.quality.(instrName)(:,:,iCond))-mean(expData_old.quality.(instrName)(:,:,iCond)),'omitnan'),...
            'Color',c(iCond,:),...
            'Marker',m{iCond});

    end

    hold off

    xlim([-1 22])
    xticks([1,4,7,10,13,16,19])
    xticklabels({'1','2','3','4','5','6','7'})
    xtickangle(0)

    ylim([-1 1])
    yticks([-1 -0.5 0 0.5 1])

    title(instrumentNames{iInstr})

end

xlabel(tl,'Pitch re F#','FontSize',pp.fsize+1,'FontWeight','normal')
ylabel(tl,'Pleasantness rating','FontSize',pp.fsize+1,'FontWeight','normal')