clear
close all

load('data/expData_norm.mat')
load('data/amps.mat')

run scripts/plot_properties.m

figure
tl = tiledlayout(1,4)
for iInstr = 1:nInstr

    instrName = instruments{iInstr};

    nexttile
    hold on

    for iCond = 1:nCond

        condName = {iCond};
        ratio = 20.*log10(amps.(instrName)(:,1,iCond)./amps.(instrName)(:,2,iCond));

        scatter(ratio,mean(expData_norm.quality.(instrName)(:,:,iCond),'omitnan'))

    end

    hold off
    title(instrumentNames{iInstr})
    xlim([-80 60])
    ylim([0 6])

end

xlabel(tl,'Ratio H1/H2 / dB')
ylabel(tl,'Sound pleasantness rating')