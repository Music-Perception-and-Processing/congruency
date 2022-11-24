clear
close all

load('/home/simon/ma/code/soundQuality/data/trueSE.mat')

run plot_properties.m

pitches = {'F#1','A#1','D2',...
    'F#2','A#2','D3',...
    'F#3','A#3','D4',...
    'F#4','A#4','D5',...
    'F#5','A#5','D6',...
    'F#6','A#6','D7',...
    'F#7'};

% setup ERB filter bank
fs = 44100;
f = 1:fs/2;
numBands = 128;
nCC = 13;
[fb,cf,bw] = designAuditoryFilterBank(fs, ...
    'FrequencyScale','erb', ...
    'FFTLength',fs-1, ...
    'NumBands',numBands, ...
    'FrequencyRange',[20 12000], ...
    'Normalization','bandwidth');

% plotting
figh1 = figure('visible', pp.visible,...
    'DefaultTextFontName', pp.fname,...
    'DefaultAxesFontName', pp.fname);
axh1 = axes('Parent', figh1);
hold(axh1, 'on');
box(axh1, 'off');

tl = tiledlayout(2,3,'TileSpacing','loose','Padding','compact');
xlabel(tl,'Frequency / kHz','FontSize',pp.fsize+1)
ylabel(tl,'Level / dB','FontSize',pp.fsize+1)
c = lines(3);
c(4,:) = [0.5 0.5 0.5];

locs = [7 10 13];
f = muspitch2freq(pitches(locs));
fNat = 20:12000;

c4 = [1 1 1];

for ii = 1:3
    nexttile
    hold on

    interpSE = interp1(cf,trueSE.ClarinetBb(7,:),fNat,'spline');
    harmonics = round((1:60).*f(ii));
    amps = zeros(1,numel(fNat));
    amps(ismember(fNat,harmonics)) = 1;
    amps = interpSE.*amps;
    idx = find(amps);

    for jj = 1:numel(idx)
        ph = plot([harmonics(jj) harmonics(jj)],[-90 amps(idx(jj))],...
            'Color',[c(1,:) 0.5],...
            'LineWidth',pp.linewidth);

    end

    if ii == 2
        pc = plot(cf,trueSE.ClarinetBb(10,:),...
            'Color',c(2,:),...
            'LineWidth',pp.linewidth,...
            'LineStyle','-');
    end

    if ii == 3
        pc = plot(cf,trueSE.ClarinetBb(13,:),...
            'Color',c(3,:),...
            'LineWidth',pp.linewidth,...
            'LineStyle','-');
    end

    pi = plot(cf,trueSE.ClarinetBb(7,:),...
        'Color',[c(1,:) c4(ii)],...
        'LineWidth',pp.linewidth);
    hold off

    if ii == 1
        legh = legend([pi ph],{'LRA-SE','H(LRA)'},...
            'Location','best');
    elseif ii == 2
        legh = legend([pc ph],{'MRA-SE','H(MRA)'},...
            'Location','best');
    else
        legh = legend([pc ph],{'HRA-SE','H(HRA)'},...
            'Location','best');
    end
    legh.ItemTokenSize = [10 1 0];
%     legh.Position(1) = legh.Position(1) - 0.0000001;
%     legh.Position(2) = legh.Position(2) - 0.0000001;

    if ii == 1
        title('A','Position',[5 4 0],'HorizontalAlignment','center')
    end

    set(gca,'FontSize',pp.fsize,...
        'LineWidth',pp.linewidth,...
        'Layer','top',...
        'XScale','log')

    xlim([20 12000])
    xticks([100 1000 10000])
    xticklabels({'0.1','1','10'})

    ylim([-90 10])
    yticks([-90 -60 -30 0])

    box off

end

locs = [10 13 16];
f = muspitch2freq(pitches(locs));
fNat = 20:12000;

for ii = 1:3
    nexttile
    hold on

    interpSE = interp1(cf,trueSE.ClarinetBb(10,:),fNat,'spline');
    harmonics = round((1:60).*f(ii));
    amps = zeros(1,numel(fNat));
    amps(ismember(fNat,harmonics)) = 1;
    amps = interpSE.*amps;
    idx = find(amps);

    for jj = 1:numel(idx)
        ph = plot([harmonics(jj) harmonics(jj)],[-90 amps(idx(jj))],...
            'Color',[c(2,:) 0.5],...
            'LineWidth',pp.linewidth);

    end

    if ii == 2
        pc = plot(cf,trueSE.ClarinetBb(13,:),...
            'Color',c(3,:),...
            'LineWidth',pp.linewidth,...
            'LineStyle','-');
    end

    if ii == 3
        pc = plot(cf,trueSE.ClarinetBb(16,:),...
            'Color',c(4,:),...
            'LineWidth',pp.linewidth,...
            'LineStyle','-');
    end

    pi = plot(cf,trueSE.ClarinetBb(10,:),...
        'Color',[c(2,:) c4(ii)],...
        'LineWidth',pp.linewidth);
    hold off

    if ii == 1
        legh = legend([pi ph],{'MRA-SE','H(MRA)'},...
            'Location','best');
    elseif ii == 2
        legh = legend([pc ph],{'HRA-SE','H(HRA)'},...
            'Location','best');
        legh.Position(1) = legh.Position(1) - 0.07;
        legh.Position(2) = legh.Position(2) - 0.196;
    else
        legh = legend([pc ph],{'F#6-SE','H(F#6)'},...
            'Location','best');
        legh.Position(1) = legh.Position(1) - 0.071;
        legh.Position(2) = legh.Position(2) - 0.196;
    end
    legh.ItemTokenSize = [10 1 0];
%     legh.Position(1) = legh.Position(1) - 0.05;
%     legh.Position(2) = legh.Position(2) - 0.05;

    if ii == 1
        title('B','Position',[5 4 0],'HorizontalAlignment','center')
    end

    set(gca,'FontSize',pp.fsize,...
        'LineWidth',pp.linewidth,...
        'Layer','top',...
        'XScale','log')

    xlim([20 12000])
    xticks([100 1000 10000])
    xticklabels({'0.1','1','10'})

    ylim([-90 10])
    yticks([-90 -60 -30 0])

    box off

end

if isfield(pp, 'figwidth')
    if ~isempty(pp.figwidth)
        set(figh1, 'PaperPositionMode', 'manual');
        set(figh1, 'PaperUnits', 'centimeters');
        set(figh1, 'PaperPosition', [0 0 pp.figwidth 7.5]);
        set(figh1, 'Toolbar', 'none')
        set(figh1, 'Menubar', 'none')
        set(figh1, 'PaperSize', [pp.figwidth 7.5])
    end
end

print(figh1, [fig_folder filesep 'clarinet_example.pdf'], '-dpdf');