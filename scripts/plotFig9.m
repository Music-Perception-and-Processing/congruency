clear
close all

run plot_properties.m

set(0,'defaultAxesFontName',pp.fname)
set(0,'defaultTextInterpreter','tex')
set(0,'defaultAxesTickLabelInterpreter','tex')
set(0,'defaultLegendInterpreter','tex')

load('data/trueSE.mat')

addpath('functions/')

pitches = {'F#1','A#1','D2',...
    'F#2','A#2','D3',...
    'F#3','A#3','D4',...
    'F#4','A#4','D5',...
    'F#5','A#5','D6',...
    'F#6','A#6','D7',...
    'F#7'};

% setup ERB filter bank
fs = 44100;
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
box(axh1, 'on');

tl = tiledlayout(3,4,'TileSpacing','loose','Padding','compact');
xlabel(tl,'Frequency / kHz','FontSize',pp.fsize+1,'Interpreter','tex','FontName',pp.fname)
ylabel(tl,'Level / dB','FontSize',pp.fsize+1,'Interpreter','tex','FontName',pp.fname)
c = lines(3);
c(4,:) = [0.5 0.5 0.5];

locs.Violin = [8 11 16];
locs.VocalAlto = [8 10  13];
locs.ClarinetBb = [7 10 13];
locs.Tuba = [2 5 8];

fNat = 20:12000;

c4 = [1 1 1];

for iInstr = 1:nInstr
    instrName = instruments{iInstr};
    f = muspitch2freq(pitches(locs.(instrName)));

    for ii = 1:3
        nexttile(iInstr+4*(ii-1))
        hold on

        interpSE = interp1(cf,trueSE.(instrName)(locs.(instrName)(1),:),fNat,'spline');
        harmonics = round((1:60).*f(ii));
        amps = zeros(1,numel(fNat));
        amps(ismember(fNat,harmonics)) = 1;
        amps = interpSE.*amps;
        idx = find(amps);

        patch([20:idx(1)+19 fliplr(20:idx(1)+19)],[interpSE(1:idx(1)) -90.*ones(1,idx(1))],[.9 .9 .9],...
            'EdgeColor',[.9 .9 .9])

        for jj = 1:numel(idx)
            ph = plot([harmonics(jj) harmonics(jj)],[-90 amps(idx(jj))],...
                'Color',[c(1,:) 0.5],...
                'LineWidth',pp.linewidth);

        end

        if ii == 2
            pc = plot(cf,trueSE.(instrName)(locs.(instrName)(2),:),...
                'Color',c(2,:),...
                'LineWidth',pp.linewidth,...
                'LineStyle','-');
        end

        if ii == 3
            pc = plot(cf,trueSE.(instrName)(locs.(instrName)(3),:),...
                'Color',c(3,:),...
                'LineWidth',pp.linewidth,...
                'LineStyle','-');
        end

        pi = plot(cf,trueSE.(instrName)(locs.(instrName)(1),:),...
            'Color',[c(1,:) c4(ii)],...
            'LineWidth',pp.linewidth+0.5);
        hold off

        %     if ii == 1
        %         legh = legend(pi,'LRSE',...
        %             'Location','southwest');
        %     elseif ii == 2
        %         legh = legend(pc,'MRSE',...
        %             'Location','southwest');
        %     else
        %         legh = legend(pc,'HRSE',...
        %             'Location','southwest');
        %     end
        %     legh.ItemTokenSize = [10 1 0];
        %     legh.Position(1) = legh.Position(1) - 0.0000001;
        %     legh.Position(2) = legh.Position(2) - 0.0000001;

        %     if ii == 1
        %         title('A','Position',[5 4 0],'HorizontalAlignment','center')
        %     end

        set(gca,'FontSize',pp.fsize,...
            'LineWidth',pp.linewidth,...
            'Layer','top',...
            'XScale','log','FontName',pp.fname)

        if ii == 1
            title(instrumentNames{iInstr})
        end

        xlim([20 12000])
        xticks([100 1000 10000])
        xticklabels({'0.1','1','10'})

%         if ~ismember(iInstr+4*(ii-1),[9 10 11 12])
%             xticklabels({})
%         end

        ylim([-90 15])
        yticks([-90 -60 -30 0])

%         if ~ismember(iInstr+4*(ii-1),[1 5 9])
%             yticklabels({})
%         end

        box on

    end

    % locs = [10 13 16];
    % f = muspitch2freq(pitches(locs));
    % fNat = 20:12000;
    %
    % for ii = 1:3
    %     nexttile
    %     hold on
    %
    %     interpSE = interp1(cf,trueSE.ClarinetBb(10,:),fNat,'spline');
    %     harmonics = round((1:60).*f(ii));
    %     amps = zeros(1,numel(fNat));
    %     amps(ismember(fNat,harmonics)) = 1;
    %     amps = interpSE.*amps;
    %     idx = find(amps);
    %
    %     for jj = 1:numel(idx)
    %         ph = plot([harmonics(jj) harmonics(jj)],[-90 amps(idx(jj))],...
    %             'Color',[c(2,:) 0.5],...
    %             'LineWidth',pp.linewidth);
    %
    %     end
    %
    %     if ii == 2
    %         pc = plot(cf,trueSE.ClarinetBb(13,:),...
    %             'Color',c(3,:),...
    %             'LineWidth',pp.linewidth,...
    %             'LineStyle','-');
    %     end
    %
    %     if ii == 3
    %         pc = plot(cf,trueSE.ClarinetBb(16,:),...
    %             'Color',c(4,:),...
    %             'LineWidth',pp.linewidth,...
    %             'LineStyle','-');
    %     end
    %
    %     pi = plot(cf,trueSE.ClarinetBb(10,:),...
    %         'Color',[c(2,:) c4(ii)],...
    %         'LineWidth',pp.linewidth);
    %     hold off
    %
    %     if ii == 1
    %         legh = legend([pi ph],{'MRA-SE','H(MRA)'},...
    %             'Location','best');
    %     elseif ii == 2
    %         legh = legend([pc ph],{'HRA-SE','H(HRA)'},...
    %             'Location','best');
    %         legh.Position(1) = legh.Position(1) - 0.0387;
    %         legh.Position(2) = legh.Position(2) - 0.188;
    %     else
    %         legh = legend([pc ph],{'F#6-SE','H(F#6)'},...
    %             'Location','best');
    %         legh.Position(1) = legh.Position(1) - 0.105;
    %         legh.Position(2) = legh.Position(2) - 0.188;
    %     end
    %     legh.ItemTokenSize = [10 1 0];
    % %     legh.Position(1) = legh.Position(1) - 0.05;
    % %     legh.Position(2) = legh.Position(2) - 0.05;
    %
    %     if ii == 1
    %         title('B','Position',[5 4 0],'HorizontalAlignment','center')
    %     end
    %
    %     set(gca,'FontSize',pp.fsize,...
    %         'LineWidth',pp.linewidth,...
    %         'Layer','top',...
    %         'XScale','log')
    %
    %     xlim([20 12000])
    %     xticks([100 1000 10000])
    %     xticklabels({'0.1','1','10'})
    %
    %     ylim([-90 10])
    %     yticks([-90 -60 -30 0])
    %
    %     box on
    %
    % end
end

if isfield(pp, 'figwidth')
    if ~isempty(pp.figwidth)
        set(figh1, 'PaperPositionMode', 'manual');
        set(figh1, 'PaperUnits', 'centimeters');
        set(figh1, 'PaperPosition', [0 0 pp.figwidth pp.figheight-2]);
        set(figh1, 'Toolbar', 'none')
        set(figh1, 'Menubar', 'none')
        set(figh1, 'PaperSize', [pp.figwidth pp.figheight-2])
    end
end

if pp.print
    print(figh1, [fig_folder filesep 'SE_F0_shift.pdf'], '-dpdf');
end
