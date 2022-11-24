clear
close all
load('/home/simon/ma/code/data/spectrum.mat')
load('/home/simon/ma/code/data/analysis_v2.mat')
load('/home/simon/ma/code/data/erbCC.mat')
load('/home/simon/ma/code/data/NoteData.mat')
addpath('/home/simon/ma/code/functions/')

% set audiowrite flag
sflag = 0;

% setup instruments
instruments = {'Violin','VocalAlto','ClarinetBb','Tuba'};
nInstr = numel(instruments);

regBounds.Violin = {'D4','E5'};
regBounds.VocalAlto = {'D4','D5'};
regBounds.ClarinetBb = {'D4','A4'};
regBounds.Tuba = {'F2','F3'};

% get spectra
for iInstr = 1:nInstr

    instrName = instruments{iInstr};
    S.(instrName) = spectrum.(instrName).x.spectrum;
    rangeF0.(instrName) = erbCC.(instrName).x.F0;
    rangeMP.(instrName) = freq2muspitch(rangeF0.(instrName));

end

% define F0s
minF0 = 'FIS1';
maxF0 = 'FIS7';
allF0 = NoteData.F0(10:end-7);
allMP = NoteData.MP(10:end-7);
nNote = numel(allF0);

% define notes in one octave
intF0 = {'FIS','AIS','D'};
offsetOct = [0 0 1];
nInt = numel(intF0);

% create F0 vector
nOct = 6;
vF0 = zeros(nOct*nInt+1,1);
nF0 = numel(vF0);
iF0 = 1;
for iOct = 1:nOct

    for iInt = 1:nInt

        vMP{iF0} = strcat(intF0{iInt},num2str(iOct+offsetOct(iInt)));
        vF0(iF0) = muspitch2freq(vMP(iF0));
        iF0 = iF0 + 1;

    end

end
vMP{nOct*nInt+1} = maxF0;
vF0(end) = muspitch2freq(maxF0);

% plot note overlap
figure, hold on
for iInstr = 1:nInstr

    instrName = instruments{iInstr};
    plot(rangeF0.(instrName),iInstr*ones(1,numel(rangeF0.(instrName))), ...
        'LineWidth',2)
    scatter(muspitch2freq(regBounds.(instrName)),iInstr*ones(1,2),80,'red', ...
        'Marker','x','LineWidth',2)
    scatter(vF0,iInstr*ones(1,numel(vF0)),20,'black')

    for ii = 1:numel(vF0)
        text(vF0(ii),iInstr+0.25,vMP{ii}, ...
            'FontSize',6,'HorizontalAlignment','center')
    end

end
hold off
ylim([0 nInstr+1])
set(gca,'XScale','log')

sampleLocs.Violin = [8 11 16];
sampleLocs.VocalAlto = [8 10 13];
sampleLocs.ClarinetBb = [7 10 13];
sampleLocs.Tuba = [2 5 8];

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

% prepare instruments' spectral envelopes
for iInstr = 1:nInstr

    instrName = instruments{iInstr};

    for iNote = 1:nF0

        idx = find(strcmp(rangeMP.(instrName),vMP{iNote}));

        if idx
            
            X.(instrName)(iNote,:) = S.(instrName)(idx,:);
            Xfilt = fb*X.(instrName)(iNote,:)';
%             Xfilt = Xfilt./rms(Xfilt);

            % max of filter bank
            Xmax = zeros(1,numBands);
            for iBand = 1:numBands
                Xmax(iBand) = max(X.(instrName)(iNote,floor(cf(iBand)-bw(iBand)/2):ceil(cf(iBand)+bw(iBand)/2)));
            end

            xx = 20.*dct(log10(Xfilt+eps),[],1);
%             xx = 20.*dct(log10(X.(instrName)(iNote,:)'+eps),[],1);
            xx(nCC+1:end) = 0;
            XX = idct(xx,[],1);

            theta = 2;
            trueSE.(instrName)(iNote,:) = trueEnvelope(Xmax,XX,nCC,theta);
%             trueSE.(instrName)(iNote,:) = trueEnvelope(X.(instrName)(iNote,:),XX,nCC,theta);

        else

            X.(instrName)(iNote,:) = zeros(1,fs/2);
            trueSE.(instrName)(iNote,:) = zeros(1,numBands);

        end

    end

    figure
    tl = tiledlayout(4,5,'TileSpacing','compact','Padding','compact');

    for iNote = 1:nF0

        nexttile
        hold on
        plot(f,20.*log10(X.(instrName)(iNote,:)))
        plot(cf,trueSE.(instrName)(iNote,:))
        hold off
        xlim([20 12000])
        ylim([-90 10])
        set(gca,'XScale','log')

    end

    title(tl,instrName)
    xlabel(tl,'Frequency / Hz')
    ylabel(tl,'Level / dB')
end

% figure,hold on
% plot(f,20.*log10(X.ClarinetBb(10,:)))
% plot(cf,trueSE.ClarinetBb(10,:))
% hold off
% xlim([20 12000])
% ylim([-90 0])
% set(gca,'XScale','log')

%%
% define synthesis parameters
t = 0:1/fs:0.5;
cosWin = tukeywin(fs/2+1,0.3)';
nHar = 32;
phase = 0;

conditions = {'low','mid','hig','all'};
nCond = numel(conditions);

allSig = [];

for iInstr = 1:nInstr

    instrName = instruments{iInstr};

    for iCond = 1:nCond

        condName = conditions{iCond};

        for iNote = 1:nF0

            F0 = vF0(iNote);
            noteName = vMP{iNote};

%             if strcmp(noteName,'AIS4')
%                 disp('stop')
%             end

            if strcmp('all',condName)
                SE = trueSE.(instrName)(iNote,:);
            else
                SE = trueSE.(instrName)(sampleLocs.(instrName)(iCond),:);
            end

%             if strcmp(instrName,'ClarinetBb') && iCond == 4 && iNote == 14
% 
%                 disp('break')
% 
%             end

            % create harmonics and corresponding amplitudes
            [h,a] = generateHarmonicsNoLim(SE,cf,20:12000,F0,nHar);

            % additive synthesis
            sig.(instrName)(iNote,:,iCond) = sum(a'.*sin(2.*pi.*h'.*t + phase)).*cosWin;

            % RMS normalization
            sig.(instrName)(iNote,:,iCond) = sig.(instrName)(iNote,:,iCond)./rms(sig.(instrName)(iNote,:,iCond));

        end

    end

end

% find max amplitude across instruments and sounds
for iInstr = 1:nInstr

    instrName = instruments{iInstr};
    A_max(iInstr) = max(abs(sig.(instrName)),[],'all');

end

maxValue = max(A_max)

for iInstr = 1:nInstr

    % normalize signals to avoid clipping
    instrName = instruments{iInstr};
    sig_norm = sig.(instrName)./maxValue;

    for iCond = 1:nCond
        
        condName = conditions{iCond};

        for iNote = 1:nF0

            F0 = vF0(iNote);
            noteName = vMP{iNote};

            foldName = strcat(instrName,'/');
            stimName = strcat(instrName,'_',condName,'_',noteName,'.wav');
            fileName = strcat('/home/simon/ma/code/soundQuality/sounds/new_sounds/',foldName,stimName);

            if sflag
                audiowrite(fileName,sig_norm(iNote,:,iCond),fs);
            end

        end

        [r,c] = size(squeeze(sig_norm(:,:,iCond)));
        testSig = reshape(squeeze(sig_norm(:,:,iCond))',1,r*c);
        stimName = strcat(instrName,'_',condName,'_test.wav');
        fileName = strcat('/home/simon/ma/code/soundQuality/sounds/new_sounds/',foldName,stimName);

        if sflag
            audiowrite(fileName,testSig,fs);
        end

        if iCond == 4
            allSig =  [allSig; squeeze(sig_norm(:,:,4))];
        end

    end

end

[r,c] = size(allSig);
randIdx = randperm(r);
testSig = reshape(allSig(randIdx,:)',1,r*c);
stimName = 'appetizer.wav';
fileName = strcat('/home/simon/ma/code/soundQuality/sounds/new_sounds/',stimName);

if sflag
    audiowrite(fileName,testSig,fs);
end

save('/home/simon/ma/code/soundQuality/data/trueSE.mat','trueSE')