clear

% Load instrument spectra
spectra = load('data/spectrum_adj.mat');

% Get instrument info
load('data/analysis_v2.mat');
classNames = fieldnames(an.fam);
nClasses = numel(classNames)-1;

% Fix Wagnertuba note name
an.range.Tubas{2,1} = 'AIS1';

% Get global note info
load('data/NoteData.mat')
nPitches = numel(NoteData.MP);

%% Make restructured spectra variable
% Initialize variables
S = nan(nPitches,22050,50);
cInstr = 1; % instrument counter
failedInstruments = {};

% Loop over classes
for iClass = 1:nClasses
    % Get class name
    className = classNames{iClass};

    % Get instruments in class
    instrNames = an.fam.(className);
    nInstr = numel(instrNames);

    % Loop over instruments in class
    for iInstr = 1:nInstr
        % Get instrument name and range
        instrName = instrNames{iInstr};
        instrRange = an.range.(className)(iInstr,:);
    
        % Get instrument range
        rangeMin = find(strcmp(NoteData.MP,instrRange{1}));
        rangeMax = find(strcmp(NoteData.MP,instrRange{2}));
        rangeIdx = rangeMin:rangeMax;

        if strcmp('Viola',instrName)
            rangeIdx = rangeIdx(rangeIdx~=71);
        end

            

        % if any(strcmp({'Viola','Wagnertuba'},instrName))
        %     nNotes = size(spectra.spectrum.(instrName).x.spectrum,1);
        %     for iNote = 1:nNotes
        %         noteName = NoteData.MP{rangeIdx(iNote)};
        %         buff = muspitch2freq(noteName)*0.9;
        %         [~,loc] = findpeaks(spectra.spectrum.(instrName).x.spectrum(iNote,:),...
        %             'MinPeakHeight',0.1,'MinPeakDistance',buff);
        %         if size(loc) == 1
        %             f0.(instrName){iNote} = freq2muspitch(loc);
        %         else
        %             f0.(instrName){iNote} = freq2muspitch(median(diff(loc)));
        %         end
        % 
        %     end
        % 
        % end

        % Write spectra to variable
        try
            S(rangeIdx,:,cInstr) = spectra.spectrum.(instrName).x.spectrum;
        catch exception
            failedInstruments = [failedInstruments;instrName];
            continue
        end

        % Update instrument counter
        cInstr = cInstr + 1;

    end

end

%% Compute spectral centroids
f = 1:22050;
sc = [];
for iInstr = 1:size(S,3)
    SC(:,iInstr) = sum(f.*S(:,:,iInstr),2)./sum(S(:,:,iInstr),2);
    sc = [sc;SC(~isnan(SC(:,iInstr)),iInstr)];
    nNotes(iInstr,1) = sum(~isnan(SC(:,iInstr)));
end

%% Compute spectral spread
ssp = [];
for iInstr = 1:size(S,3)
    SSP(:,iInstr) = sqrt(sum((f - SC(:,iInstr)).^2.*S(:,:,iInstr),2)./sum(S(:,:,iInstr),2));
    ssp = [ssp;SSP(~isnan(SSP(:,iInstr)),iInstr)];
end

%% Compute spectral skewness
ssk = [];
for iInstr = 1:size(S,3)
    SSK(:,iInstr) = sum((f - SC(:,iInstr)).^3.*S(:,:,iInstr),2)./(SSP(:,iInstr).^3.*sum(S(:,:,iInstr),2));
    ssk = [ssk;SSK(~isnan(SSK(:,iInstr)),iInstr)];
end

%% Compute spectral Kurtosis
sku = [];
for iInstr = 1:size(S,3)
    SKU(:,iInstr) = sum((f - SC(:,iInstr)).^4.*S(:,:,iInstr),2)./(SSP(:,iInstr).^4.*sum(S(:,:,iInstr),2));
    sku = [sku;SKU(~isnan(SKU(:,iInstr)),iInstr)];
end

%% Compute correlations
load('data/pca_data_v5.mat');

for iDim = 1:12
    [SCcorr(iDim,1),SCcorr(iDim,2)] = corr(pcaData.score(:,iDim),sc,'type','Spearman');
    [F0corr(iDim,1),F0corr(iDim,2)] = corr(pcaData.score(:,iDim),pcaData.F0,'type','Spearman');
    [SSPcorr(iDim,1),SSPcorr(iDim,2)] = corr(pcaData.score(:,iDim),ssp,'type','Spearman');
    [SSKcorr(iDim,1),SSKcorr(iDim,2)] = corr(pcaData.score(:,iDim),ssk,'type','Spearman');
    [SKUcorr(iDim,1),SKUcorr(iDim,2)] = corr(pcaData.score(:,iDim),sku,'type','Spearman');
end

%%
writetable(sc_correlation,'/Users/simon/congruency/data/sc_correlation.csv', 'WriteRowNames', false,...
    'Delimiter', ',');