clear
close all

% define path of data
allFiles = dir('/home/simon/ma/code/soundQuality/data/all/');

% get all file names
allNames = {allFiles(~[allFiles.isdir]).name};

% get number of participants
nPart = numel(allNames);

% define instruments
instruments = {'Violin','VocalAlto','ClarinetBb','Tuba'};
nInstr = numel(instruments);

% define experiments
experiments = {'quality','plausib'};
nExp = numel(experiments);

% define conditions
conditions = {'low','mid','hig','all'};
nCond = numel(conditions);

% define maximum of stimuli per condition and experiment
nNote = 19;

% preallocate data matrices for all instruments
for iInstr = 1:nInstr
    instrName = instruments{iInstr};
    expData_raw.quality.(instrName) = zeros(nPart,nNote,nCond);
    expData_raw.plausib.(instrName) = zeros(nPart,nNote,nCond);
    expData_norm.quality.(instrName) = zeros(nPart,nNote,nCond);
    expData_norm.plausib.(instrName) = zeros(nPart,nNote,nCond);
    rowNum.quality.(instrName) = zeros(nPart,nNote,nCond);
    rowNum.plausib.(instrName) = zeros(nPart,nNote,nCond);
end

% define pitch names
pitches = {'FIS1','AIS1','D2',...
    'FIS2','AIS2','D3',...
    'FIS3','AIS3','D4',...
    'FIS4','AIS4','D5',...
    'FIS5','AIS5','D6',...
    'FIS6','AIS6','D7',...
    'FIS7'};
nPitch = numel(pitches);

% read table
for iPart = 1:nPart
    filename = strcat('/home/simon/ma/code/soundQuality/data/all/',allNames{iPart});
    T = sortrows(readtable(filename));

    for iExp = 1:nExp
        expName = experiments{iExp};

        rawData = [];

        for iInstr = 1:nInstr
            instrName = instruments{iInstr};
    
            for iCond = 1:nCond
                condName = conditions{iCond};
        
                for iPitch = 1:nPitch
                    pitchName = pitches{iPitch};
    
                    stimName = strcat(instrName,'_',condName,'_',pitchName);
                    idx = find(strcmp(T.stim1,stimName));
    
                    if ~isempty(idx)
                        if numel(idx) == 2
                            expData_raw.(expName).(instrName)(iPart,iPitch,iCond) = T.response(idx(iExp));
                            rowNum.(expName).(instrName)(iPart,iPitch,iCond) = T.rowNo(idx(iExp));
                        elseif numel(idx) == 1
                            expData_raw.quality.(instrName)(iPart,iPitch,iCond) = T.response(idx);
                            expData_raw.plausib.(instrName)(iPart,iPitch,iCond) = nan;
                            rowNum.quality.(instrName)(iPart,iPitch,iCond) = T.rowNo(idx);
                            rowNum.plausib.(instrName)(iPart,iPitch,iCond) = nan;
                        end
                    else
                        expData_raw.(expName).(instrName)(iPart,iPitch,iCond) = nan;
                        rowNum.(expName).(instrName)(iPart,iPitch,iCond) = nan;
                    end
    
                end
    
            end

%             expData.(expName).(instrName)(iPart,:,:) = rescale(expData.(expName).(instrName)(iPart,:,:),1,6);
            rawData(:,:,iInstr) = squeeze(expData_raw.(expName).(instrName)(iPart,:,:))';
            
        end

        normData = rescale(rawData,1,6);

        for iInstr = 1:nInstr

            instrName = instruments{iInstr};
            expData_norm.(expName).(instrName)(iPart,:,:) = normData(:,:,iInstr)';

        end

    end

end

save('/home/simon/ma/code/soundQuality/data/expData_raw.mat','expData_raw')
save('/home/simon/ma/code/soundQuality/data/expData_norm.mat','expData_norm')