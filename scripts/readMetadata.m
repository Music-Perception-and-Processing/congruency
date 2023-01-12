clear
close all

% male = 1, female = 2, divers = 3

%% Pilot study
% 
% nPart = 7;
% partAge = [27;23;21;26;22;29;25];
% partGender = [2;1;2;1;1;1;1];


%% Main study

% define path of data
allFiles = dir('data/responses/');

% get all file names
allNames = {allFiles(~[allFiles.isdir]).name};

% get number of participants
nPart = numel(allNames);

% read table
for iPart = 1:nPart
    filename = strcat('data/responses/',allNames{iPart});
    opts = detectImportOptions(filename);
    T = sortrows(readtable(filename));

    metadata(iPart,1) = T.response(find(strcmp(T.condition1,'age')));
    metadata(iPart,2) = T.responseCode(find(strcmp(T.condition1,'gender')));
    metadata(iPart,3:17) = T.responseCode(find(strcmp(T.condition1,'msiPerception11')):find(strcmp(T.condition1,'msiTraining37')));

end

age = metadata(:,1);
meanAge = mean(age,'omitnan')
stdAge = std(age,'omitnan')

nMale = sum(metadata(:,2)==1)
nFemale = sum(metadata(:,2)==2)
nDivers = sum(metadata(:,2)==3)

% save('/home/simon/ma/code/soundQuality/data/expData.mat','expData')