% plot properties
pp.visible = 'off'; % 'on'/'off', if figures only need to be saved
pp.print = 1; % 0 / 1 if figures need to be printed
pp.fsize = 10;
% pp.fsize = 9;
pp.figwidth = 15;   % width in cm
pp.figheight = 10;  % height in cm
pp.linewidth = 1;
pp.markersize = 5;
pp.scatterMsize = 10;
pp.fname = 'Times New Roman';
% pp.fname = 'Arial';

if strcmp(pp.fname,'Arial')
    fig_folder = 'figures/arial/';
elseif strcmp(pp.fname,'Times New Roman')
    fig_folder = 'figures/tnr/';
end

% define instruments
instruments = {'Violin','VocalAlto','ClarinetBb','Tuba'};
nInstr = numel(instruments);

instrumentNames = {'Violin','Alto voice','Clarinet','Tuba'};

% define experiments
experiments = {'quality','plausib'};
nExp = numel(experiments);

% define conditions
conditions = {'low','mid','hig','all'};
nCond = numel(conditions);

conditionNames = {'low','mid','high','all'};

% instrument registers
r.Violin = [8 11 16];
r.VocalAlto = [8 10 13];
r.ClarinetBb = [7 10 13];
r.Tuba = [2 5 8];