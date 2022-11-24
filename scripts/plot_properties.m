% plot properties
pp.visible = 'off';  % 'on'/'off', if figures only need to be saved
pp.fsize = 10;
pp.figwidth = 15;   % width in cm
pp.figheight = 10;  % height in cm
pp.linewidth = 1;
pp.markersize = 5;
pp.scatterMsize = 10;
pp.fname = 'Helvetica';

fig_folder = './plots';

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