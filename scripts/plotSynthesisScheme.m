clear
close all
clc

addpath("functions/");
fig_folder = 'figures/';

%% Color palette
hex_colors  = {'#0072BD', '#79ABDA', ...
               '#D95319', '#EF9C7E', ...
               '#EDB120', '#F5D18E'};
c1 = hex2rgb(hex_colors([1 3 5])); 
c2 = hex2rgb(hex_colors([2 4 6]));

%% Load clarinet example sounds

% Select audio folder
projectFolder = pwd;
audioFolder = strcat(projectFolder,filesep,'sounds',filesep,'clarinet_example');

% Get all audio file names
audioFiles = dir(audioFolder);
fileNames = {audioFiles(~[audioFiles.isdir]).name};
fileNames = fileNames([1 3]);

% Read audio
for iFile = 1:numel(fileNames)
    [x, fs] = audioread(strcat(audioFolder,filesep,fileNames{iFile}));

    % Select only the left channel
    s(:,iFile) = x(:,2);

end


s_plot = s;
if numel(fileNames) == 2
    % Adjust wave 1
    s_plot(:,1) = 0;
    s_plot(1:fs, 1) = s(0.2*fs:fs+0.2*fs-1, 1);

    s_plot(:,2) = 0;
    s_plot(0.2*fs:fs+0.2*fs-1, 2) = s(1:fs, 2);
elseif numel(fileNames) == 3
    % Adjust wave 3
    s_plot(:,3) = 0;
    s_plot(0.2*fs:fs+0.2*fs-1, 3) = s(1:fs, 3);
end

% Plot wave form
[figh, pp] = plotHelper;

T = size(s_plot,1)/2;
p = plot(s_plot(1:T,:)./max(s_plot), ...
    "LineWidth", 0.5);

for ii = 1:size(p, 1)
    p(ii).Color = c1(ii,:);
end

set(gca, "FontSize", pp.fsize, ...
    "LineWidth", pp.linewidth, ...
    "Layer", "top");

% x-axis
xlim([0 fs])
xticks([0 0.25*fs 0.5*fs 0.75*fs fs])
xticklabels({'0', '0.25', '0.5', '0.75', '1'})
xlabel('Time / s', ...
    "FontSize", pp.fsize+1, ...
    "Interpreter", "tex", ...
    "FontName", pp.fname);

% y-axis
ylim([-1.1 1.1])
ylabel('Amplitude', ...
    "FontSize", pp.fsize, ...
    "Interpreter", "tex", ...
    "FontName", pp.fname);

title('Waveform', "FontSize", pp.fsize)

chH = get(gca,'Children');
set(gca,'Children',[chH(end);chH(1:end-1)]);

fig_name = 'clarinet_waveform';
printHelper(figh, pp, fig_folder, fig_name)

%% Compute spectra

% Take absolute magnitude spectrum
S = abs(fft(s, fs));
% Take one-sided spectrum
S = S(1:end/2, :);
% Normalize
S = S./max(S);
% Set dynamic range
S(20*log10(S) < -90) = 0;

% Plot spectrum
[figh, pp] = plotHelper;

semilogx(20*log10(S), ...
    "LineWidth", 1);

set(gca, "FontSize", pp.fsize, ...
    "LineWidth", pp.linewidth, ...
    "Layer", "top", ...
    "XScale", "log");

% x-axis
xlim([20 12000])
xticks([100, 1000, 10000])
xticklabels({'100', '1k', '10k'})
xlabel('Frequency / Hz', ...
    "FontSize", pp.fsize+1, ...
    "Interpreter", "tex", ...
    "FontName", pp.fname);

% y-axis
ylim([-80 1]);
ylabel('Level / dB', ...
    "FontSize", pp.fsize, ...
    "Interpreter", "tex", ...
    "FontName", pp.fname);

title('Spectrum', "FontSize", pp.fsize)

chH = get(gca,'Children');
set(gca,'Children',[chH(end);chH(1:end-1)]);

fig_name = 'clarinet_spectrum';
printHelper(figh, pp, fig_folder, fig_name)

%% Manipulation of spectrum
S_m = S;
[a, idx] = max(S_m);
for ii = 1:size(S_m,2)
    S_m(1:idx(ii),ii) = a(ii);
end

% Plot spectrum
[figh, pp] = plotHelper;

semilogx(20*log10(S_m), ...
    "LineWidth", pp.linewidth);

set(gca, "FontSize", pp.fsize, ...
    "LineWidth", pp.linewidth, ...
    "Layer", "top", ...
    "XScale", "log");

% x-axis
xlim([20 12000])
xticks([100, 1000, 10000])
xticklabels({'100', '1k', '10k'})
xlabel('Frequency / Hz', ...
    "FontSize", pp.fsize+1, ...
    "Interpreter", "tex", ...
    "FontName", pp.fname);

% y-axis
ylim([-80 1]);
ylabel('Level / dB', ...
    "FontSize", pp.fsize, ...
    "Interpreter", "tex", ...
    "FontName", pp.fname);

title('Manipulated spectrum', "FontSize", pp.fsize)

chH = get(gca,'Children');
set(gca,'Children',[chH(end);chH(1:end-1)]);

fig_name = 'clarinet_spectrum_manipulated';
printHelper(figh, pp, fig_folder, fig_name)

%% Plot both spectra in one figure
[figh, pp] = plotHelper;

p = semilogx(20*log10(S_m), ...
    "LineWidth", pp.linewidth);

for ii = 1:size(p, 1)
    p(ii).Color = c1(ii,:);
end

for ii = 1:size(S, 2)
    [a, idx] = max(S(:,ii));
    low_energy_spectrum = S(1:idx,ii);

    pl = semilogx(1:idx, 20*log10(low_energy_spectrum), ...
        "LineWidth", pp.linewidth, ...
        "Color", c2(ii,:));

end

set(gca, "FontSize", pp.fsize, ...
    "LineWidth", pp.linewidth, ...
    "Layer", "top", ...
    "XScale", "log");

% x-axis
xlim([20 12000])
xticks([100, 1000, 10000])
xticklabels({'100', '1k', '10k'})
xlabel('Frequency / Hz', ...
    "FontSize", pp.fsize+1, ...
    "Interpreter", "tex", ...
    "FontName", pp.fname);

% y-axis
ylim([-80 10]);
ylabel('Level / dB', ...
    "FontSize", pp.fsize, ...
    "Interpreter", "tex", ...
    "FontName", pp.fname);

title('Spectrum (manipulated)', "FontSize", pp.fsize)

chH = get(gca,'Children');
set(gca,'Children',[chH(end);chH(1:end-1)]);

fig_name = 'clarinet_spectrum_both';
printHelper(figh, pp, fig_folder, fig_name)

%% Compute erb magnitudes

% Create auditory filter bank
[fb, cf, bw] = designAuditoryFilterBank(fs, ...
    "FrequencyScale", "erb", ...
    "FFTLength", fs-1, ...
    "NumBands", 128, ...
    "FrequencyRange", [20 12000], ...
    "Normalization", "bandwidth");

% Compute ERB-magnitude spectrum
S_erb = fb*S_m;
S_erb_og = fb*S;

% Alternative computation
for ii = 1:size(S_m,2)
    S_erb_alt(:,ii) = computeERBfilterBank(S_m(:,ii), fb);
end

% Normalize ERB spectrum
S_erb = S_erb ./ rms(S_erb);
% S_erb_og = S_erb_og ./ rms(S_erb_og);
S_erb_alt = S_erb_alt ./ rms(S_erb_alt);

[figh, pp] = plotHelper;

for ii = 1:size(S_erb, 2)
    plot(cf, 20.*log10(S_erb(:,ii)), ...
        "LineWidth", pp.linewidth, ...
        "Color", c1(ii,:));
    pl = plot(cf, 20.*log10(S_erb_og(:,ii)), ...
        "LineWidth", pp.linewidth, ...
        "Color", c2(ii,:));
end

set(gca, "FontSize", pp.fsize, ...
    "LineWidth", pp.linewidth, ...
    "Layer", "top", ...
    "XScale", "log");

% x-axis
xlim([20 12000])
xticks([100, 1000, 10000])
xticklabels({'100', '1k', '10k'})
xlabel('Frequency / Hz', ...
    "FontSize", pp.fsize+1, ...
    "Interpreter", "tex", ...
    "FontName", pp.fname);

% y-axis
ylim([-80 20]);
ylabel('Level / dB', ...
    "FontSize", pp.fsize, ...
    "Interpreter", "tex", ...
    "FontName", pp.fname);

title('ERB-magnitude spectrum', "FontSize", pp.fsize)

chH = get(gca,'Children');
set(gca,'Children',[chH(end);chH(1:end-1)]);

fig_name = 'clarinet_erbMagn_spectrum_both';
printHelper(figh, pp, fig_folder, fig_name)

%% Compute cepstrum

% Compute cepstrum
ss = dct(log10(S_erb + eps), [], 1);
ss_og = dct(log10(S_erb_og + eps), [], 1);

% Smooth cepstrum
ss_smooth = ss;
ss_smooth(14:end, :) = 0;

% Plot cepstral coefficients
[figh, pp] = plotHelper;

tl = tiledlayout(2, 1, "TileSpacing", "none", "Padding" ,"loose");

nexttile
b1 = bar(ss_smooth(1:13, :));

for ii = 1:size(b1, 2)
    b1(ii).FaceColor = c1(ii,:);
    b1(ii).EdgeColor = c1(ii,:);
    b1(ii).BaseLine.Visible = "off";
end

hold on

plot([1.5 1.5], [-40 30], ...
    "Color", [0.4 0.4 0.4], ...
    "LineWidth", pp.linewidth, ...
    "LineStyle", "--");

plot([0 14], [0 0], ...
    "Color", [0 0 0], ...
    "LineWidth", 0.5);

% x-axis
xlim([0 14])
xticks(1:13)
xticklabels({'', '', '', '', '', '', '', '', '', '', '', ''})
xtickangle(0)

% y-axis
ylim([-35 20])
ylabel(tl, 'Level', ...
    "FontSize", pp.fsize, ...
    "Interpreter", "tex", ...
    "FontName", pp.fname);

set(gca, "FontSize", pp.fsize, ...
    "LineWidth", pp.linewidth, ...
    "Layer", "top");

title('ERB-frequency cepstral coefficients', "FontSize", pp.fsize)

nexttile;
b2 = bar(ss_og(1:13, :));

for ii = 1:size(b2, 2)
    b2(ii).FaceColor = c2(ii,:);
    b2(ii).EdgeColor = c2(ii,:);
    b2(ii).BaseLine.Visible = "off";
end

hold on

plot([0 14], [0 0], ...
    "Color", [0 0 0], ...
    "LineWidth", 0.5);

% x-axis
xlim([0 14])
xticks(1:13)
xticklabels({'1' ,'2', '', '4', '', '6', '', '8', '', '10', '', '12', ''})
xtickangle(0)
xlabel(tl, 'Cepstral coefficient number', ...
    "FontSize", pp.fsize+1, ...
    "Interpreter", "tex", ...
    "FontName", pp.fname);

% y-axis
ylim([-40 15])
ylabel(tl, 'Level', ...
    "FontSize", pp.fsize, ...
    "Interpreter", "tex", ...
    "FontName", pp.fname);

set(gca, "FontSize", pp.fsize, ...
    "LineWidth", pp.linewidth, ...
    "Layer", "top");

fig_name = 'clarinet_cepstral_coefficients';
printHelper(figh, pp, fig_folder, fig_name)

%% Compute spectral envelope
SS = 20.*idct(ss_smooth, [], 1);

[figh, pp] = plotHelper;

plot(cf, SS, ...
    "LineWidth", pp.linewidth);

set(gca, "FontSize", pp.fsize, ...
    "LineWidth", pp.linewidth, ...
    "Layer", "top", ...
    "XScale", "log");

% x-axis
xlim([20 12000])
xticks([100 1000 10000])
xticklabels({'100', '1k', '10k'})
xlabel('Frequency / Hz', ...
    "FontSize", pp.fsize+1, ...
    "Interpreter", "tex", ...
    "FontName", pp.fname);

% y-axis
ylim([-80 20])
ylabel('Level / dB', ...
    "FontSize", pp.fsize, ...
    "Interpreter", "tex", ...
    "FontName", pp.fname);

title('Spectral envelope', "FontSize", pp.fsize)

chH = get(gca,'Children');
set(gca,'Children',[chH(end);chH(1:end-1)]);

fig_name = 'clarinet_spectral_envelope';
printHelper(figh, pp, fig_folder, fig_name)

%% True envelope computation
% Compute cepstrum
ss_og = dct(log10(S_erb_og + eps), [], 1);
% Smooth cepstrum
ss_smooth_og = ss_og;
ss_smooth_og(14:end, :) = 0;

SS_og = 20.*idct(ss_smooth_og, [], 1);

% Get max in each filter
S_max = zeros(128, size(S, 2));
for ii = 1:128
    S_max(ii,:) = max(S(floor(cf(ii) - bw(ii)/2):ceil(cf(ii) + bw(ii)/2), :));
end

theta = 2;
nCC = 13;
for ii = 1:size(S, 2)
    [SE(:,ii), CC(:,ii)] = trueEnvelope2(S_max(:,ii), SS_og(:,ii), nCC, theta);
end

%% Move into component space

% Load pca data
load("data/pca_data_v5.mat");

% Project clarinet examples on component space
score = (ss_smooth(2:13,:)' - pcaData.mu) * pcaData.coeff;

% Old scores
clarinet_notes = {'FIS3', 'FIS4', 'FIS5'};
clarinet_notes = {'FIS3', 'FIS5'};
clarinet_idx = strcmp(pcaData.instrLabel, 'ClarinetBb');
clarinet_f0 = pcaData.F0(clarinet_idx);
clarinet_mp = freq2muspitch(clarinet_f0);
clarinet_scores = pcaData.score(clarinet_idx,1:2);

note_idx = zeros(44,1);
for ii = 1:numel(clarinet_notes)
    idx = strcmp(clarinet_mp, clarinet_notes{ii});
    note_idx = note_idx + idx;
    score_old(:,ii) = clarinet_scores(strcmp(clarinet_mp, clarinet_notes{ii}),:);
end
note_idx = logical(note_idx);

[figh, pp] = plotHelper;

scatter(pcaData.score(:,1), ...
    pcaData.score(:,2), ...
    pp.markersize, ...
    [0.7 0.7 0.7], ...
    "Marker", ".");
scatter(score_old(1,:), ...
    score_old(2,:), ...
    pp.markersize, ...
    c1(1:numel(clarinet_notes), :), ...
    "filled", ...
    "Marker", "o");

set(gca, "FontSize", pp.fsize, ...
    "LineWidth", pp.linewidth, ...
    "Layer", "top");

% x-axis
xlim([-9 8])
xlabel('Principal component 1', ...
    "FontSize", pp.fsize, ...
    "Interpreter", "tex", ...
    "FontName", pp.fname);

% y-axis
ylim([-6 7])
ylabel('Principal component 2', ...
    "FontSize", pp.fsize, ...
    "Interpreter", "tex", ...
    "FontName", pp.fname);

title('Component space', ...
    "FontSize", pp.fsize)

fig_name = 'clarinet_pca';
printHelper(figh, pp, fig_folder, fig_name)

%% Extracted envelopes
clarinet_notes = {'FIS3', 'FIS4', 'FIS5'};
clarinet_notes = {'FIS3', 'FIS5'};
clarinet_idx = strcmp(pcaData.instrLabel, 'ClarinetBb');
clarinet_f0 = pcaData.F0(clarinet_idx);
clarinet_mp = freq2muspitch(clarinet_f0);
clarinet_scores = pcaData.score(clarinet_idx,1:2);

for ii = 1:numel(clarinet_notes)
    clarinet_score(:,ii) = clarinet_scores(strcmp(clarinet_mp, clarinet_notes{ii}),:);
end

ss_hat = clarinet_score' * pcaData.coeff(:,1:2)' + pcaData.mu;
ss_hat = ss_hat';
ss_hat = [-20*ones(1,numel(clarinet_notes)); ss_hat; zeros(128-13,numel(clarinet_notes))];

SS_hat = 20.*idct(ss_hat);

[figh, pp] = plotHelper;

plot(cf, SS_hat, ...
    "LineWidth", pp.linewidth);

set(gca, "FontSize", pp.fsize, ...
    "LineWidth", pp.linewidth, ...
    "Layer", "top", ...
    "XScale", "log");

% x-axis
xlim([20 12000])
xticks([100 1000 10000])
xticklabels({'100', '1k', '10k'})
xlabel('Frequency / Hz', ...
    "FontSize", pp.fsize+1, ...
    "Interpreter", "tex", ...
    "FontName", pp.fname);

% y-axis
ylim([-80 20])
ylabel('Level / dB', ...
    "FontSize", pp.fsize, ...
    "Interpreter", "tex", ...
    "FontName", pp.fname);

title('Reconstructed spectral envelope', "FontSize", pp.fsize)

chH = get(gca,'Children');
set(gca,'Children',[chH(end);chH(1:end-1)]);

fig_name = 'clarinet_spectral_envelope_pca';
printHelper(figh, pp, fig_folder, fig_name)


%% Plot CCs for Exp 1

% Plot cepstral coefficients
[figh, pp] = plotHelper;

bp2 = bar(ss_hat(2:13,:));

for ii = 1:size(bp2, 2)
    bp2(ii).FaceColor = c1(ii,:);
    bp2(ii).EdgeColor = c1(ii,:);
    bp2(ii).BaseLine.Visible = "off";
end

hold on

plot([-1 14], [0 0], ...
    "Color", [0 0 0], ...
    "LineWidth", 0.5);

set(gca, "FontSize", pp.fsize, ...
    "LineWidth", pp.linewidth, ...
    "Layer", "top");

% x-axis
xlim([-1 13])
xticks(1:12)
xticklabels({'2', '' '4', '', '6', '', '8', '', '10', '', '12', ''})
xtickangle(0)
xlabel('Cepstral coefficient number', ...
    "FontSize", pp.fsize+1, ...
    "Interpreter", "tex", ...
    "FontName", pp.fname);

% y-axis
ylim([-5 25])
ylabel('Level', ...
    "FontSize", pp.fsize, ...
    "Interpreter", "tex", ...
    "FontName", pp.fname);

title('Exp 1: EFCCs 2-13', "FontSize", pp.fsize)

fig_name = 'clarinet_cepstral_coefficients_exp1';
printHelper(figh, pp, fig_folder, fig_name)

%% Plot CCs for Exp 2

% Plot cepstral coefficients
[figh, pp] = plotHelper;

% Convert decibels
CC_plot = CC ./ 20;

bp1 = bar(CC_plot(1:13,:));

for ii = 1:size(bp1, 2)
    bp1(ii).FaceColor = c2(ii,:);
    bp1(ii).EdgeColor = c2(ii,:);
    bp1(ii).BaseLine.Visible = "off";
end

hold on

plot([0 14], [0 0], ...
    "Color", [0 0 0], ...
    "LineWidth", 0.5);

set(gca, "FontSize", pp.fsize, ...
    "LineWidth", pp.linewidth, ...
    "Layer", "top");

% x-axis
xlim([0 14])
xticks(1:13)
xticklabels({'1', '2', '', '4', '', '6', '', '8', '', '10', '', '12', ''})
xtickangle(0)
xlabel('Cepstral coefficient number', ...
    "FontSize", pp.fsize+1, ...
    "Interpreter", "tex", ...
    "FontName", pp.fname);

% y-axis
ylim([-35 15])
ylabel('Level', ...
    "FontSize", pp.fsize, ...
    "Interpreter", "tex", ...
    "FontName", pp.fname);

title('Exp 2: EFCCs 1-13', "FontSize", pp.fsize)

fig_name = 'clarinet_cepstral_coefficients_exp2';
printHelper(figh, pp, fig_folder, fig_name)

%% Plot spectral envelopes for exp 1
clarinet_notes = {'FIS3', 'FIS4', 'FIS5'};
clarinet_notes = {'FIS3', 'FIS5'};
clarinet_idx = strcmp(pcaData.instrLabel, 'ClarinetBb');
clarinet_f0 = pcaData.F0(clarinet_idx);
clarinet_mp = freq2muspitch(clarinet_f0);
clarinet_scores = pcaData.score(clarinet_idx,1:2);

for ii = 1:numel(clarinet_notes)
    clarinet_score(:,ii) = clarinet_scores(strcmp(clarinet_mp, clarinet_notes{ii}),:);
end

ss_hat = clarinet_score' * pcaData.coeff(:,1:2)' + pcaData.mu;
ss_hat = ss_hat';
ss_hat = [-20*ones(1,numel(clarinet_notes)); ss_hat; zeros(128-13,numel(clarinet_notes))];

SS_hat = 20.*idct(ss_hat);

[figh, pp] = plotHelper;


for ii = 1:size(SS_hat, 2)
    pl2 = plot(cf, SS_hat(:,ii), ...
        "LineWidth", pp.linewidth, ...
        "LineStyle", "-", ...
        "Color", c1(ii,:));
end
chH = get(gca,'Children');
chH(1:size(SS_hat, 2)) = flipud(chH(1:size(SS_hat, 2)));
set(gca,'Children',chH);

for ii = 1:size(SS, 2)
    pl3 = plot(cf, SS(:,ii), ...
        "LineWidth", pp.linewidth, ...
        "LineStyle", ":", ...
        "Color", c1(ii,:));
end
chH = get(gca,'Children');
chH(1:size(SS, 2)) = flipud(chH(1:size(SS, 2)));
set(gca,'Children',chH);


set(gca, "FontSize", pp.fsize, ...
    "LineWidth", pp.linewidth, ...
    "Layer", "top", ...
    "XScale", "log");

% x-axis
xlim([20 12000])
xticks([100 1000 10000])
xticklabels({'100', '1k', '10k'})
xlabel('Frequency / Hz', ...
    "FontSize", pp.fsize+1, ...
    "Interpreter", "tex", ...
    "FontName", pp.fname);

% y-axis
ylim([-80 20])
ylabel('Level / dB', ...
    "FontSize", pp.fsize, ...
    "Interpreter", "tex", ...
    "FontName", pp.fname);

title('Exp 1: spectral envelopes', "FontSize", pp.fsize)

fig_name = 'clarinet_spectral_envelope_exp1';
printHelper(figh, pp, fig_folder, fig_name)

%% Plot spectral envelopes for exp 2
[figh, pp] = plotHelper;

for ii = 1:size(SS_og, 2)
    pl1 = plot(cf, SS_og(:,ii), ...
        "LineWidth", pp.linewidth, ...
        "LineStyle", ":", ...
        "Color", c2(ii,:));
end
chH = get(gca,'Children');
chH(1:size(SS_og, 2)) = flipud(chH(1:size(SS_og, 2)));
set(gca,'Children',chH);

for ii = 1:size(SE, 2)
    pl3 = plot(cf, SE(:,ii), ...
        "LineWidth", pp.linewidth, ...
        "LineStyle", "-", ...
        "Color", c2(ii,:));
end
chH = get(gca,'Children');
chH(1:size(SE, 2)) = flipud(chH(1:size(SE, 2)));
set(gca,'Children',chH);

set(gca, "FontSize", pp.fsize, ...
    "LineWidth", pp.linewidth, ...
    "Layer", "top", ...
    "XScale", "log");

% x-axis
xlim([20 12000])
xticks([100 1000 10000])
xticklabels({'100', '1k', '10k'})
xlabel('Frequency / Hz', ...
    "FontSize", pp.fsize+1, ...
    "Interpreter", "tex", ...
    "FontName", pp.fname);


% y-axis
ylim([-80 10])
ylabel('Level / dB', ...
    "FontSize", pp.fsize, ...
    "Interpreter", "tex", ...
    "FontName", pp.fname);

title('Exp 2: spectral envelopes', "FontSize", pp.fsize)

fig_name = 'clarinet_spectral_envelope_exp2';
printHelper(figh, pp, fig_folder, fig_name)

%% Plot synthesized waveforms for exp 2
audioFolder = strcat(projectFolder,filesep,'sounds',filesep,'exp2',filesep,'ClarinetBb');
fileNames = {'ClarinetBb_all_FIS3.wav', 'ClarinetBb_all_FIS4.wav', 'ClarinetBb_all_FIS5.wav'};
fileNames = {'ClarinetBb_all_FIS3.wav', 'ClarinetBb_all_FIS5.wav'};

% Read audio
for iFile = 1:numel(fileNames)
    [s_synth(:,iFile), fs] = audioread(strcat(audioFolder,filesep,fileNames{iFile}));

end

s_plot_synth = s_synth(1:end-1,:);
s_plot_synth = [s_plot_synth; zeros(fs/2,numel(fileNames))];
if numel(fileNames) == 2
    % Adjust wave 2
    s_plot_synth(:,2) = 0;
    s_plot_synth(fs*0.5:end-1, 2) = s_synth(1:end-1, 2);
elseif numel(fileNames) == 3
    % Adjust wave 3
    s_plot_synth(:,3) = 0;
    s_plot_synth(fs*0.5:end-1, 3) = s_synth(1:end-1, 3);
end

% Plot wave form
[figh, pp] = plotHelper;

T = size(s_plot_synth,1);
p = plot(s_plot_synth(1:T,:)./max(s_plot_synth), ...
    "LineWidth", 0.5);

for ii = 1:size(p, 1)
    p(ii).Color = c2(ii,:);
end

set(gca, "FontSize", pp.fsize, ...
    "LineWidth", pp.linewidth, ...
    "Layer", "top");

% x-axis
xlim([0 fs])
xticks([0 0.25*fs 0.5*fs 0.75*fs fs])
xticklabels({'0', '0.25', '0.5', '0.75', '1'})
xlabel('Time / s', ...
    "FontSize", pp.fsize+1, ...
    "Interpreter", "tex", ...
    "FontName", pp.fname);

% y-axis
ylim([-1.1 1.1])
ylabel('Amplitude', ...
    "FontSize", pp.fsize, ...
    "Interpreter", "tex", ...
    "FontName", pp.fname);

title('Waveform', "FontSize", pp.fsize)

chH = get(gca,'Children');
set(gca,'Children',[chH(end);chH(1:end-1)]);

fig_name = 'clarinet_waveform_synthesized_exp2';
printHelper(figh, pp, fig_folder, fig_name)

%% Plot synthesized waveforms for exp 1
audioFolder = strcat(projectFolder,filesep, ...
    'sounds',filesep, ...
    'exp1',filesep, ...
    'sound_quality_rating',filesep, ...
    'congruentF0');
fileNames = {'congruentF0_34.wav', 'congruentF0_38.wav', 'congruentF0_37.wav'};
fileNames = {'congruentF0_34.wav', 'congruentF0_37.wav'};

% Read audio
for iFile = 1:numel(fileNames)
    [s_synth(:,iFile), fs] = audioread(strcat(audioFolder,filesep,fileNames{iFile}));

end

s_plot_synth = s_synth(1:end-1,:);
s_plot_synth = [s_plot_synth; zeros(fs/2,numel(fileNames))];
if numel(fileNames) == 2
    % Adjust wave 2
    s_plot_synth(:,2) = 0;
    s_plot_synth(fs*0.5:end-1, 2) = s_synth(1:end-1, 2);
elseif numel(fileNames) == 3
    % Adjust wave 3
    s_plot_synth(:,3) = 0;
    s_plot_synth(fs*0.5:end-1, 3) = s_synth(1:end-1, 3);
end

% Plot wave form
[figh, pp] = plotHelper;

T = size(s_plot_synth,1);
p = plot(s_plot_synth(1:T,:)./max(s_plot_synth), ...
    "LineWidth", 0.5);

for ii = 1:size(p, 1)
    p(ii).Color = c1(ii,:);
end

set(gca, "FontSize", pp.fsize, ...
    "LineWidth", pp.linewidth, ...
    "Layer", "top");

% x-axis
xlim([0 fs])
xticks([0 0.25*fs 0.5*fs 0.75*fs fs])
xticklabels({'0', '0.25', '0.5', '0.75', '1'})
xlabel('Time / s', ...
    "FontSize", pp.fsize+1, ...
    "Interpreter", "tex", ...
    "FontName", pp.fname);

% y-axis
ylim([-1.1 1.1])
ylabel('Amplitude', ...
    "FontSize", pp.fsize, ...
    "Interpreter", "tex", ...
    "FontName", pp.fname);

title('Waveform', "FontSize", pp.fsize)

chH = get(gca,'Children');
set(gca,'Children',[chH(end);chH(1:end-1)]);

fig_name = 'clarinet_waveform_synthesized_exp1';
printHelper(figh, pp, fig_folder, fig_name)