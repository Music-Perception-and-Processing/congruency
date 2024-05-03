clear
close all
clc

run plot_properties.m

set(0, "defaultAxesFontName", pp.fname);
set(0, "defaultTextInterpreter", "tex");
set(0, "defaultAxesTickLabelInterpreter", "tex");
set(0, "defaultLegendInterpreter", "tex");

load("data/pca_data_v5.mat");

figh = figure("Visible", pp.visible, ...
    "DefaultTextFontName", pp.fname, ...
    "DefaultAxesFontName", pp.fname);
axh = axes("Parent", figh);
hold(axh, "on");
box(axh, "on");

% Variables
score = pcaData.score;
coeff = pcaData.coeff;
mu = pcaData.mu;

fs = 44100;
numBands = 128;
[~, fERB, ~] = designAuditoryFilterBank(fs, ...
    "FrequencyScale", "erb", ...
    "FFTLength", fs-1, ...
    "NumBands", numBands, ...
    "FrequencyRange", [20 12000], ...
    "Normalization", "bandwidth");

x = [-10 -5 0 5 10];
x = fliplr(x);

y = [-10 -5 0 5 10];

[Y, X] = meshgrid(x, y);
datap = [X(:), Y(:)];
nDP = length(datap(:, 1));

CC1 = -19.95 + (-1.1851).*datap(:,1)...
    + (-0.3267).*datap(:,2)...
    + (-0.0851).*datap(:,1).*datap(:,2);

reX = datap * coeff(:, 1:2)' + mu;
reEFCC = [CC1 .* ones(nDP,1), reX, zeros(nDP, 115)] .* 20;
reSE = idct(reEFCC, [], 2);
reSEamp = 10.^(reSE./20);
reSEamp = reSEamp./max(reSEamp, [], 2);

% Compute spectral centroids
for iDP = 1:nDP
    SC(iDP) = sum(fERB .* reSEamp(iDP, :)) / sum(reSEamp(iDP, :));
end

% Plotting
tl = tiledlayout(5, 5, "TileSpacing", "compact", "Padding", "compact");
xlabel(tl, 'Frequency / kHz', ...
    "FontSize", pp.fsize+1, ...
    "Interpreter", "tex", ...
    "FontName", pp.fname);
ylabel(tl, 'Level / dB', ...
    "FontSize", pp.fsize+1, ...
    "Interpreter", "tex", ...
    "FontName", pp.fname);

fMin = 100;
fMax = 1500;
colorValues = logspace(log10(fMin), log10(fMax), 1000);
colors = turbo(numel(colorValues));

for iDP = 1:nDP
    h(iDP) = nexttile(tl);

    % Find correct color value for plot
    [minVal, minIdx] = min(abs(colorValues - SC(iDP)));

    plot(fERB, reSE(iDP, :), ...
        "LineWidth", pp.linewidth, ...
        "Color", colors(minIdx, :));

    % SC text
    text(28, -70, ['SC   = ', num2str(round(SC(iDP))), ' Hz'], ...
        "FontSize", pp.fsize-5, ...
        "FontWeight", "normal");

    % PC1 text
    text(28, -30, ['PC1 = ', num2str(datap(iDP, 1))], ...
        "FontSize", pp.fsize-5, ...
        "FontWeight", "normal");

    % PC2 text
    text(28, -50, ['PC2 = ', num2str(datap(iDP, 2))], ...
        "FontSize", pp.fsize-5, ...
        "FontWeight", "normal");

    xlim([20 12000])
    if ismember(iDP, [21, 22, 23, 24, 25])
        xticks([100 1000 10000]);
        xticklabels({'0.1', '1', '10'});
    else
        xticks([100 1000 10000]);
        xticklabels({});
    end

    ylim([-90 30])
    if ismember(iDP, [1, 6, 11, 16, 21])
        yticks([-90 -60 -30 0 30]);
    else
        yticks([-90 -60 -30 0 30]);
        yticklabels({});
    end

    set(gca, "FontSize", pp.fsize, ...
        "LineWidth", pp.linewidth, ...
        "Layer", "top", ...
        "XScale", "log");

end

set(h, "Colormap", turbo);
set(gca, "ColorScale", "log");
clim([100 1500])
cbh = colorbar(h(end), ...
    "Ticks", [100 500 1000 1500], ...
    "TickLabels", {'0.1', '0.5', '1', '1.5'}, ...
    "LineWidth", pp.linewidth, ...
    "FontSize", pp.fsize-2);
cbh.Layout.Tile = "east";
cbh.Label.String = 'Spectral centroid (SC) / kHz';
cbh.Label.FontSize = pp.fsize+1;
cbh.Label.Interpreter = "tex";
set(cbh, "TickLabelInterpreter", "tex");


if isfield(pp, "figwidth")
    if ~isempty(pp.figwidth)
        set(figh, "PaperPositionMode", "manual");
        set(figh, "PaperUnits", "centimeters");
        set(figh, "PaperPosition", [0 0 pp.figwidth pp.figheight]);
        set(figh, "Toolbar", "none");
        set(figh, "Menubar", "none");
        set(figh, "PaperSize", [pp.figwidth pp.figheight]);

    end

end

if pp.print
    print(figh, [fig_folder, filesep, 'SEonPCSpace.pdf'], "-dpdf");
end