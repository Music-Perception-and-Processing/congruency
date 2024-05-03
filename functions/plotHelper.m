function [figh, pp] = plotHelper
%PLOTHELPER plots the input figure in a nice way
%   Detailed explanation goes here

% Basic plot settings
pp.visible = 'off';
pp.print = 1;
pp.fsize = 20;
pp.figwidth = 15;
pp.figheight = 10;
pp.linewidth = 1;
pp.markersize = 100;
pp.fname = 'Times New Roman';

% More plot settings
set(0, "defaultAxesFontName", pp.fname);
set(0, "defaultTextInterpreter", "tex");
set(0, "defaultAxesTickLabelInterpreter", "tex");
set(0, "defaultLegendInterpreter", "tex");

figh = figure("Visible", pp.visible, ...
    "DefaultTextFontName", pp.fname, ...
    "DefaultAxesFontName", pp.fname);
axh = axes("Parent", figh);
hold(axh, "on");
box(axh,"on");

end

