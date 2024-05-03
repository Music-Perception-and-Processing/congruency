function printHelper(figh, pp, fig_folder, fig_name)
%PRINTHELPER prints the given figure in a nice way
%   Detailed explanation goes here

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
    print(figh, [fig_folder, filesep, fig_name, '.pdf'], "-dpdf");
end

end

