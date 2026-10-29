function file = save_figure(fig, name)
%SAVE_FIGURE Export a figure as a PNG into the project's output folder.
%   file = SAVE_FIGURE(fig, name) writes output/<name>.png and returns the
%   full path. The output folder is created if needed.

    root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    outDir = fullfile(root, 'output');
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end
    file = fullfile(outDir, [name '.png']);
    print(fig, file, '-dpng', '-r130');
    fprintf('  saved %s\n', file);
end
