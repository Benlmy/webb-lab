function fig = new_figure(name, sizePx)
%NEW_FIGURE Create a figure that stays hidden when running headless.
%   fig = NEW_FIGURE(name, [width height]) opens a figure window, or an
%   invisible figure when JWST_HEADLESS is set to 1 in the environment
%   (run_all and the tests do this so batch runs never pop up windows).

    if nargin < 2
        sizePx = [1000 600];
    end
    visible = 'on';
    if strcmp(getenv('JWST_HEADLESS'), '1')
        visible = 'off';
    end
    fig = figure('Name', name, 'NumberTitle', 'off', 'Color', 'w', ...
                 'Visible', visible, 'Position', [80 80 sizePx]);
end
