%% RUN_ALL  Run every JWST demo and save the figures to output/
% Works in MATLAB (R2019b or newer) and GNU Octave (7 or newer).
%   >> run_all
% Set the environment variable JWST_HEADLESS=1 to keep figure windows hidden.

run(fullfile(fileparts(mfilename('fullpath')), 'startup.m'));

demos = {
    'demo_mirror_psf'
    'demo_segment_alignment'
    'demo_wavefront_error'
    'demo_l2_halo_orbit'
    'demo_sunshield_thermal'
    'demo_why_infrared'
};

fprintf('James Webb Space Telescope: physics demos\n');
fprintf('=========================================\n');
for iDemo = 1:numel(demos)
    tStart = tic;
    run(demos{iDemo});
    fprintf('  (%s finished in %.1f s)\n', demos{iDemo}, toc(tStart));
    if strcmp(getenv('JWST_HEADLESS'), '1')
        close all;
    end
end
fprintf('\nAll demos complete. Figures are in the output/ folder.\n');
