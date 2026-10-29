% STARTUP Put the JWST project folders on the MATLAB/Octave path.
%   MATLAB runs this automatically when started in this folder. In Octave,
%   or from anywhere else, run it once by hand: >> run('path/to/startup.m')

jwstRoot = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(jwstRoot, 'src')));
addpath(fullfile(jwstRoot, 'demos'));
addpath(fullfile(jwstRoot, 'tests'));
clear jwstRoot
