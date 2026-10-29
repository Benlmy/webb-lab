function c = jwst_constants()
%JWST_CONSTANTS Physical and mission constants used throughout the project.
%   c = JWST_CONSTANTS() returns a struct of SI-unit constants describing the
%   James Webb Space Telescope, its environment at Sun-Earth L2, and a few
%   reference values (Hubble) used for comparison.

    % Fundamental
    c.sigma   = 5.670374419e-8;     % Stefan-Boltzmann constant [W m^-2 K^-4]
    c.h       = 6.62607015e-34;     % Planck constant [J s]
    c.kB      = 1.380649e-23;       % Boltzmann constant [J/K]
    c.cLight  = 299792458;          % speed of light [m/s]
    c.AU      = 1.495978707e11;     % astronomical unit [m]
    c.solarConstant = 1361;         % solar irradiance at 1 AU [W/m^2]
    c.arcsecPerRad  = 180/pi*3600;

    % Sun-Earth/Moon system for the circular restricted three-body problem
    c.muSun   = 1.32712440018e20;   % GM of the Sun [m^3/s^2]
    c.muEarthMoon = 4.0350323550e14;% GM of Earth + Moon [m^3/s^2]
    c.mu      = c.muEarthMoon / (c.muSun + c.muEarthMoon); % mass ratio
    c.yearDays = 365.256363;        % sidereal year [days]

    % Optical telescope element
    c.segmentFlatToFlat = 1.32;     % hexagonal segment size, flat to flat [m]
    c.segmentGap        = 0.007;    % edge-to-edge gap between segments [m]
    c.numSegments       = 18;
    c.secondaryDiameter = 0.74;     % secondary mirror diameter [m]
    c.strutWidth        = 0.08;     % projected width of a secondary support strut [m]
    c.strutAnglesDeg    = [90 240 300]; % strut directions as seen from the pupil centre
    c.focalLength       = 131.4;    % effective focal length of the OTE [m]
    c.wfeRequirement    = 150e-9;   % total wavefront error budget at NIRCam [m RMS]
    c.diffractionLimitWavelength = 2e-6; % diffraction-limited at and beyond 2 um [m]

    % Spectral coverage
    c.lambdaMin = 0.6e-6;           % [m]
    c.lambdaMax = 28.5e-6;          % [m]

    % Sunshield
    c.sunshieldLength = 21.197;     % [m]
    c.sunshieldWidth  = 14.162;     % [m]
    c.sunshieldLayers = 5;

    % Hubble, for comparison
    c.hubbleDiameter = 2.4;         % [m]
    c.hubbleArea     = 4.5;         % collecting area [m^2]
end
