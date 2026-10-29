function F = nircam_filters()
%NIRCAM_FILTERS NIRCam wide-band filters used for high-redshift surveys.
%   F = NIRCAM_FILTERS() returns a struct array with fields name, pivot,
%   lo and hi (pivot wavelength and approximate half-power band edges, in
%   micrometres). Band edges are rounded values from the JWST user
%   documentation and are adequate for planning-level calculations.

    data = {
        'F090W', 0.902, 0.795, 1.005
        'F115W', 1.154, 1.013, 1.282
        'F150W', 1.501, 1.331, 1.668
        'F200W', 1.989, 1.755, 2.227
        'F277W', 2.762, 2.423, 3.132
        'F356W', 3.568, 3.135, 3.981
        'F444W', 4.408, 3.881, 4.982
    };
    F = struct('name', data(:, 1), 'pivot', data(:, 2), 'lo', data(:, 3), 'hi', data(:, 4));
end
