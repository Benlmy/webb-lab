function B = planck_radiance(lambda, T)
%PLANCK_RADIANCE Blackbody spectral radiance.
%   B = PLANCK_RADIANCE(lambda, T) returns the spectral radiance
%   [W m^-2 sr^-1 m^-1] at wavelength lambda [m] and temperature T [K].
%   lambda and T are broadcast against each other.

    h = 6.62607015e-34;
    c = 299792458;
    kB = 1.380649e-23;
    B = 2*h*c^2 ./ lambda.^5 ./ expm1(h*c ./ (lambda .* kB .* T));
end
