function out = jwst_psf(P, lambda, opd, pad)
%JWST_PSF Monochromatic point-spread function by Fraunhofer diffraction.
%   out = JWST_PSF(P, lambda) computes the PSF of pupil P (see JWST_PUPIL)
%   at wavelength lambda [m]. out = JWST_PSF(P, lambda, opd) adds an
%   optical path difference map [m]. pad is the zero-padding factor
%   (default 4), which sets the image sampling to lambda/(pad*6.6 m).
%
%   Fields of out:
%     psf         image normalised so an unaberrated pupil has peak 1
%     strehl      peak of psf (Strehl ratio)
%     theta       angular coordinate vector of the image [arcsec]
%     pixelScale  angular size of an image pixel [arcsec]
%     rmsWfe      RMS wavefront error over the clear aperture [m]

    if nargin < 3 || isempty(opd)
        opd = zeros(size(P.mask));
    end
    if nargin < 4 || isempty(pad)
        pad = 4;
    end

    N = size(P.mask, 1);
    M = pad * N;

    field = zeros(M, M);
    field(1:N, 1:N) = P.mask .* exp(1i * 2*pi * opd / lambda);
    img = abs(fftshift(fft2(field))).^2;

    % The unaberrated peak equals the squared number of clear pixels.
    img = img / nnz(P.mask)^2;

    radPerPix = lambda / (M * P.dx);
    arcsecPerRad = 180/pi*3600;

    out.psf = img;
    out.strehl = max(img(:));
    out.pixelScale = radPerPix * arcsecPerRad;
    out.theta = ((0:M-1) - floor(M/2)) * out.pixelScale;
    w = opd(P.mask);
    out.rmsWfe = sqrt(mean((w - mean(w)).^2));
end
