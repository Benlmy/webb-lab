function out = jwst_psf_mft(P, lambda, opd, fovArcsec, nPix)
%JWST_PSF_MFT PSF over a chosen field of view via a matrix Fourier transform.
%   out = JWST_PSF_MFT(P, lambda, opd, fovArcsec, nPix) evaluates the
%   Fraunhofer integral directly on an nPix-by-nPix grid spanning
%   fovArcsec, so small fields can be sampled finely without the huge
%   zero-padded FFT that JWST_PSF would need. Normalisation and output
%   fields match JWST_PSF (psf, strehl, theta, pixelScale, rmsWfe).

    if nargin < 3 || isempty(opd)
        opd = zeros(size(P.mask));
    end
    if nargin < 4 || isempty(fovArcsec)
        fovArcsec = 2;
    end
    if nargin < 5 || isempty(nPix)
        nPix = 256;
    end

    arcsecPerRad = 180/pi*3600;
    theta = ((0:nPix-1) - floor(nPix/2)) * (fovArcsec / nPix);
    thetaRad = theta / arcsecPerRad;

    xs = P.x(1, :);
    ys = P.y(:, 1)';
    Ax = exp(-1i * 2*pi / lambda * (thetaRad(:) * xs));
    Ay = exp(-1i * 2*pi / lambda * (thetaRad(:) * ys));

    field = P.mask .* exp(1i * 2*pi * opd / lambda);
    E = Ay * field * Ax.';

    out.psf = abs(E).^2 / nnz(P.mask)^2;
    out.strehl = max(out.psf(:));
    out.theta = theta;
    out.pixelScale = fovArcsec / nPix;
    w = opd(P.mask);
    out.rmsWfe = sqrt(mean((w - mean(w)).^2));
end
