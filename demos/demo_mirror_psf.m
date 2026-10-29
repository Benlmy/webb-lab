%% DEMO_MIRROR_PSF  JWST's segmented primary mirror and its diffraction pattern
% Builds the 18-segment pupil, then computes the point-spread function
% (the image of a single star) at three NIRCam wavelengths. The hexagonal
% segment edges produce six bright diffraction spikes, and the vertical
% secondary-mirror strut adds two fainter horizontal ones: the familiar
% eight-pointed JWST star.

c = jwst_constants();
P = jwst_pupil(512, c);

fprintf('\n== Primary mirror ==\n');
fprintf('Segments:              %d\n', size(P.seg.center, 1));
fprintf('Clear collecting area: %.1f m^2 (Hubble: %.1f m^2, x%.1f)\n', ...
        P.area, c.hubbleArea, P.area / c.hubbleArea);

%% Pupil map
fig1 = new_figure('JWST pupil', [760 700]);
imagesc(P.x(1, :), P.y(:, 1), P.mask);
axis image; set(gca, 'YDir', 'normal');
colormap(gca, [0.08 0.09 0.13; 0.93 0.75 0.30]);
hold on;
for k = 1:size(P.seg.center, 1)
    text(P.seg.center(k, 1), P.seg.center(k, 2), P.seg.label{k}, ...
         'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', [0.2 0.1 0], ...
         'BackgroundColor', [0.93 0.75 0.30], 'Margin', 1);
end
hold off;
xlabel('x [m]'); ylabel('y [m]');
title(sprintf('JWST entrance pupil: 18 segments, %.1f m^2 clear area', P.area));
save_figure(fig1, 'mirror_pupil');

%% PSFs at three wavelengths, shown over the same 2" patch of sky
lambdas = [0.9 2.0 4.4] * 1e-6;
fig2 = new_figure('JWST PSF', [1350 500]);
for i = 1:numel(lambdas)
    out = jwst_psf_mft(P, lambdas(i), [], 2.0, 300);

    subplot(1, 3, i);
    imagesc(out.theta, out.theta, log10(max(out.psf, 1e-6)));
    axis image; set(gca, 'YDir', 'normal');
    colormap(gca, hot); caxis([-4.5 0]);
    xlabel('arcsec'); ylabel('arcsec');
    title(sprintf('\\lambda = %.1f \\mum   (\\lambda/D = %.0f mas)', ...
          lambdas(i)*1e6, lambdas(i)/6.5*c.arcsecPerRad*1e3));

    % FWHM from a horizontal cut through the peak
    mid = floor(size(out.psf, 1)/2) + 1;
    cut = out.psf(mid, :);
    above = find(cut >= 0.5 * max(cut));
    fwhm = (above(end) - above(1) + 1) * out.pixelScale;
    fprintf('PSF at %.1f um: FWHM ~ %.0f mas\n', lambdas(i)*1e6, fwhm*1e3);
end
subplot(1, 3, 1);
text(-0.95, -0.88, 'log scale, 4.5 decades', 'Color', 'w', 'FontSize', 8);
save_figure(fig2, 'psf_wavelengths');
