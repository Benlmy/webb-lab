%% DEMO_SEGMENT_ALIGNMENT  From 18 star images to one: mirror commissioning
% After launch, each of the 18 segments formed its own image of the same
% star. The commissioning team identified each spot, tilted the segments to
% stack the spots on top of each other, and then adjusted segment pistons
% until the segments acted as a single 6.5 m mirror. This demo reproduces
% those three stages at a wavelength of 2 um.

c = jwst_constants();
P = jwst_pupil(512, c);
lambda = 2e-6;
n = size(P.seg.center, 1);
rng(7);

% Stage 1: segments tilted so each spot lands at its own sky position, in a
% pattern mirroring the segment layout (like the February 2022 image).
spreadArcsecPerMeter = 1.3;
tilt = spreadArcsecPerMeter / c.arcsecPerRad * P.seg.center;
opd1 = jwst_segment_wavefront(P, 3e-6 * randn(n, 1), tilt(:, 1), tilt(:, 2));

% Stage 2: spots stacked, but segments still out of phase by ~1 um.
opd2 = jwst_segment_wavefront(P, 1.0e-6 * randn(n, 1), ...
                              0.01/c.arcsecPerRad * randn(n, 1), ...
                              0.01/c.arcsecPerRad * randn(n, 1));

% Stage 3: fine phased, leaving ~20 nm RMS piston error.
opd3 = jwst_segment_wavefront(P, 20e-9 * randn(n, 1));

stages = {opd1, opd2, opd3};
names = {'1. Segment image array', '2. Image stacking', '3. Fine phasing'};
fov = [9 2.4 2.4];

fprintf('\n== Segment alignment at %.1f um ==\n', lambda*1e6);
fig = new_figure('Segment alignment', [1300 470]);
for i = 1:3
    out = jwst_psf_mft(P, lambda, stages{i}, fov(i), 300);
    subplot(1, 3, i);
    imagesc(out.theta, out.theta, log10(max(out.psf, 1e-6)));
    axis image; set(gca, 'YDir', 'normal');
    colormap(gca, hot); caxis([-4.5 0]);
    xlabel('arcsec'); ylabel('arcsec');
    title(sprintf('%s\nStrehl %.3f, WFE %.0f nm RMS', names{i}, out.strehl, out.rmsWfe*1e9));
    fprintf('%-26s Strehl = %.3f\n', names{i}, out.strehl);
end
save_figure(fig, 'segment_alignment');
