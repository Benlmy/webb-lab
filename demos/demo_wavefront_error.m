%% DEMO_WAVEFRONT_ERROR  How mirror errors erode image quality
% Applies random piston errors to the 18 segments, measures the resulting
% Strehl ratio (peak brightness relative to a perfect mirror) and compares
% it with the Marechal approximation S ~ exp(-(2*pi*sigma/lambda)^2).
% JWST's budget of 150 nm RMS wavefront error gives S ~ 0.8 at 2 um, the
% conventional definition of a diffraction-limited telescope.

c = jwst_constants();
P = jwst_pupil(256, c);
lambda = c.diffractionLimitWavelength;
n = size(P.seg.center, 1);
rng(42);

pistonRms = linspace(0, 300e-9, 13);
trials = 8;
wfe = zeros(numel(pistonRms), trials);
strehl = zeros(numel(pistonRms), trials);
for i = 1:numel(pistonRms)
    for j = 1:trials
        opd = jwst_segment_wavefront(P, pistonRms(i) * randn(n, 1));
        out = jwst_psf(P, lambda, opd, 4);
        wfe(i, j) = out.rmsWfe;
        strehl(i, j) = out.strehl;
    end
end

sig = linspace(0, max(wfe(:)), 200);
marechal = exp(-(2*pi*sig/lambda).^2);

fig = new_figure('Wavefront error', [900 560]);
plot(wfe(:)*1e9, strehl(:), 'o', 'MarkerSize', 5, ...
     'MarkerFaceColor', [0.25 0.45 0.85], 'MarkerEdgeColor', 'none');
hold on;
plot(sig*1e9, marechal, 'k-', 'LineWidth', 1.6);
yl = [0 1.05];
plot([1 1] * c.wfeRequirement*1e9, yl, '--', 'Color', [0.85 0.3 0.1], 'LineWidth', 1.4);
text(c.wfeRequirement*1e9 + 5, 0.95, 'JWST budget: 150 nm', 'Color', [0.85 0.3 0.1]);
hold off; grid on; ylim(yl);
xlabel('RMS wavefront error [nm]');
ylabel('Strehl ratio at 2 \mum');
legend('Monte Carlo (random segment pistons)', 'Marechal approximation', 'Location', 'southwest');
title('Image quality vs. segment phasing error');
save_figure(fig, 'strehl_vs_wfe');

fprintf('\n== Wavefront error budget ==\n');
fprintf('Marechal Strehl at 150 nm RMS, 2 um: %.2f\n', exp(-(2*pi*c.wfeRequirement/lambda)^2));
