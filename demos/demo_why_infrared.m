%% DEMO_WHY_INFRARED  The science case: cold telescope, big mirror, red light
% Three panels:
%   1. Blackbody spectra: cool objects (dust, planets) and the telescope
%      itself glow in the infrared, so the telescope must be colder than
%      its targets.
%   2. Angular resolution 1.22*lambda/D for JWST (6.5 m) and Hubble (2.4 m).
%      JWST's larger mirror keeps infrared images sharp.
%   3. Expansion of the universe stretches the Lyman break of the first
%      galaxies from the ultraviolet into NIRCam's infrared filters.

c = jwst_constants();
lam = logspace(log10(0.2e-6), log10(300e-6), 500);

%% Blackbody spectra (each normalised to its peak)
sources = {
    'Sun-like star (5800 K)',     5800
    'Habitable planet (288 K)',    288
    'Sunshield hot side (354 K)',  354
    'Telescope mirrors (40 K)',     40
};
fig = new_figure('Why infrared', [1400 460]);
subplot(1, 3, 1);
cols = [0.95 0.65 0.1; 0.2 0.6 0.3; 0.85 0.25 0.15; 0.25 0.45 0.85];
for k = 1:size(sources, 1)
    B = planck_radiance(lam, sources{k, 2});
    semilogx(lam*1e6, B / max(B), 'LineWidth', 1.6, 'Color', cols(k, :));
    hold on;
end
yl = [0 1.75];
patch([c.lambdaMin c.lambdaMax c.lambdaMax c.lambdaMin]*1e6, [yl(1) yl(1) yl(2) yl(2)], ...
      [0.5 0.5 0.5], 'FaceAlpha', 0.12, 'EdgeColor', 'none');
hold off; grid on; ylim([0 1.75]); xlim([0.2 300]);
xlabel('wavelength [\mum]'); ylabel('normalised radiance');
legend([sources(:, 1); {'JWST coverage 0.6-28.5 \mum'}], 'Location', 'north', 'FontSize', 8);
title('Blackbody spectra');

%% Resolution vs wavelength
subplot(1, 3, 2);
lamR = linspace(0.6e-6, 28.5e-6, 300);
resJ = 1.22 * lamR / 6.5 * c.arcsecPerRad;
resH = 1.22 * lamR / c.hubbleDiameter * c.arcsecPerRad;
hubbleIR = lamR <= 2.5e-6;
loglog(lamR*1e6, resJ, 'LineWidth', 1.8, 'Color', [0.9 0.6 0.1]); hold on;
loglog(lamR(hubbleIR)*1e6, resH(hubbleIR), 'LineWidth', 1.8, 'Color', [0.3 0.4 0.8]);
hold off; grid on; xlim([0.5 30]);
xlabel('wavelength [\mum]'); ylabel('diffraction limit 1.22\lambda/D [arcsec]');
legend('JWST (6.5 m)', 'Hubble (2.4 m, to 2.5 \mum)', 'Location', 'northwest');
title('Sharpness of the diffraction limit');

%% Lyman-break redshift and NIRCam filters
subplot(1, 3, 3);
F = nircam_filters();
z = linspace(5, 25, 300);
lamObs = 0.12157 * (1 + z);
hold on;
for k = 1:numel(F)
    patch([5 25 25 5], [F(k).lo F(k).lo F(k).hi F(k).hi], [0.4 0.55 0.85], ...
          'FaceAlpha', 0.10 + 0.05*mod(k, 2), 'EdgeColor', 'none');
    text(5.4, F(k).pivot, F(k).name, 'FontSize', 8, 'Color', [0.15 0.25 0.55]);
end
plot(z, lamObs, 'k-', 'LineWidth', 1.8);
hold off; box on; grid on; xlim([5 25]); ylim([0.7 3.25]);
xlabel('redshift z'); ylabel('observed Lyman break [\mum]');
title('First galaxies drop out of blue filters');

save_figure(fig, 'why_infrared');

fprintf('\n== Lyman-break dropouts in NIRCam ==\n');
fprintf('  z    break[um]  invisible in   visible in\n');
for zz = [6 8 10 12 14 16 20]
    [drop, det, lo] = lyman_break_filter(zz);
    if isempty(drop), drop = '-'; end
    if isempty(det), det = '-'; end
    fprintf('%4d   %7.2f    %-12s   %s\n', zz, lo, drop, det);
end
