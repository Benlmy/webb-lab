%% DEMO_SUNSHIELD_THERMAL  Five layers from 350 K down to ~36 K
% The tennis-court-sized sunshield blocks sunlight, Earthshine and
% Moonshine so the telescope can cool passively to below 50 K, cold enough
% that its own infrared glow does not swamp faint galaxies. Each layer is
% coated with aluminium (low emissivity), and the layers are spread apart
% at the edges so heat radiates out sideways instead of passing inward.
%
% This is a one-dimensional radiative model: it captures the physics that
% sets the temperature drop, not the full 3D thermal design.

c = jwst_constants();
base = sunshield_equilibrium();
area = c.sunshieldLength * c.sunshieldWidth;

fprintf('\n== Sunshield thermal balance ==\n');
fprintf('Absorbed sunlight on layer 1: %.0f W/m^2 (%.0f kW over %.0f m^2)\n', ...
        base.absorbed, base.absorbed * area / 1e3, area);
for k = 1:numel(base.T)
    fprintf('  Layer %d: %6.1f K (%6.1f C)\n', k, base.T(k), base.T(k) - 273.15);
end
fprintf('Heat leaking from the cold side: %.2f W in total\n', base.heatLeak * area);

nLayers = 1:8;
Tcold = zeros(size(nLayers));
for i = nLayers
    r = sunshield_equilibrium(struct('nLayers', i));
    Tcold(i) = r.T(end);
end

angles = 0:5:60;
Thot = zeros(size(angles));
Tcold5 = zeros(size(angles));
for i = 1:numel(angles)
    r = sunshield_equilibrium(struct('sunAngleDeg', angles(i)));
    Thot(i) = r.T(1);
    Tcold5(i) = r.T(end);
end

fig = new_figure('Sunshield', [1350 460]);
hotToCold = [linspace(0.9, 0.2, 5)' linspace(0.3, 0.45, 5)' linspace(0.15, 0.9, 5)'];

subplot(1, 3, 1);
hold on;
for k = 1:numel(base.T)
    rectangle('Position', [k-0.35 0 0.7 base.T(k)], 'FaceColor', hotToCold(k, :), ...
              'EdgeColor', 'none');
end
hold off; box on; grid on; xlim([0.4 numel(base.T) + 0.6]); xlabel('layer (1 faces the Sun)'); ylabel('temperature [K]');
title('Steady-state layer temperatures');
for k = 1:numel(base.T)
    text(k, base.T(k) + 12, sprintf('%.0f K', base.T(k)), 'HorizontalAlignment', 'center');
end
ylim([0 max(base.T) * 1.15]);

subplot(1, 3, 2);
semilogy(nLayers, Tcold, 'o-', 'LineWidth', 1.6, 'Color', [0.2 0.45 0.85], ...
         'MarkerFaceColor', [0.2 0.45 0.85]);
hold on;
semilogy(5, Tcold(5), 'p', 'MarkerSize', 16, 'MarkerFaceColor', [0.95 0.7 0.1], ...
         'MarkerEdgeColor', 'k');
hold off; grid on;
xlabel('number of layers'); ylabel('cold-side temperature [K]');
title('Why five layers?');
text(5.3, Tcold(5) * 1.4, sprintf('JWST: %.0f K', Tcold(5)));

subplot(1, 3, 3);
plot(angles, Thot, 'LineWidth', 1.6, 'Color', hotToCold(1, :)); hold on;
plot(angles, Tcold5, 'LineWidth', 1.6, 'Color', hotToCold(5, :)); hold off;
grid on; ylim([0 400]);
xlabel('Sun angle from sunshield normal [deg]'); ylabel('temperature [K]');
legend('layer 1 (Sun side)', 'layer 5 (cold side)', 'Location', 'east');
title('Effect of pointing (Sun incidence)');

save_figure(fig, 'sunshield_thermal');
