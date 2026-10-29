%% DEMO_L2_HALO_ORBIT  JWST's halo orbit around Sun-Earth L2
% Models the Sun and Earth-Moon system as the circular restricted
% three-body problem, locates L2 about 1.5 million km beyond Earth, and
% computes a halo orbit similar in size to JWST's (reaching ~800,000 km
% from L2, about six months per revolution) with Richardson's approximation
% followed by differential correction. The orbit is unstable, which is why JWST burns
% its thrusters for station-keeping about every three weeks; the last panel
% shows how a 1 cm/s velocity error grows when left uncorrected.

c = jwst_constants();
mu = c.mu;
L = c.AU / 1e3;                           % length unit [km]
Tunit = c.yearDays / (2*pi);              % time unit [days]
V = L / (Tunit * 86400);                  % velocity unit [km/s]

[xL2, gamma] = cr3bp_l2(mu);
fprintf('\n== Sun-Earth L2 halo orbit ==\n');
fprintf('L2 lies %.0f km beyond Earth\n', gamma * L);

AzKm = 400e3;
[S0, Tguess] = cr3bp_halo_guess(mu, AzKm / L, true);
[S0, T, info] = cr3bp_halo_correct(mu, S0);
[t, X] = cr3bp_propagate(mu, S0, T, 600);

rel = (X(:, 1:3) - [xL2 0 0]) * L;        % km relative to L2
earth = [(1 - mu - xL2) * L, 0, 0];
lambdaMax = max(abs(eig(info.monodromy)));
eFold = T * Tunit / log(lambdaMax);

fprintf('Differential correction converged in %d iterations (residual %.1e)\n', ...
        info.iterations, info.residual);
fprintf('Period:            %.1f days (first guess %.1f days)\n', T*Tunit, Tguess*Tunit);
fprintf('Extent along x:    %.0f km\n', max(rel(:,1)) - min(rel(:,1)));
fprintf('Extent along y:    %.0f km\n', max(rel(:,2)) - min(rel(:,2)));
fprintf('Extent along z:    %.0f km\n', max(rel(:,3)) - min(rel(:,3)));
fprintf('Closure error:     %.2f km after one revolution\n', norm(X(end,1:3) - X(1,1:3)) * L);
fprintf('Largest monodromy eigenvalue %.0f -> errors grow e-fold every %.1f days\n', ...
        lambdaMax, eFold);

%% Uncorrected perturbation growth
dv = 1e-5 / V;                            % 1 cm/s along y
nOut = 800;
[tp, Xn] = cr3bp_propagate(mu, S0, 2*T, nOut);
[~, Xp] = cr3bp_propagate(mu, S0 + [0 0 0 0 dv 0]', 2*T, nOut);
sepKm = sqrt(sum((Xp(:,1:3) - Xn(:,1:3)).^2, 2)) * L;

%% Figures
fig = new_figure('L2 halo orbit', [1300 820]);
blue = [0.20 0.45 0.85];
rk = rel / 1e3;                           % thousands of km, for readable axes
ek = earth / 1e3;

subplot(2, 2, 1);
plot3(rk(:,1), rk(:,2), rk(:,3), 'Color', blue, 'LineWidth', 1.6); hold on;
plot3(0, 0, 0, 'k+', 'MarkerSize', 10, 'LineWidth', 1.5);
plot3(ek(1), 0, 0, 'o', 'MarkerFaceColor', [0.1 0.5 0.9], 'MarkerEdgeColor', 'none');
hold off; grid on; axis equal; view(-35, 22);
xlabel('x [10^3 km]'); ylabel('y [10^3 km]'); zlabel('z [10^3 km]');
title(sprintf('Halo orbit around L2 (period %.0f days)', T*Tunit));
legend('JWST halo orbit', 'L2', 'Earth', 'Location', 'northeast');

projections = {[2 3], [1 2], [1 3]};
labels = {'y', 'z'; 'x', 'y'; 'x', 'z'};
titles = {'As seen from Earth (y-z)', 'Ecliptic view (x-y)', 'Side view (x-z)'};
for k = 1:2
    subplot(2, 2, k + 1);
    a = projections{k};
    plot(rk(:,a(1)), rk(:,a(2)), 'Color', blue, 'LineWidth', 1.6); hold on;
    plot(0, 0, 'k+', 'MarkerSize', 10, 'LineWidth', 1.5);
    if k == 2
        plot(ek(1), 0, 'o', 'MarkerFaceColor', [0.1 0.5 0.9], 'MarkerEdgeColor', 'none');
        text(ek(1), 0, '  Earth');
    end
    hold off; grid on; axis equal;
    xlabel([labels{k,1} ' [10^3 km]']); ylabel([labels{k,2} ' [10^3 km]']);
    title(titles{k});
end

subplot(2, 2, 4);
semilogy(tp * Tunit, max(sepKm, 1e-3), 'Color', [0.85 0.3 0.1], 'LineWidth', 1.6);
grid on;
xlabel('days since perturbation'); ylabel('deviation from halo [km]');
title('Unstable: a 1 cm/s error left uncorrected');

save_figure(fig, 'l2_halo_orbit');
