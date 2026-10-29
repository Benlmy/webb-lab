function [S0, T, info] = cr3bp_halo_correct(mu, S0, opts)
%CR3BP_HALO_CORRECT Differential correction of a symmetric halo orbit.
%   [S0, T, info] = CR3BP_HALO_CORRECT(mu, S0) refines an initial state
%   [x0 0 z0 0 vy0 0] on the x-z plane so that the trajectory crosses the
%   plane y = 0 again perpendicularly (vx = vz = 0). By the mirror symmetry
%   of the CR3BP that crossing is the half period, so the orbit closes.
%   z0 is held fixed while x0 and vy0 are corrected.
%
%   opts (optional struct): tol (default 1e-11), maxIter (default 30).
%   info has fields iterations, residual and monodromy (6x6 state
%   transition matrix over one full period).

    if nargin < 3
        opts = struct();
    end
    if ~isfield(opts, 'tol'),     opts.tol = 1e-11; end
    if ~isfield(opts, 'maxIter'), opts.maxIter = 30; end

    S0 = S0(:);
    converged = false;
    for iter = 1:opts.maxIter
        [tHalf, Sh, PhiHalf] = half_period(mu, S0);
        residual = norm([Sh(4) Sh(6)]);
        if residual < opts.tol
            converged = true;
            break;
        end

        dS = cr3bp_eom(0, Sh, mu);
        % Sensitivity of (vx, vz) at the crossing to (x0, vy0), including
        % the shift in crossing time needed to stay on y = 0.
        A = [PhiHalf(4,1) PhiHalf(4,5); PhiHalf(6,1) PhiHalf(6,5)] ...
            - (1/Sh(5)) * [dS(4); dS(6)] * [PhiHalf(2,1) PhiHalf(2,5)];
        delta = -A \ [Sh(4); Sh(6)];
        S0(1) = S0(1) + delta(1);
        S0(5) = S0(5) + delta(2);
    end

    if ~converged
        error('cr3bp_halo_correct:noConvergence', ...
              'Differential correction did not converge (residual %.3g).', residual);
    end

    T = 2 * tHalf;
    [~, ~, Phi] = cr3bp_propagate(mu, S0, T);

    info.iterations = iter;
    info.residual = residual;
    info.monodromy = Phi;
end

function [tHalf, Sh, PhiHalf] = half_period(mu, S0)
    S = [S0; reshape(eye(6), 36, 1)];
    f = @(t, s) cr3bp_eom(t, s, mu);
    base = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
    % Step away from the starting point, which already lies on y = 0,
    % before arming the crossing event.
    [~, Sa] = ode45(f, [0 0.5], S, base);
    [~, ~, te, Se] = ode45(f, [0.5 10], Sa(end, :)', odeset(base, 'Events', @y_crossing));
    if isempty(te)
        error('cr3bp_halo_correct:noCrossing', 'Trajectory never returned to y = 0.');
    end
    tHalf = te(1);
    Sh = Se(1, 1:6)';
    PhiHalf = reshape(Se(1, 7:42), 6, 6);
end

function [value, isterminal, direction] = y_crossing(~, s)
    value = s(2);
    isterminal = 1;
    direction = 0;
end
