function [t, S, Phi] = cr3bp_propagate(mu, S0, tf, nOut)
%CR3BP_PROPAGATE Integrate a CR3BP trajectory (and its STM) over [0, tf].
%   [t, S] = CR3BP_PROPAGATE(mu, S0, tf) returns times t and states S
%   (rows of [x y z vx vy vz]). With nOut, the solution is sampled at nOut
%   evenly spaced times. The third output is the 6x6 state transition
%   matrix from 0 to tf, which is only integrated when requested.

    S0 = S0(:);
    if nargin < 4 || isempty(nOut)
        span = [0 tf];
    else
        span = linspace(0, tf, nOut);
    end
    opt = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
    f = @(tt, s) cr3bp_eom(tt, s, mu);

    if nargout < 3
        [t, S] = ode45(f, span, S0, opt);
        return;
    end

    [t, X] = ode45(f, span, [S0; reshape(eye(6), 36, 1)], opt);
    S = X(:, 1:6);
    Phi = reshape(X(end, 7:42), 6, 6);
end
