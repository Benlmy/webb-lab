function [S0, Tguess] = cr3bp_halo_guess(mu, Az, northern)
%CR3BP_HALO_GUESS Richardson third-order initial guess for an L2 halo orbit.
%   [S0, Tguess] = CR3BP_HALO_GUESS(mu, Az) returns an approximate initial
%   state [x y z vx vy vz] (non-dimensional, rotating frame) on the x-z
%   plane crossing of an L2 halo orbit with out-of-plane amplitude Az
%   (non-dimensional, i.e. a fraction of the Sun-Earth distance), and the
%   approximate period. Set northern = false for the southern family.
%
%   Reference: D. L. Richardson, "Analytic construction of periodic orbits
%   about the collinear points", Celestial Mechanics 22 (1980).

    if nargin < 3
        northern = true;
    end

    [xL2, gamma] = cr3bp_l2(mu);

    cn = @(n) (1/gamma^3) * ((-1)^n*mu + (-1)^n*(1-mu)*gamma^(n+1)/(1+gamma)^(n+1));
    c2 = cn(2); c3 = cn(3); c4 = cn(4);

    lam = sqrt((2 - c2 + sqrt((c2-2)^2 + 4*(c2-1)*(1+2*c2))) / 2);
    k = (lam^2 + 1 + 2*c2) / (2*lam);
    d1 = 3*lam^2/k * (k*(6*lam^2 - 1) - 2*lam);
    d2 = 8*lam^2/k * (k*(11*lam^2 - 1) - 2*lam);

    a21 = 3*c3*(k^2 - 2) / (4*(1 + 2*c2));
    a22 = 3*c3 / (4*(1 + 2*c2));
    a23 = -3*c3*lam/(4*k*d1) * (3*k^3*lam - 6*k*(k - lam) + 4);
    a24 = -3*c3*lam/(4*k*d1) * (2 + 3*k*lam);
    b21 = -3*c3*lam/(2*d1) * (3*k*lam - 4);
    b22 = 3*c3*lam/d1;
    d21 = -c3/(2*lam^2);

    a31 = -9*lam/(4*d2) * (4*c3*(k*a23 - b21) + k*c4*(4 + k^2)) ...
          + (9*lam^2 + 1 - c2)/(2*d2) * (3*c3*(2*a23 - k*b21) + c4*(2 + 3*k^2));
    a32 = -1/d2 * (9*lam/4*(4*c3*(k*a24 - b22) + k*c4) ...
          + 3/2*(9*lam^2 + 1 - c2)*(c3*(k*b22 + d21 - 2*a24) - c4));
    b31 = 3/(8*d2) * (8*lam*(3*c3*(k*b21 - 2*a23) - c4*(2 + 3*k^2)) ...
          + (9*lam^2 + 1 + 2*c2)*(4*c3*(k*a23 - b21) + k*c4*(4 + k^2)));
    b32 = 1/d2 * (9*lam*(c3*(k*b22 + d21 - 2*a24) - c4) ...
          + 3/8*(9*lam^2 + 1 + 2*c2)*(4*c3*(k*a24 - b22) + k*c4));
    d31 = 3/(64*lam^2) * (4*c3*a24 + c4);
    d32 = 3/(64*lam^2) * (4*c3*(a23 - d21) + c4*(4 + k^2));

    den = 2*lam*(lam*(1 + k^2) - 2*k);
    s1 = 1/den * (3/2*c3*(2*a21*(k^2 - 2) - a23*(k^2 + 2) - 2*k*b21) ...
         - 3/8*c4*(3*k^4 - 8*k^2 + 8));
    s2 = 1/den * (3/2*c3*(2*a22*(k^2 - 2) + a24*(k^2 + 2) + 2*k*b22 + 5*d21) ...
         + 3/8*c4*(12 - k^2));
    l1 = -3/2*c3*(2*a21 + a23 + 5*d21) - 3/8*c4*(12 - k^2) + 2*lam^2*s1;
    l2 = 3/2*c3*(a24 - 2*a22) + 9/8*c4 + 2*lam^2*s2;
    Delta = lam^2 - c2;

    az = Az / gamma;
    ax = sqrt((-Delta - l2*az^2) / l1);
    omega = 1 + s1*ax^2 + s2*az^2;
    dm = 1;
    if ~northern
        dm = -1;
    end

    tau = 0;
    x = a21*ax^2 + a22*az^2 - ax*cos(tau) + (a23*ax^2 - a24*az^2)*cos(2*tau) ...
        + (a31*ax^3 - a32*ax*az^2)*cos(3*tau);
    z = dm*az*cos(tau) + dm*d21*ax*az*(cos(2*tau) - 3) ...
        + dm*(d32*az*ax^2 - d31*az^3)*cos(3*tau);
    ydot = lam*omega*(k*ax*cos(tau) + 2*(b21*ax^2 - b22*az^2)*cos(2*tau) ...
        + 3*(b31*ax^3 - b32*ax*az^2)*cos(3*tau));

    % Richardson's coordinates are centred on L2 and scaled by gamma.
    S0 = [xL2 + gamma*x; 0; gamma*z; 0; gamma*ydot; 0];
    Tguess = 2*pi / (lam*omega);
end
