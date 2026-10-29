function dS = cr3bp_eom(~, S, mu)
%CR3BP_EOM Equations of motion of the circular restricted three-body problem.
%   dS = CR3BP_EOM(t, S, mu) returns the time derivative of the state in the
%   rotating, non-dimensional Sun-Earth frame (Sun at x = -mu, Earth at
%   x = 1-mu). S is either a 6-element state [x y z vx vy vz] or a
%   42-element vector that also carries the 6x6 state transition matrix
%   (column-major) after the state.

    x = S(1); y = S(2); z = S(3);
    vx = S(4); vy = S(5); vz = S(6);

    r1 = sqrt((x + mu)^2 + y^2 + z^2);
    r2 = sqrt((x - 1 + mu)^2 + y^2 + z^2);
    r13 = r1^3; r23 = r2^3;

    ax = 2*vy + x - (1-mu)*(x+mu)/r13 - mu*(x-1+mu)/r23;
    ay = -2*vx + y - (1-mu)*y/r13 - mu*y/r23;
    az = -(1-mu)*z/r13 - mu*z/r23;

    if numel(S) == 6
        dS = [vx; vy; vz; ax; ay; az];
        return;
    end

    r15 = r1^5; r25 = r2^5;
    Uxx = 1 - (1-mu)/r13 - mu/r23 + 3*(1-mu)*(x+mu)^2/r15 + 3*mu*(x-1+mu)^2/r25;
    Uyy = 1 - (1-mu)/r13 - mu/r23 + 3*(1-mu)*y^2/r15 + 3*mu*y^2/r25;
    Uzz =   - (1-mu)/r13 - mu/r23 + 3*(1-mu)*z^2/r15 + 3*mu*z^2/r25;
    Uxy = 3*(1-mu)*(x+mu)*y/r15 + 3*mu*(x-1+mu)*y/r25;
    Uxz = 3*(1-mu)*(x+mu)*z/r15 + 3*mu*(x-1+mu)*z/r25;
    Uyz = 3*(1-mu)*y*z/r15 + 3*mu*y*z/r25;

    A = [zeros(3) eye(3);
         Uxx Uxy Uxz  0 2 0;
         Uxy Uyy Uyz -2 0 0;
         Uxz Uyz Uzz  0 0 0];

    Phi = reshape(S(7:42), 6, 6);
    dPhi = A * Phi;
    dS = [vx; vy; vz; ax; ay; az; dPhi(:)];
end
