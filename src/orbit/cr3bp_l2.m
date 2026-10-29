function [xL2, gamma] = cr3bp_l2(mu)
%CR3BP_L2 Location of the L2 Lagrange point.
%   [xL2, gamma] = CR3BP_L2(mu) returns the non-dimensional x coordinate of
%   L2 in the rotating frame and its distance gamma beyond the smaller
%   primary, from the quintic gamma^5 + (3-mu)gamma^4 + (3-2mu)gamma^3
%   - mu gamma^2 - 2mu gamma - mu = 0.

    g = (mu/3)^(1/3);
    for k = 1:50
        f  = g^5 + (3-mu)*g^4 + (3-2*mu)*g^3 - mu*g^2 - 2*mu*g - mu;
        df = 5*g^4 + 4*(3-mu)*g^3 + 3*(3-2*mu)*g^2 - 2*mu*g - 2*mu;
        step = f / df;
        g = g - step;
        if abs(step) < 1e-15
            break;
        end
    end
    gamma = g;
    xL2 = 1 - mu + gamma;
end
