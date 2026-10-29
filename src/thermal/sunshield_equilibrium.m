function res = sunshield_equilibrium(p)
%SUNSHIELD_EQUILIBRIUM Steady-state temperatures of a multi-layer sunshield.
%   res = SUNSHIELD_EQUILIBRIUM(p) solves the radiative energy balance of n
%   parallel membranes per unit area. Layer 1 faces the Sun. Each pair of
%   adjacent layers exchanges heat as grey parallel plates, but only a
%   fraction p.coupling of the radiation leaving an inner face reaches the
%   next layer: the rest escapes to space through the open, V-shaped edges
%   between layers. That edge venting is what makes the JWST design work.
%
%   Any field of p may be omitted (defaults in brackets):
%     nLayers        number of membranes                       [5]
%     solarFlux      irradiance on the sun-facing layer W/m^2   [1361/1.01^2]
%     sunAngleDeg    angle between Sun line and layer normal    [0]
%     alphaSun       solar absorptance of layer 1 front face    [0.55]
%     epsHot         IR emissivity of layer 1 front face        [0.80]
%     epsInner       emissivity of the inner aluminised faces   [0.03]
%     epsCold        emissivity of the last layer's back face   [0.03]
%     coupling       fraction of inner-face radiation reaching
%                    the neighbouring layer                     [0.34]
%
%   res has fields T (n x 1 temperatures [K]), heatLeak (W/m^2 radiated
%   from the cold side), absorbed (W/m^2) and params.

    if nargin < 1
        p = struct();
    end
    def = struct('nLayers', 5, 'solarFlux', 1361/1.01^2, 'sunAngleDeg', 0, ...
                 'alphaSun', 0.55, 'epsHot', 0.80, 'epsInner', 0.03, ...
                 'epsCold', 0.03, 'coupling', 0.34);
    names = fieldnames(def);
    for k = 1:numel(names)
        if ~isfield(p, names{k})
            p.(names{k}) = def.(names{k});
        end
    end

    sigma = 5.670374419e-8;
    n = p.nLayers;
    absorbed = p.alphaSun * p.solarFlux * max(cosd(p.sunAngleDeg), 0);

    epsF = [p.epsHot; repmat(p.epsInner, n-1, 1)];   % sun-facing side
    epsB = [repmat(p.epsInner, n-1, 1); p.epsCold];  % space-facing side
    if n == 1
        epsB = p.epsCold;
    end

    % Effective emissivity between back of layer i and front of layer i+1.
    epsPair = 1 ./ (1 ./ epsB(1:n-1) + 1 ./ epsF(2:n) - 1);

    % Face emissivities to deep space: outer faces radiate fully, inner
    % faces only through the open edges.
    vent = 1 - p.coupling;
    emitF = vent * epsF;  emitF(1) = epsF(1);
    emitB = vent * epsB;  emitB(n) = epsB(n);

    % The balance is linear in u = T^4: K u = q.
    K = diag(emitF + emitB);
    for i = 1:n-1
        g = p.coupling * epsPair(i);
        K(i, i) = K(i, i) + g;
        K(i+1, i+1) = K(i+1, i+1) + g;
        K(i, i+1) = K(i, i+1) - g;
        K(i+1, i) = K(i+1, i) - g;
    end
    q = [absorbed; zeros(n-1, 1)] / sigma;
    T = (K \ q).^0.25;

    res.T = T;
    res.heatLeak = sigma * emitB(n) * T(n)^4;
    res.absorbed = absorbed;
    res.params = p;
end
