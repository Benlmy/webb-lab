function P = jwst_pupil(N, c)
%JWST_PUPIL Sampled entrance pupil of JWST.
%   P = JWST_PUPIL(N) samples the pupil on an N-by-N grid spanning 6.6 m.
%   The pupil includes the 18 hexagonal segments with their gaps, the
%   secondary mirror obscuration and the three secondary support struts.
%
%   Fields of P:
%     mask      NxN logical transmission map
%     segmentId NxN map: k inside segment k, 0 elsewhere (struts included)
%     x, y      NxN pupil coordinates [m]
%     dx        grid spacing [m]
%     area      clear collecting area [m^2]
%     seg       output of JWST_SEGMENT_LAYOUT

    if nargin < 1 || isempty(N)
        N = 512;
    end
    if nargin < 2
        c = jwst_constants();
    end

    seg = jwst_segment_layout(c);
    span = 6.6;
    dx = span / N;
    v = ((0:N-1) - (N-1)/2) * dx;
    [x, y] = meshgrid(v, v);

    segmentId = zeros(N, N);
    R = seg.circumR;
    a = seg.apothem;
    for k = 1:size(seg.center, 1)
        ax = abs(x - seg.center(k, 1));
        ay = abs(y - seg.center(k, 2));
        inside = (ay <= a) & (sqrt(3)*ax + ay <= sqrt(3)*R);
        segmentId(inside) = k;
    end

    blocked = hypot(x, y) <= c.secondaryDiameter/2;
    for ang = c.strutAnglesDeg
        u = [cosd(ang) sind(ang)];
        along = x*u(1) + y*u(2);
        across = abs(-x*u(2) + y*u(1));
        blocked = blocked | (along > 0 & across <= c.strutWidth/2);
    end

    mask = segmentId > 0 & ~blocked;
    segmentId(~mask) = 0;

    P.mask = mask;
    P.segmentId = segmentId;
    P.x = x;
    P.y = y;
    P.dx = dx;
    P.area = nnz(mask) * dx^2;
    P.seg = seg;
end
