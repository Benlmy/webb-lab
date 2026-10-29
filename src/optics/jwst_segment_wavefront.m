function opd = jwst_segment_wavefront(P, piston, tipX, tiltY)
%JWST_SEGMENT_WAVEFRONT Optical path difference from rigid-body segment errors.
%   opd = JWST_SEGMENT_WAVEFRONT(P, piston, tipX, tiltY) builds the OPD map
%   [m] on the pupil grid P (see JWST_PUPIL). Each input is an 18-element
%   vector with one value per segment:
%     piston  OPD offset of the whole segment [m]
%     tipX    OPD slope along x across the segment [m/m = rad]
%     tiltY   OPD slope along y across the segment [m/m = rad]
%   A slope of s radians moves that segment's image by s radians on the sky.
%   Missing inputs default to zero.

    n = size(P.seg.center, 1);
    if nargin < 2 || isempty(piston), piston = zeros(n, 1); end
    if nargin < 3 || isempty(tipX),   tipX   = zeros(n, 1); end
    if nargin < 4 || isempty(tiltY),  tiltY  = zeros(n, 1); end

    opd = zeros(size(P.mask));
    for k = 1:n
        in = P.segmentId == k;
        cx = P.seg.center(k, 1);
        cy = P.seg.center(k, 2);
        opd(in) = piston(k) + tipX(k) * (P.x(in) - cx) + tiltY(k) * (P.y(in) - cy);
    end
end
