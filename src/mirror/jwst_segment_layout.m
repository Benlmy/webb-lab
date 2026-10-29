function seg = jwst_segment_layout(c)
%JWST_SEGMENT_LAYOUT Centres and labels of the 18 primary mirror segments.
%   seg = JWST_SEGMENT_LAYOUT() lays out two rings of flat-topped hexagons
%   around an empty central position, matching the JWST primary mirror.
%
%   Fields of seg:
%     center   18x2 segment centres [m]
%     circumR  circumradius of one segment (centre to vertex) [m]
%     apothem  centre-to-flat distance of one segment [m]
%     label    18x1 cell of names (A1-A6 inner ring, B1-B6 and C1-C6 outer ring)
%     ring     18x1 ring index (1 or 2)

    if nargin < 1
        c = jwst_constants();
    end

    pitch = c.segmentFlatToFlat + c.segmentGap;  % centre-to-centre of neighbours
    Rgrid = pitch / sqrt(3);                     % circumradius of the lattice cell

    centers = zeros(0, 2);
    rings = zeros(0, 1);
    angles = zeros(0, 1);
    for q = -2:2
        for r = -2:2
            s = -q - r;
            d = max(abs([q r s]));
            if d == 0 || d > 2
                continue;
            end
            x = 1.5 * Rgrid * q;
            y = sqrt(3) * Rgrid * (r + q/2);
            centers(end+1, :) = [x y]; %#ok<AGROW>
            rings(end+1, 1) = d;       %#ok<AGROW>
            angles(end+1, 1) = mod(atan2(x, y), 2*pi); %#ok<AGROW> clockwise from +y
        end
    end

    [~, order] = sortrows([rings angles]);
    centers = centers(order, :);
    rings = rings(order);

    % JWST naming: A = inner ring, B = outer ring "corner" segments that sit
    % radially beyond an A segment, C = outer ring segments between B's.
    label = cell(18, 1);
    nA = 0; nB = 0; nC = 0;
    for k = 1:numel(rings)
        if rings(k) == 1
            nA = nA + 1;
            label{k} = sprintf('A%d', nA);
        else
            dist = norm(centers(k, :));
            if abs(dist - 2*pitch) < 1e-6
                nB = nB + 1;
                label{k} = sprintf('B%d', nB);
            else
                nC = nC + 1;
                label{k} = sprintf('C%d', nC);
            end
        end
    end

    seg.center  = centers;
    seg.circumR = c.segmentFlatToFlat / sqrt(3);
    seg.apothem = c.segmentFlatToFlat / 2;
    seg.label   = label;
    seg.ring    = rings;
end
