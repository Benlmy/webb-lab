function [dropout, detect, lambdaObs] = lyman_break_filter(z)
%LYMAN_BREAK_FILTER NIRCam filters bracketing the Lyman break of a galaxy.
%   [dropout, detect, lambdaObs] = LYMAN_BREAK_FILTER(z) redshifts the
%   Lyman-alpha break (121.6 nm) to lambdaObs [um]. Neutral hydrogen
%   absorbs nearly all light blueward of the break, so the galaxy vanishes
%   ("drops out") in filters below it. dropout is the reddest filter lying
%   entirely blueward of the break and detect is the bluest filter lying
%   entirely redward of it. Either is '' when no NIRCam filter qualifies.

    lambdaObs = 0.12157 * (1 + z);
    F = nircam_filters();
    hi = [F.hi];
    lo = [F.lo];

    below = find(hi < lambdaObs, 1, 'last');
    above = find(lo > lambdaObs, 1, 'first');

    dropout = '';
    detect = '';
    if ~isempty(below), dropout = F(below).name; end
    if ~isempty(above), detect = F(above).name; end
end
