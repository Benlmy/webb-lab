function ok = run_tests()
%RUN_TESTS Self-checking test suite for the JWST project.
%   ok = RUN_TESTS() runs every check, prints a pass/fail line for each and
%   returns true when all pass. Plain assertions keep it portable between
%   MATLAB and GNU Octave.

    root = fileparts(fileparts(mfilename('fullpath')));
    run(fullfile(root, 'startup.m'));

    tests = {
        @test_segment_layout
        @test_pupil_area
        @test_perfect_psf
        @test_mft_matches_fft
        @test_marechal
        @test_tilt_moves_image
        @test_l2_distance
        @test_halo_orbit_closes
        @test_sunshield_gradient
        @test_sunshield_energy_balance
        @test_planck_wien
        @test_lyman_break
    };

    nFail = 0;
    for k = 1:numel(tests)
        name = func2str(tests{k});
        try
            tests{k}();
            fprintf('  PASS  %s\n', name);
        catch err
            nFail = nFail + 1;
            fprintf('  FAIL  %s: %s\n', name, err.message);
        end
    end
    fprintf('%d of %d tests passed\n', numel(tests) - nFail, numel(tests));
    ok = nFail == 0;
end

function test_segment_layout()
    seg = jwst_segment_layout();
    assert(size(seg.center, 1) == 18, 'expected 18 segments');
    assert(sum(seg.ring == 1) == 6 && sum(seg.ring == 2) == 12, 'ring sizes wrong');
    assert(numel(unique(seg.label)) == 18, 'segment labels must be unique');
    nB = sum(strncmp(seg.label, 'B', 1));
    nC = sum(strncmp(seg.label, 'C', 1));
    assert(nB == 6 && nC == 6, 'expected six B and six C segments');
end

function test_pupil_area()
    P = jwst_pupil(512);
    % JWST quotes ~25.4 m^2; ideal sharp-edged hexagons give slightly more.
    assert(P.area > 25 && P.area < 27.2, sprintf('area %.2f m^2 out of range', P.area));
end

function test_perfect_psf()
    P = jwst_pupil(128);
    out = jwst_psf(P, 2e-6);
    assert(abs(out.strehl - 1) < 1e-10, 'unaberrated Strehl must be 1');
    [~, idx] = max(out.psf(:));
    [r, cidx] = ind2sub(size(out.psf), idx);
    assert(abs(out.theta(r)) < 1e-12 && abs(out.theta(cidx)) < 1e-12, 'peak not at origin');
end

function test_mft_matches_fft()
    P = jwst_pupil(128);
    rng(3);
    opd = jwst_segment_wavefront(P, 80e-9 * randn(18, 1));
    a = jwst_psf(P, 2e-6, opd, 4);
    mid = floor(size(a.psf, 1)/2) + 1;
    half = 16;
    idx = mid-half:mid+half-1;
    b = jwst_psf_mft(P, 2e-6, opd, 2*half*a.pixelScale, 2*half);
    err = max(max(abs(a.psf(idx, idx) - b.psf)));
    assert(err < 1e-9, sprintf('MFT and FFT disagree by %.2g', err));
end

function test_marechal()
    P = jwst_pupil(256);
    rng(1);
    opd = jwst_segment_wavefront(P, 60e-9 * randn(18, 1));
    out = jwst_psf(P, 2e-6, opd);
    predicted = exp(-(2*pi*out.rmsWfe/2e-6)^2);
    assert(abs(out.strehl - predicted) < 0.03, ...
           sprintf('Strehl %.3f vs Marechal %.3f', out.strehl, predicted));
end

function test_tilt_moves_image()
    P = jwst_pupil(256);
    shiftArcsec = 0.5;
    s = shiftArcsec / (180/pi*3600);
    % Segment tilts act about each segment's centre, so a global tilt also
    % needs the matching per-segment pistons.
    opd = jwst_segment_wavefront(P, s * P.seg.center(:, 1), s * ones(18, 1), []);
    out = jwst_psf(P, 2e-6, opd, 4);
    [~, idx] = max(out.psf(:));
    [~, col] = ind2sub(size(out.psf), idx);
    assert(abs(out.theta(col) - shiftArcsec) <= out.pixelScale, ...
           sprintf('image moved %.3f", expected %.3f"', out.theta(col), shiftArcsec));
end

function test_l2_distance()
    c = jwst_constants();
    [~, gamma] = cr3bp_l2(c.mu);
    km = gamma * c.AU / 1e3;
    assert(km > 1.49e6 && km < 1.52e6, sprintf('L2 distance %.0f km', km));
end

function test_halo_orbit_closes()
    c = jwst_constants();
    S0 = cr3bp_halo_guess(c.mu, 200e6 / c.AU);
    [S0, T] = cr3bp_halo_correct(c.mu, S0);
    [~, X] = cr3bp_propagate(c.mu, S0, T);
    closureKm = norm(X(end, 1:3) - S0(1:3)') * c.AU / 1e3;
    days = T * c.yearDays / (2*pi);
    assert(closureKm < 50, sprintf('orbit misses closure by %.1f km', closureKm));
    assert(days > 170 && days < 190, sprintf('period %.1f days', days));
end

function test_sunshield_gradient()
    r = sunshield_equilibrium();
    assert(all(diff(r.T) < 0), 'temperatures must fall layer by layer');
    assert(r.T(1) > 300 && r.T(1) < 400, sprintf('hot side %.1f K', r.T(1)));
    assert(r.T(end) < 50, sprintf('cold side %.1f K', r.T(end)));
    r6 = sunshield_equilibrium(struct('nLayers', 6));
    assert(r6.T(end) < r.T(end), 'more layers should run colder');
end

function test_sunshield_energy_balance()
    r = sunshield_equilibrium();
    p = r.params;
    sigma = 5.670374419e-8;
    vent = 1 - p.coupling;
    n = p.nLayers;
    epsF = [p.epsHot; repmat(p.epsInner, n-1, 1)];
    epsB = [repmat(p.epsInner, n-1, 1); p.epsCold];
    emit = vent * (epsF + epsB);
    emit(1) = epsF(1) + vent * epsB(1);
    emit(n) = vent * epsF(n) + epsB(n);
    radiated = sum(sigma * emit .* r.T.^4);
    assert(abs(radiated - r.absorbed) / r.absorbed < 1e-9, 'energy is not conserved');
end

function test_planck_wien()
    lam = linspace(1e-6, 30e-6, 20000);
    [~, i] = max(planck_radiance(lam, 288));
    wien = 2.897771955e-3 / 288;
    assert(abs(lam(i) - wien) / wien < 1e-3, 'peak disagrees with Wien''s law');
end

function test_lyman_break()
    [drop, det] = lyman_break_filter(12);
    assert(strcmp(drop, 'F115W') && strcmp(det, 'F200W'), ...
           sprintf('z=12 gave %s / %s', drop, det));
    [drop, det] = lyman_break_filter(7.5);
    assert(strcmp(drop, 'F090W') && strcmp(det, 'F150W'), ...
           sprintf('z=7.5 gave %s / %s', drop, det));
end
