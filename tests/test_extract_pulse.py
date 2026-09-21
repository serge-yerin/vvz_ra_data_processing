"""extract_pulse must reproduce the original per-channel shift loop exactly."""

import numpy as np

from dspz_pipeline.gui.show_pulse import extract_pulse
from dspz_pipeline.io.dmt import compute_dm_delays


def _reference(dat_ucd, dt, time_res, sh_t, nsf):
    wofsg = dat_ucd.shape[0]
    pulse = np.zeros((wofsg, 2 * nsf + 1), dtype=np.float64)
    for j in range(wofsg):
        stsp = int(round(dt[j] / time_res)) + sh_t + nsf
        stsp = max(0, min(stsp, dat_ucd.shape[1] - 2 * nsf - 1))
        end = stsp + 2 * nsf + 1
        if end <= dat_ucd.shape[1]:
            pulse[j, :] = dat_ucd[j, stsp:end]
    return pulse


def test_matches_reference_for_real_delays():
    wofsg, n_time, nsf = 4096, 3000, 50
    time_res = 64 * 8192.0 / 66_000_000.0
    dat = np.random.default_rng(0).standard_normal((wofsg, n_time)).astype(np.float32)
    dt = compute_dm_delays(12.872, 33.0, 16.5, wofsg, 16.5 / wofsg)
    for sh_t in (-49, -5, 0, 3, 49):
        np.testing.assert_array_equal(
            extract_pulse(dat, dt, time_res, sh_t, nsf),
            _reference(dat, dt, time_res, sh_t, nsf),
        )


def test_clamps_when_delay_exceeds_window():
    wofsg, n_time, nsf = 64, 150, 50          # delays far beyond the array
    dat = np.random.default_rng(1).standard_normal((wofsg, n_time)).astype(np.float32)
    dt = compute_dm_delays(200.0, 33.0, 16.5, wofsg, 16.5 / wofsg)
    np.testing.assert_array_equal(
        extract_pulse(dat, dt, 0.00794, 0, nsf),
        _reference(dat, dt, 0.00794, 0, nsf),
    )
