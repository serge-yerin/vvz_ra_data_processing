"""erov and the DSPZ frame decoder must stay bit-identical to the originals."""

import numpy as np
import pytest

from dspz_pipeline.cleaning.robust_stats import erov
from dspz_pipeline.config import DECODE_SCALE, EXPONENT_MASK, MANTISSA_MASK
from dspz_pipeline.io.jds_reader import decode_dspz_frame


def _erov_reference(data):
    m = np.array(data, dtype=np.float64).ravel().copy()
    if m.size == 0:
        return 0.0, 0.0
    if m.size == 1:
        return float(m[0]), 0.0
    sr = np.mean(m)
    count = m.size
    while True:
        srkv = np.std(m, ddof=1) if m.size > 1 else 0.0
        count_pr = count
        sr_pr = sr
        mask = np.abs(m - sr) <= srkv * 3.0
        m = m * mask
        m_in = np.nonzero(m)[0]
        count = m_in.size
        if count == 0:
            return 0.0, 0.0
        m = m[m_in]
        sr = np.mean(m)
        ster = np.std(m, ddof=1) if m.size > 1 else 0.0
        if sr == 0.0:
            break
        if abs(sr_pr / sr - 1.0) < 1e-5 or count_pr == count:
            break
    return sr, ster


@pytest.mark.parametrize("seed", range(20))
def test_erov_matches_reference(seed):
    rng = np.random.default_rng(seed)
    n = int(rng.integers(2, 3000))
    data = rng.standard_normal(n) * rng.uniform(0.1, 50) + rng.uniform(-10, 10)
    data[rng.integers(0, n, size=n // 20)] *= 50          # outliers
    data[rng.integers(0, n, size=n // 50)] = 0.0          # exact zeros (IDL where() quirk)
    assert erov(data) == _erov_reference(data)


@pytest.mark.parametrize("data", [np.array([]), np.array([3.5]), np.zeros(10), np.ones(10),
                                  np.array([1.0, 1.0, 1.0, 1e9])])
def test_erov_edge_cases(data):
    assert erov(data) == _erov_reference(data)


def _decode_reference(raw, wofsg, nofs, avrs, mode):
    mant = (raw & MANTISSA_MASK).astype(np.float64)
    expn = (raw & EXPONENT_MASK).astype(np.int32)
    data = mant / np.power(2.0, expn) * DECODE_SCALE / avrs
    if mode == 0:
        imdat = data[0, :, :] - data[1, :, :]
    elif mode == 1:
        imdat = data[0, :, :]
    else:
        imdat = data[1, :, :]
    return imdat.reshape(wofsg, nofs)


@pytest.mark.parametrize("mode", [0, 1, 2])
def test_decode_matches_reference(mode):
    wofsg, nofs, avrs = 128, 256, 64
    rng = np.random.default_rng(mode)
    raw_file_layout = rng.integers(0, 2**32, size=(nofs, wofsg, 2), dtype=np.uint32)
    raw = raw_file_layout.transpose(2, 1, 0)           # non-contiguous view, as the reader passes it
    out = decode_dspz_frame(raw, wofsg, nofs, avrs, mode)
    ref = _decode_reference(raw.copy(), wofsg, nofs, avrs, mode)
    np.testing.assert_array_equal(out, ref)
    assert out.flags.c_contiguous
