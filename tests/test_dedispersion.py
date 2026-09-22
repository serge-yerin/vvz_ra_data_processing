"""Dedispersion must stay bit-identical to the original memmap-column loop."""

import numpy as np

from dspz_pipeline.analysis.dedispersion import dedisperse, ind_search
from dspz_pipeline.config import DM_TOTAL_STEPS, HEAD_OFFSET, HEADER_SIZE_BYTES
from dspz_pipeline.io.dmt import read_dmt
from dspz_pipeline.io.jds_reader import JdsHeader, write_ucd_header


def _reference(ucd_data, shifts, picsize):
    """Original algorithm: ucd_data is (picsize, wofsg), columns are channels."""
    wofsg = ucd_data.shape[1]
    dedispersed = np.zeros(picsize, dtype=np.float64)
    for ch in range(wofsg):
        s = int(shifts[ch])
        if s >= picsize:
            continue
        dedispersed[:picsize - s] += ucd_data[s:, ch].astype(np.float64)
    return dedispersed.astype(np.float32)


def test_dedisperse_matches_reference():
    picsize, wofsg = 5000, 64
    rng = np.random.default_rng(0)
    data_tp = rng.standard_normal((picsize, wofsg)).astype(np.float32)   # file layout
    shifts = np.sort(rng.integers(0, 6000, size=wofsg))                 # some >= picsize
    out = dedisperse(np.ascontiguousarray(data_tp.T), shifts)
    np.testing.assert_array_equal(out, _reference(data_tp, shifts, picsize))
    assert out.dtype == np.float32


def _write_tiny_ucd(path, data_tp, wofsg):
    hdr = JdsHeader()
    hdr.sdspp[14 + HEAD_OFFSET] = wofsg
    hdr.sdspp[15 + HEAD_OFFSET] = 64            # avrs -> time resolution
    write_ucd_header(path, hdr, 1024)
    with open(path, "ab") as fh:
        fh.write(data_tp.astype(np.float32).tobytes(order="C"))


def test_ind_search_end_to_end(tmp_path):
    picsize, wofsg = 3000, 64
    rng = np.random.default_rng(1)
    data_tp = rng.standard_normal((picsize, wofsg)).astype(np.float32)
    ucd = tmp_path / "tiny.ucd"
    _write_tiny_ucd(ucd, data_tp, wofsg)

    dmt_path = ind_search(ucd, 12.872)
    acc, n_dm, n_t = read_dmt(dmt_path)
    assert (n_dm, n_t) == (DM_TOTAL_STEPS, picsize)
    assert np.isfinite(acc).all()
    # row 25 (central DM) must equal the reference computed with the same delays
    from dspz_pipeline.io.dmt import compute_dm_delays
    from dspz_pipeline.io.jds_reader import read_ucd_header
    hdr = read_ucd_header(ucd)
    dt = compute_dm_delays(12.872, 33.0, 16.5, wofsg, 16.5 / wofsg)
    shifts = np.round(dt / hdr.time_res_s).astype(np.int64)
    np.testing.assert_array_equal(acc[25], _reference(data_tp, shifts, picsize))
