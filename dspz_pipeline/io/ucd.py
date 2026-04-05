"""
Shared reader for .ucd observation files (cleaned DSPZ data).

The .ucd format has the same 1024-byte header as .jds files, followed
by float32 spectrograms in Fortran (column-major) order.

This module provides functions used by both stage 1 and stage 2
(indsearch) of the pipeline.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np

from dspz_pipeline.config import (
    HEADER_SIZE_BYTES,
    HEAD_OFFSET,
    NYQUIST_BANDWIDTH_MHZ,
    TOTAL_FFT_CHANNELS,
    FFT_LENGTH,
    CLOCK_FREQ_HZ,
    SDSPP_COUNT,
    HEADER_TEXT_FIELDS,
)
from dspz_pipeline.utils import decode_text


def read_ucd_header_as_dict(path: str | Path) -> dict:
    """Read and parse a .ucd file header, returning a dict.

    This is the common interface used by indsearch and other tools
    that need header parameters as a plain dict.

    Parameters
    ----------
    path : str or Path
        Path to the .ucd file.

    Returns
    -------
    dict with keys:
        sname, stime, sgmtt, ssysn, splace, sdesc  — decoded strings
        sdspp     — raw uint32 array, length 128
        nofs      — spectra per frame
        fmin, fmax — frequency range [MHz]
        wofsg     — number of frequency channels
        avrs      — number of averages
        time_res  — time resolution [seconds]
    """
    with open(path, "rb") as f:
        sname = f.read(32)
        stime = f.read(32)
        sgmtt = f.read(32)
        ssysn = f.read(32)
        ssyst = f.read(32)
        splace = f.read(96)
        sdesc = f.read(256)
        sdspp = np.frombuffer(f.read(SDSPP_COUNT * 4), dtype="<u4")

    nofs = int(sdspp[10 + HEAD_OFFSET])
    fmin = float(sdspp[12 + HEAD_OFFSET]) * NYQUIST_BANDWIDTH_MHZ / TOTAL_FFT_CHANNELS
    fmax = float(sdspp[13 + HEAD_OFFSET]) * NYQUIST_BANDWIDTH_MHZ / TOTAL_FFT_CHANNELS
    wofsg = int(sdspp[14 + HEAD_OFFSET])
    avrs = int(sdspp[15 + HEAD_OFFSET])
    time_res = avrs * FFT_LENGTH / CLOCK_FREQ_HZ

    return {
        "sname": decode_text(sname),
        "stime": decode_text(stime),
        "sgmtt": decode_text(sgmtt),
        "ssysn": decode_text(ssysn),
        "splace": decode_text(splace),
        "sdesc": decode_text(sdesc),
        "sdspp": sdspp,
        "nofs": nofs,
        "fmin": fmin,
        "fmax": fmax,
        "wofsg": wofsg,
        "avrs": avrs,
        "time_res": time_res,
    }


def count_frames(file_path: str | Path, wofsg: int, nofs: int) -> int:
    """Return the number of complete data frames in a .ucd file.

    Replicates IDL::

        n_kadr = floor((fstatus.size - 1024) / (4L * nofs * wofsg))
    """
    file_size = os.path.getsize(file_path)
    return int((file_size - HEADER_SIZE_BYTES) // (4 * nofs * wofsg))


def read_frame(
    f,
    frame_idx: int,
    wofsg: int,
    nofs: int,
) -> np.ndarray:
    """Read one float32 data frame from an open .ucd file.

    Replicates IDL::

        datdspz = assoc(1, fltarr(wofsg, nofs), 1024)
        data    = datdspz[j]

    Parameters
    ----------
    f : binary file object
    frame_idx : int
        Zero-based frame index.
    wofsg : int
        Number of frequency channels.
    nofs : int
        Number of time samples per frame.

    Returns
    -------
    ndarray, shape (wofsg, nofs), dtype float32
    """
    offset = HEADER_SIZE_BYTES + frame_idx * wofsg * nofs * 4
    f.seek(offset)
    raw = np.frombuffer(f.read(wofsg * nofs * 4), dtype="<f4")
    return raw.reshape((wofsg, nofs), order="F")
