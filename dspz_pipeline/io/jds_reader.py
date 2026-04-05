"""
Reader for DSPZ .jds binary data files.

Translated from the file-reading portion of IDL ``process_survey.pro``.

The .jds format consists of:
  - A 1024-byte header containing text metadata and a 128-element UINT32
    parameter array (``Sdspp``).
  - A sequence of data frames, each containing ``UINT32[2, wofsg, nofs]``
    values packed in a custom floating-point encoding (5-bit exponent in
    the low bits, 27-bit mantissa in the high bits).

This module provides:
  - :func:`read_jds_header` — parse the header and return a dataclass.
  - :func:`decode_dspz_frame` — decode one raw UINT32 frame to float64.
  - :class:`JdsFile` — a context-managed reader that yields decoded frames.
  - :func:`write_ucd_header` — write a modified .jds header to a .ucd file.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator

import numpy as np

from dspz_pipeline.config import (
    CLOCK_FREQ_HZ,
    DECODE_SCALE,
    EXPONENT_MASK,
    FFT_LENGTH,
    HEAD_OFFSET,
    HEADER_SIZE_BYTES,
    HEADER_TEXT_FIELDS,
    MANTISSA_MASK,
    NYQUIST_BANDWIDTH_MHZ,
    SDSPP_COUNT,
    TOTAL_FFT_CHANNELS,
)
from dspz_pipeline.utils import decode_text


@dataclass
class JdsHeader:
    """Parsed content of a .jds file header."""

    sname: str = ""
    stime: str = ""
    sgmtt: str = ""
    ssysn: str = ""
    ssyst: str = ""
    splace: str = ""
    sdesc: str = ""
    sdspp: np.ndarray = field(default_factory=lambda: np.zeros(SDSPP_COUNT, dtype=np.uint32))

    # Derived convenience fields (populated by read_jds_header)
    data_mode: int = 0          # 0=waveform, 1=spectra, 2=correlation
    fmin_mhz: float = 0.0       # low edge of observing band (MHz)
    fmax_mhz: float = 0.0       # high edge of observing band (MHz)
    wofsg: int = 0              # number of frequency channels
    avrs: int = 0               # number of FFT averages per spectrum
    time_res_s: float = 0.0     # time resolution (seconds per spectrum)
    nofs: int = 0               # spectra per frame (stored in header or default)


def read_jds_header(path: str | Path) -> JdsHeader:
    """Read and parse a .jds file header.

    Parameters
    ----------
    path : str or Path
        Path to the .jds file.

    Returns
    -------
    JdsHeader
        Populated header dataclass.
    """
    hdr = JdsHeader()
    with open(path, "rb") as fh:
        raw = fh.read(HEADER_SIZE_BYTES)

    offset = 0
    for fname, nbytes in HEADER_TEXT_FIELDS.items():
        setattr(hdr, fname, decode_text(raw[offset: offset + nbytes]))
        offset += nbytes

    # Sdspp: 128 × UINT32 = 512 bytes
    hdr.sdspp = np.frombuffer(raw[offset: offset + SDSPP_COUNT * 4], dtype=np.uint32).copy()

    # Derive instrument parameters using HeadOffset=16
    ho = HEAD_OFFSET
    hdr.data_mode = int(hdr.sdspp[8 + ho])
    fmin_chan = int(hdr.sdspp[12 + ho])
    fmax_chan = int(hdr.sdspp[13 + ho])
    hdr.fmin_mhz = fmin_chan * NYQUIST_BANDWIDTH_MHZ / TOTAL_FFT_CHANNELS
    hdr.fmax_mhz = fmax_chan * NYQUIST_BANDWIDTH_MHZ / TOTAL_FFT_CHANNELS
    hdr.wofsg = int(hdr.sdspp[14 + ho])
    hdr.avrs = int(hdr.sdspp[15 + ho])
    hdr.time_res_s = hdr.avrs * FFT_LENGTH / CLOCK_FREQ_HZ
    # nofs may be stored at sdspp[10+ho]; use it if set, else default
    stored_nofs = int(hdr.sdspp[10 + ho])
    hdr.nofs = stored_nofs if stored_nofs > 0 else 1024

    return hdr


def read_ucd_header(path: str | Path) -> JdsHeader:
    """Read a .ucd file header (same 1024-byte format as .jds).

    Parameters
    ----------
    path : str or Path

    Returns
    -------
    JdsHeader
    """
    return read_jds_header(path)


def decode_dspz_frame(
    raw: np.ndarray,
    wofsg: int,
    nofs: int,
    avrs: int,
    mode: int = 1,
) -> np.ndarray:
    """Decode a raw UINT32 data frame to float64.

    Parameters
    ----------
    raw : np.ndarray, dtype uint32, shape (2, wofsg, nofs)
        Raw data frame as read from the .jds file.
    wofsg : int
        Number of frequency channels.
    nofs : int
        Number of time samples (spectra) per frame.
    avrs : int
        Averaging count (from header).
    mode : int
        Channel selection: 0 = ch0-ch1, 1 = ch0 only, 2 = ch1 only.

    Returns
    -------
    imdat : np.ndarray, dtype float64, shape (wofsg, nofs)
        Decoded spectrogram.
    """
    # Decode custom float: bits 5-31 = mantissa, bits 0-4 = exponent
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


class JdsFile:
    """Context-managed reader for .jds data files.

    Usage::

        with JdsFile("data.jds") as jds:
            print(jds.header)
            for frame_idx, imdat in jds.frames(mode=1):
                process(imdat)
    """

    def __init__(self, path: str | Path, nofs: int = 1024):
        self.path = Path(path)
        self.nofs = nofs
        self._fh = None
        self.header: JdsHeader | None = None
        self.nframe: int = 0

    def __enter__(self):
        self.header = read_jds_header(self.path)
        self._fh = open(self.path, "rb")

        file_size = os.path.getsize(self.path)
        wofsg = self.header.wofsg
        frame_bytes = 4 * self.nofs * wofsg * 2  # UINT32 × 2 channels × wofsg × nofs
        self.nframe = (file_size - HEADER_SIZE_BYTES) // frame_bytes
        return self

    def __exit__(self, *exc):
        if self._fh:
            self._fh.close()

    def frames(self, mode: int = 1) -> Iterator[tuple[int, np.ndarray]]:
        """Yield ``(frame_index, decoded_spectrogram)`` for each frame.

        Parameters
        ----------
        mode : int
            Channel selection (0 = difference, 1 = ch0, 2 = ch1).

        Yields
        ------
        frame_idx : int
            Zero-based frame index.
        imdat : np.ndarray, shape (wofsg, nofs)
            Decoded float64 spectrogram.
        """
        wofsg = self.header.wofsg
        avrs = self.header.avrs
        nofs = self.nofs
        frame_bytes = 4 * nofs * wofsg * 2

        self._fh.seek(HEADER_SIZE_BYTES)
        for i in range(self.nframe):
            raw_bytes = self._fh.read(frame_bytes)
            if len(raw_bytes) < frame_bytes:
                break
            raw = np.frombuffer(raw_bytes, dtype=np.uint32).reshape(nofs, wofsg, 2)
            # IDL layout is (2, wofsg, nofs) in column-major = (nofs, wofsg, 2)
            # in C-order.  We need (2, wofsg, nofs):
            raw = raw.transpose(2, 1, 0).copy()
            yield i, decode_dspz_frame(raw, wofsg, nofs, avrs, mode)


def write_ucd_header(
    out_path: str | Path,
    jds_header: JdsHeader,
    nofs: int,
) -> None:
    """Write a .ucd file header (1024 bytes, modified copy of .jds header).

    The only change from the original .jds header is that ``sdspp[10 + HEAD_OFFSET]``
    is set to *nofs* (the number of spectra per frame).

    Parameters
    ----------
    out_path : str or Path
        Path to the output .ucd file (opened in write-binary mode).
    jds_header : JdsHeader
        Parsed header from the first .jds file.
    nofs : int
        Value to store as the nofs parameter.
    """
    # Reconstruct the raw 1024-byte header
    header_bytes = bytearray(HEADER_SIZE_BYTES)
    offset = 0
    for fname, nbytes in HEADER_TEXT_FIELDS.items():
        text = getattr(jds_header, fname).encode("ascii", errors="replace")
        header_bytes[offset: offset + len(text)] = text[:nbytes]
        offset += nbytes

    # Patch nofs into sdspp
    sdspp = jds_header.sdspp.copy()
    sdspp[10 + HEAD_OFFSET] = np.uint32(nofs)

    sdspp_bytes = sdspp.tobytes()
    header_bytes[offset: offset + len(sdspp_bytes)] = sdspp_bytes

    with open(out_path, "wb") as fh:
        fh.write(bytes(header_bytes))
