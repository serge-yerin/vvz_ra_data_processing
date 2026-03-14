"""
Binary file reader for UTR-2 .ucd observation files.

Header structure (1024 bytes total) matches IDL HStr:
    Sname  : BYTARR(32)    — source name
    Stime  : BYTARR(32)    — local time string
    Sgmtt  : BYTARR(32)    — UTC time string
    Ssysn  : BYTARR(32)    — record identifier
    Ssyst  : BYTARR(32)    — system name
    Splace : BYTARR(96)    — place description
    Sdesc  : BYTARR(256)   — data descriptor
    Sdspp  : ULONARR(128)  — unsigned 32-bit parameter array (512 bytes)
"""

import os
import numpy as np

HEADER_SIZE = 1024
HEAD_OFFSET = 16        # IDL: HeadOffset = 16
FCLK        = 66.0 / 2.0   # MHz — ADC clock / 2


def _decode(b):
    return b.decode("ascii", errors="replace").rstrip("\x00").strip()


def read_header(f):
    """
    Read and parse the 1024-byte file header.

    Parameters
    ----------
    f : file object opened in binary mode, positioned at byte 0

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
    sname  = f.read(32)
    stime  = f.read(32)
    sgmtt  = f.read(32)
    ssysn  = f.read(32)
    ssyst  = f.read(32)   # read but not used (matches IDL Ssyst field)
    splace = f.read(96)
    sdesc  = f.read(256)
    sdspp  = np.frombuffer(f.read(128 * 4), dtype="<u4")   # 512 bytes, little-endian uint32

    # Parameter extraction — indices match IDL: HStr.Sdspp[N + HeadOffset]
    nofs    = int(sdspp[10 + HEAD_OFFSET])                          # spectra per frame
    fmin    = float(sdspp[12 + HEAD_OFFSET]) * FCLK / 8192.0       # MHz
    fmax    = float(sdspp[13 + HEAD_OFFSET]) * FCLK / 8192.0       # MHz
    wofsg   = int(sdspp[14 + HEAD_OFFSET])                         # frequency channels
    avrs    = int(sdspp[15 + HEAD_OFFSET])                         # number of averages
    time_res = avrs * 8192.0 / 2.0 / (FCLK * 1e6)                 # seconds

    return {
        "sname"   : _decode(sname),
        "stime"   : _decode(stime),
        "sgmtt"   : _decode(sgmtt),
        "ssysn"   : _decode(ssysn),
        "splace"  : _decode(splace),
        "sdesc"   : _decode(sdesc),
        "sdspp"   : sdspp,
        "nofs"    : nofs,
        "fmin"    : fmin,
        "fmax"    : fmax,
        "wofsg"   : wofsg,
        "avrs"    : avrs,
        "time_res": time_res,
    }


def count_frames(file_path, wofsg, nofs):
    """
    Return the number of complete data frames in the file.

    Replicates IDL:
        n_kadr = floor((fstatus.size - 1024) / (4L * nofs * wofsg))
    """
    file_size = os.path.getsize(file_path)
    return int((file_size - HEADER_SIZE) // (4 * nofs * wofsg))


def read_frame(f, frame_idx, wofsg, nofs):
    """
    Read one data frame from the open file.

    Replicates IDL:
        datdspz = assoc(1, fltarr(wofsg, nofs), 1024)
        data    = datdspz[j]

    IDL fltarr(wofsg, nofs) stores elements in Fortran (column-major) order:
    element [i_freq, i_time] is at memory offset i_freq + i_time * wofsg.
    Reshaping the flat bytes with order='F' reproduces this layout, so
    the returned array satisfies:
        frame[i_freq, i_time]  ==  IDL data[i_freq, i_time]

    Parameters
    ----------
    f         : binary file object
    frame_idx : int  zero-based frame index
    wofsg     : int  number of frequency channels
    nofs      : int  number of time samples per frame

    Returns
    -------
    ndarray, shape (wofsg, nofs), dtype float32
    """
    offset = HEADER_SIZE + frame_idx * wofsg * nofs * 4
    f.seek(offset)
    raw = np.frombuffer(f.read(wofsg * nofs * 4), dtype="<f4")
    return raw.reshape((wofsg, nofs), order="F")
