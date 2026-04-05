"""
Shared I/O for .dmt files and DM delay computation.

The .dmt file format:
  - 4 bytes: int32  DMstepnumb  (number of DM trial values, typically 51)
  - 4 bytes: int32  picsize     (number of time samples)
  - float32 array of shape (DMstepnumb, picsize) in Fortran (column-major) order
"""

from __future__ import annotations

import struct
from pathlib import Path

import numpy as np

from dspz_pipeline.config import (
    DM_HALF_STEPS,
    DM_STEP_SIZE,
    DM_TOTAL_STEPS,
)


def compute_dm_delays(
    dm: float,
    fmax_mhz: float,
    fmin_mhz: float,
    tot_chan: int,
    df_mhz: float,
) -> np.ndarray:
    """Compute per-channel dispersion delays in seconds.

    Translates IDL ``defDMforDSPZ`` from showpulse.pro::

        dt[i] = (DM / 2.4103) * ((1e4 / (fmin + df*(i+1))^2) - (1e4 / fmax^2))

    The delay is relative to the highest frequency channel (``fmax``),
    which has zero delay.  Lower frequency channels have positive delays.

    Parameters
    ----------
    dm : float
        Dispersion measure in pc cm^-3.
    fmax_mhz : float
        Upper edge of the band in MHz.
    fmin_mhz : float
        Lower edge of the band in MHz.
    tot_chan : int
        Total number of frequency channels.
    df_mhz : float
        Channel width in MHz (= bandwidth / tot_chan).

    Returns
    -------
    dt : np.ndarray, shape (tot_chan,), dtype float64
        Delay in seconds for each frequency channel.
    """
    i = np.arange(tot_chan, dtype=np.float64)
    freq = fmin_mhz + df_mhz * (i + 1)
    dt = (dm / 2.4103) * (1.0e4 / freq ** 2 - 1.0e4 / fmax_mhz ** 2)
    return dt


def read_dmt(path: str | Path) -> tuple[np.ndarray, int, int]:
    """Read a .dmt file.

    Parameters
    ----------
    path : str or Path
        Path to the ``.dmt`` file.

    Returns
    -------
    acc_dm : np.ndarray, shape (dm_stepnumb, picsize), dtype float32
        Dedispersed time series for each DM trial.
    dm_stepnumb : int
        Number of DM trial steps.
    picsize : int
        Number of time samples.
    """
    with open(path, "rb") as fh:
        dm_stepnumb = struct.unpack("<i", fh.read(4))[0]
        picsize = struct.unpack("<i", fh.read(4))[0]
        acc_dm = np.fromfile(fh, dtype=np.float32, count=dm_stepnumb * picsize)
        acc_dm = acc_dm.reshape(dm_stepnumb, picsize, order="F")
    return acc_dm, dm_stepnumb, picsize


def write_dmt(
    path: str | Path,
    acc_dm: np.ndarray,
    dm_stepnumb: int,
    picsize: int,
) -> None:
    """Write a .dmt file.

    Parameters
    ----------
    path : str or Path
        Output file path.
    acc_dm : np.ndarray, shape (dm_stepnumb, picsize)
        Dedispersed time series for each DM trial.
    dm_stepnumb : int
        Number of DM trial steps.
    picsize : int
        Number of time samples.
    """
    with open(path, "wb") as fh:
        fh.write(struct.pack("<i", dm_stepnumb))
        fh.write(struct.pack("<i", picsize))
        fh.write(acc_dm.astype(np.float32).tobytes(order="F"))
