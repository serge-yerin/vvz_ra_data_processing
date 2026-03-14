"""
Dispersion measure delay computation.
Replicates IDL procedure defDMforDSPZ.
"""

import numpy as np


def compute_dm_delays(dm, fmax, fmin, tot_chan, df):
    """
    Compute dispersion-measure time delays for each frequency channel.

    IDL formula (loop i = 0 .. tot_chan-1):
        dt[i] = (DM / 2.4103) * (1e4 / (fmin + df*(i+1))^2 - 1e4 / fmax^2)

    Note: the channel frequency uses i+1 (1-based), so frequencies run
    from fmin+df to fmin+df*tot_chan = fmax.

    Parameters
    ----------
    dm       : float  Dispersion measure [pc/cm^3]
    fmax     : float  Maximum frequency [MHz]
    fmin     : float  Minimum frequency [MHz]
    tot_chan : int    Number of frequency channels
    df       : float  Frequency step per channel [MHz]

    Returns
    -------
    dt : ndarray, shape (tot_chan,), dtype float32
        Time delay for each channel relative to fmax [seconds].
        Positive values mean the channel arrives later (lower frequency).
    """
    i = np.arange(tot_chan, dtype=np.float64)
    freq = fmin + df * (i + 1)          # IDL: fmin + df*(i+1), channel 0 → fmin+df
    dt = (dm / 2.4103) * (1e4 / freq ** 2 - 1e4 / fmax ** 2)
    return dt.astype(np.float32)
