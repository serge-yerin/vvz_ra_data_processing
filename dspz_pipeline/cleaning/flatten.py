"""
Per-frequency-channel normalisation (flattening) of spectrograms.

Translated from IDL ``flatten.pro``.

For each frequency channel the routine computes the robust mean via
:func:`~dspz_pipeline.cleaning.robust_stats.erov` and normalises the channel
as ``(data - mean) / mean``.  This removes the bandpass shape so that
all channels have comparable fluctuation levels for subsequent RFI
flagging.
"""

import numpy as np

from dspz_pipeline.cleaning.robust_stats import erov


def flatten(imdat: np.ndarray) -> None:
    """Flatten a 2-D spectrogram in-place.

    Parameters
    ----------
    imdat : np.ndarray, shape (wofsg, nofs)
        Spectrogram array (frequency x time).  Modified **in-place**:
        each frequency row is replaced by ``(row - mean) / mean``
        where *mean* is the robust (sigma-clipped) mean of that row.
    """
    wofsg = imdat.shape[0]

    for j in range(wofsg):
        row = imdat[j, :].copy().astype(np.float64)
        mt, _st = erov(row)
        if mt != 0.0:
            imdat[j, :] = (row - mt) / mt
        else:
            # Avoid division by zero — leave row as-is
            imdat[j, :] = row - mt
