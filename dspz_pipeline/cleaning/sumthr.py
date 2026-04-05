"""
SumThreshold RFI flagging algorithm.

Translated from IDL ``SumThr.pro``.

The SumThreshold method flags runs of *M* consecutive good pixels whose
sum exceeds a specified threshold.  It is applied repeatedly with
increasing window sizes to catch RFI of various widths.  See
Offringa et al. (2010, MNRAS 405, 155) for the underlying idea.
"""

import numpy as np


def sumthr(
    x: np.ndarray,
    m: int,
    thr: float,
    y: np.ndarray,
) -> None:
    """Flag sequences whose running sum exceeds a threshold.

    Parameters
    ----------
    x : np.ndarray, shape (N,)
        1-D array of normalised values (zero-mean, unit-sigma expected).
    m : int
        Window size — number of consecutive samples to sum.
    thr : float
        Threshold **per sample** — the running sum is compared against
        ``thr * m``.
    y : np.ndarray, shape (N,), dtype uint8 or similar
        Mask of good pixels (1 = good, 0 = bad).  **Modified in-place.**
        Only pixels currently marked good participate in the test; newly
        flagged pixels are set to 0.
    """
    w = np.where(y == 1)[0]
    if w.size == 0:
        return

    xx = x[w]
    nxx = xx.size
    if nxx < m:
        return

    # Running sum of M consecutive good values
    xxa = np.zeros(nxx - m + 1, dtype=np.float64)
    for i in range(m):
        xxa += xx[i: i + nxx - m + 1]

    p = np.ones(nxx, dtype=np.uint8)
    wa = np.where(xxa > thr * m)[0]
    if wa.size > 0:
        for i in range(m):
            p[wa + i] = 0

    # Write back into the original mask at the positions of good pixels
    y[w] = p
