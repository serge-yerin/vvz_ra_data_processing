"""
Robust statistical estimators using iterative sigma-clipping.

Translated from IDL files:
  - erov.pro   — iterative 3-sigma clipping for mean and stddev
  - background.pro — iterative 3-sigma clipping for background estimation

Both routines converge by repeatedly discarding outliers beyond 3 sigma
and recomputing the statistics until the mean stabilises.
"""

import numpy as np


def erov(data: np.ndarray) -> tuple[float, float]:
    """Iterative 3-sigma clipped mean and standard deviation.

    Translates IDL ``erov.pro``.  Starting from the full array, the routine
    removes values more than 3 sigma from the current mean, then recomputes.
    Iteration stops when the relative change in the mean is < 1e-5 or when
    the surviving sample count does not change.

    Parameters
    ----------
    data : np.ndarray
        1-D array of values (will not be modified).

    Returns
    -------
    mt : float
        Robust mean (clipped).
    st : float
        Robust standard deviation (clipped).
    """
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

        # Zero out values beyond 3-sigma (IDL: m = m * (abs(m-sr) le srkv*3))
        mask = np.abs(m - sr) <= srkv * 3.0
        m = m * mask

        # Keep only non-zero survivors (IDL: m_in = where(m, count))
        # NOTE: IDL "where(m)" returns indices where m != 0.
        m_in = np.nonzero(m)[0]
        count = m_in.size

        if count == 0:
            return 0.0, 0.0

        m = m[m_in]
        sr = np.mean(m)
        ster = np.std(m, ddof=1) if m.size > 1 else 0.0

        # Convergence check (IDL: abs(sr_pr/sr - 1) lt 1e-5 or count_pr eq count)
        if sr == 0.0:
            break
        if abs(sr_pr / sr - 1.0) < 1e-5 or count_pr == count:
            break

    return sr, ster


def background(
    tab: np.ndarray,
    positive: bool = False,
) -> tuple[float, float, int]:
    """Iterative 3-sigma clipped background and fluctuation estimate.

    Translates IDL ``background.pro``.  Similar to :func:`erov` but returns
    the background level, its 1-sigma fluctuation, and the number of values
    used in the final estimate.

    Parameters
    ----------
    tab : np.ndarray
        1-D array of intensities (will not be modified).
    positive : bool, optional
        If True, only positive values in *tab* are considered.

    Returns
    -------
    fon : float
        Background level (robust mean).
    sigma : float
        Fluctuation level (robust standard deviation).
    ny : int
        Number of values used in the final estimate.
    """
    tab2 = np.asarray(tab, dtype=np.float64).ravel().copy()

    if positive:
        tab2 = tab2[tab2 > 0.0]
        if tab2.size == 0:
            return 0.0, -1.0, 0

    if tab2.size <= 1:
        return float(tab2[0]) if tab2.size == 1 else 0.0, 0.0, tab2.size

    fon = np.mean(tab2)
    sigma = np.std(tab2, ddof=1)  # IDL STDDEV uses sample stddev (N-1)
    ny = 1

    while True:
        test = np.where(np.abs(tab2 - fon) < 3.0 * sigma)[0]
        if test.size <= 1:
            break
        ny = test.size
        moy = np.mean(tab2[test])
        sigma = np.std(tab2[test], ddof=1)
        if moy == fon:
            break
        fon = moy

    return fon, sigma, ny
