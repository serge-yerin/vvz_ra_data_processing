"""
RFI (Radio Frequency Interference) cleaning routines.

Translated from IDL files:
  - adr_cleaning.pro  — top-level cleaning driver
  - launch_sumthr.pro — runs SumThreshold in both directions
  - patrol.pro        — narrow-band and wide-band RFI flagging

The cleaning pipeline:
1. Remove NaN / Inf values
2. Flatten (normalise per frequency channel)
3. SumThreshold flagging (time and frequency directions)
4. Narrow-band flagging (anomalous channels)
5. Wide-band flagging (anomalous time steps)
6. Apply mask to data
"""

import numpy as np

from dspz_pipeline.cleaning.flatten import flatten
from dspz_pipeline.cleaning.robust_stats import background, erov
from dspz_pipeline.cleaning.sumthr import sumthr


# ------------------------------------------------------------------ #
#  launch_sumthr  (from launch_sumthr.pro)
# ------------------------------------------------------------------ #

def launch_sumthr(data: np.ndarray, p: np.ndarray) -> None:
    """Apply SumThreshold RFI flagging in both directions.

    Translates IDL ``launch_sumthr.pro``.

    Parameters
    ----------
    data : np.ndarray, shape (wofsg, nofs)
        Flattened spectrogram (frequency x time).  **Not modified.**
    p : np.ndarray, shape (wofsg, nofs), dtype uint8
        Mask (1 = good, 0 = bad).  **Modified in-place.**
    """
    # Transpose to (nofs, wofsg) — IDL rotate(data,4) on 2-D = transpose
    x = data.T.copy()
    p_t = p.T.copy()

    nt, nf = x.shape
    p1 = p_t.copy()
    p2 = p_t.copy()

    # --- Columns (vertical = along frequency for each time step) ---------- #
    m_vert = np.array([2, 8, 16, 128, 256])
    thr_vert = 10.0 / (1.5 ** (np.log(m_vert) / np.log(2)))

    for i in range(nt):
        fon, sigma, _ny = background(x[i, :])
        y = p2[i, :].copy()
        xx = (x[i, :] - fon) / sigma if sigma != 0 else x[i, :] - fon
        for j in range(len(m_vert)):
            sumthr(xx, int(m_vert[j]), float(thr_vert[j]), y)
        p2[i, :] = y

    # --- Rows (horizontal = along time for each frequency channel) -------- #
    m_horiz = np.array([1, 2, 4, 8, 64])
    thr_horiz = 10.0 / (1.5 ** (np.log(m_horiz) / np.log(2)))

    for i in range(nf):
        fon, sigma, _ny = background(x[:, i])
        y = p1[:, i].copy()
        xx = (x[:, i] - fon) / sigma if sigma != 0 else x[:, i] - fon
        for j in range(len(m_horiz)):
            sumthr(xx, int(m_horiz[j]), float(thr_horiz[j]), y)
        p1[:, i] = y

    # Combine: pixel is good only if good in both passes
    p_t[:] = p1 * p2

    # Transpose back into caller's mask
    p[:] = p_t.T


# ------------------------------------------------------------------ #
#  patrol  (from patrol.pro)
# ------------------------------------------------------------------ #

def patrol(
    imdat: np.ndarray,
    flatdat: np.ndarray,
    med: float,
    p: np.ndarray,
    *,
    wid: bool = False,
    nar: bool = False,
) -> None:
    """Detect narrow-band and wide-band RFI by channel/time-step statistics.

    Translates IDL ``patrol.pro``.

    Parameters
    ----------
    imdat : np.ndarray, shape (wofsg, nofs)
        Original (un-flattened) spectrogram.  **Not modified.**
    flatdat : np.ndarray, shape (wofsg, nofs)
        Flattened spectrogram.
    med : float
        Median value (unused in current IDL — kept for signature compatibility).
    p : np.ndarray, shape (wofsg, nofs), dtype uint8
        Mask (1 = good, 0 = bad).  **Modified in-place.**
    wid : bool
        If True, flag entire time steps whose mean intensity is an outlier.
    nar : bool
        If True, flag entire frequency channels whose stddev/mean ratio
        is an outlier.
    """
    wofsg, nofs = imdat.shape
    l2 = 4.0  # sigma threshold for flagging

    # --- Narrow-band: flag anomalous frequency channels ------------------- #
    if nar:
        s_sk = np.zeros(wofsg, dtype=np.float64)
        m_sk = np.zeros(wofsg, dtype=np.float64)

        for j in range(wofsg):
            good = np.where(p[j, :] == 1)[0]
            if good.size > 0:
                op = imdat[j, good].copy().astype(np.float64)
                mt, st = erov(op)
                s_sk[j] = st
                m_sk[j] = mt

        # Replace zero means with median of all means
        med_msk = np.median(m_sk)
        m_sk[m_sk == 0] = med_msk

        # Ratio of stddev to mean per channel
        op = s_sk / m_sk
        mt, st = erov(op)

        outliers = np.where(np.abs(op - mt) > st * l2)[0]
        if outliers.size > 0:
            p[outliers, :] = 0

    # --- Wide-band: flag anomalous time steps ----------------------------- #
    if wid:
        op = np.zeros(nofs, dtype=np.float64)
        for j in range(nofs):
            w = np.where(p[:, j] == 1)[0]
            if w.size > 0:
                op[j] = np.sum(imdat[w, j]) / np.sum(p[:, j])

        mt, st = erov(op)
        outliers = np.where(np.abs(op - mt) > st * l2)[0]
        if outliers.size > 0:
            p[:, outliers] = 0


# ------------------------------------------------------------------ #
#  adr_cleaning  (from adr_cleaning.pro)
# ------------------------------------------------------------------ #

def adr_cleaning(
    data: np.ndarray,
    *,
    pat: bool = True,
    nar: bool = True,
    wid: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """Clean a 2-D spectrogram by removing RFI.

    Translates IDL ``adr_cleaning.pro``.

    Parameters
    ----------
    data : np.ndarray, shape (wofsg, nofs)
        Spectrogram (frequency x time).  **Modified in-place** — on return
        it contains the cleaned, flattened, masked data.
    pat : bool
        Run SumThreshold flagging.
    nar : bool
        Run narrow-band flagging.
    wid : bool
        Run wide-band flagging.

    Returns
    -------
    data : np.ndarray
        The same array passed in, now cleaned.
    mask : np.ndarray, dtype uint8
        Binary mask (1 = good, 0 = flagged).
    """
    if data.ndim != 2:
        raise ValueError("Data should be 2D!")

    wofsg, nofs = data.shape
    err = False

    # --- NaN / Inf removal ------------------------------------------------ #
    mask = np.ones((wofsg, nofs), dtype=np.uint8)
    bad = ~np.isfinite(data)
    if np.any(bad):
        mask[bad] = 0
        data[bad] = np.median(data[np.isfinite(data)])

    if data.min() == data.max():
        err = True

    # Keep a copy of the original (un-flattened) data for patrol
    imdat_orig = data.copy()

    # --- Flatten (per-channel normalisation) ------------------------------ #
    flatarr = data.copy().astype(np.float64)
    flatten(flatarr)

    med = 0.0

    # --- SumThreshold flagging -------------------------------------------- #
    if pat:
        launch_sumthr(flatarr, mask)

    # --- Narrow / wide-band flagging -------------------------------------- #
    if wid or nar:
        patrol(imdat_orig, flatarr, med, mask, wid=wid, nar=nar)

    # --- Apply mask ------------------------------------------------------- #
    data[:] = flatarr * mask

    if err:
        data[:] = 0.0

    return data, mask
