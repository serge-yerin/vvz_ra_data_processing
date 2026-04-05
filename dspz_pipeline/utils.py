"""Shared utility functions used across the DSPZ pipeline."""

import numpy as np


def decode_text(raw: bytes) -> str:
    """Decode a fixed-width header string, stripping nulls and whitespace.

    Used by both .jds and .ucd header readers.
    """
    return raw.split(b"\x00")[0].decode("ascii", errors="replace").strip()


def smooth_edge(arr: np.ndarray, width: int) -> np.ndarray:
    """Uniform (box-car) smoothing with edge truncation.

    Matches IDL ``smooth(..., /EDGE_TRUNCATE)``.

    Parameters
    ----------
    arr : np.ndarray
        1-D input array.
    width : int
        Smoothing window width.

    Returns
    -------
    np.ndarray
        Smoothed array of the same length.
    """
    if width <= 1:
        return arr.copy()
    n = arr.size
    out = np.empty(n, dtype=arr.dtype)
    half = width // 2
    cs = np.concatenate(([0.0], np.cumsum(arr.astype(np.float64))))
    for i in range(n):
        lo = max(0, i - half)
        hi = min(n - 1, i + half)
        out[i] = (cs[hi + 1] - cs[lo]) / (hi - lo + 1)
    return out
