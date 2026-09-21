"""Shared utility functions used across the DSPZ pipeline."""

import numpy as np


def decode_text(raw: bytes) -> str:
    """Decode a fixed-width header string, stripping nulls and whitespace.

    Used by both .jds and .ucd header readers.
    """
    return raw.split(b"\x00")[0].decode("ascii", errors="replace").strip()


def smooth_edge(arr: np.ndarray, width: int) -> np.ndarray:
    """Uniform (box-car) smoothing with edge truncation.

    Matches IDL ``smooth(..., /EDGE_TRUNCATE)``: every output sample is the
    mean of the input over ``[i - width//2, i + width//2]`` clipped to the
    array bounds.  Vectorised with a cumulative sum; results are bit-identical
    to the original per-element loop.

    Parameters
    ----------
    arr : np.ndarray
        1-D input array, or an N-D array smoothed along its last axis.
    width : int
        Smoothing window width.

    Returns
    -------
    np.ndarray
        Smoothed array of the same shape and dtype.
    """
    if width <= 1:
        return arr.copy()
    n = arr.shape[-1]
    half = width // 2
    cs = np.cumsum(arr.astype(np.float64), axis=-1)
    cs = np.concatenate((np.zeros(arr.shape[:-1] + (1,)), cs), axis=-1)
    i = np.arange(n)
    lo = np.maximum(0, i - half)
    hi = np.minimum(n - 1, i + half)
    out = (cs[..., hi + 1] - cs[..., lo]) / (hi - lo + 1)
    return out.astype(arr.dtype)
