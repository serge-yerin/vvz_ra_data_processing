"""smooth_edge must stay bit-identical to the original per-element loop
(IDL smooth(..., /EDGE_TRUNCATE))."""

import numpy as np
import pytest

from dspz_pipeline.utils import smooth_edge


def _reference(arr, width):
    """The original scalar implementation, kept as the oracle."""
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


@pytest.mark.parametrize("width", [1, 2, 3, 4, 5, 8, 512, 5000])
@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_matches_reference_exactly(width, dtype):
    rng = np.random.default_rng(width)
    arr = rng.standard_normal(3001).astype(dtype)
    out = smooth_edge(arr, width)
    ref = _reference(arr, width)
    assert out.dtype == ref.dtype
    np.testing.assert_array_equal(out, ref)


def test_width_larger_than_array():
    arr = np.arange(10, dtype=np.float64)
    np.testing.assert_array_equal(smooth_edge(arr, 100), _reference(arr, 100))


def test_width_one_returns_copy():
    arr = np.arange(5, dtype=np.float32)
    out = smooth_edge(arr, 1)
    np.testing.assert_array_equal(out, arr)
    assert out is not arr


def test_2d_smooths_along_last_axis():
    rng = np.random.default_rng(0)
    arr = rng.standard_normal((6, 200)).astype(np.float32)
    out = smooth_edge(arr, 4)
    assert out.shape == arr.shape
    for j in range(6):
        np.testing.assert_array_equal(out[j], _reference(arr[j], 4))


def test_is_fast_enough_for_the_gui():
    import time
    arr = np.random.default_rng(0).standard_normal((51, 131072)).astype(np.float32)
    t = time.perf_counter()
    smooth_edge(arr, 4)
    smooth_edge(arr, 512)
    assert time.perf_counter() - t < 1.0
