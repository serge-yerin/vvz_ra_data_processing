"""Tests for the level-of-detail image helpers used by the Cleaned-data viewer."""

import numpy as np
import pytest

from dspz_pipeline.gui.lod_image import pool_window


@pytest.fixture
def data():
    # 8 rows (freq channels) x 20 cols (time samples), unique values
    return np.arange(8 * 20, dtype=np.float32).reshape(8, 20)


def test_small_window_returned_at_full_resolution(data):
    pooled, extent = pool_window(data, 3, 9, 2, 6, max_cols=100, max_rows=100)
    np.testing.assert_array_equal(pooled, data[2:6, 3:9])
    assert extent == (3, 9, 2, 6)


def test_max_pooling_keeps_single_spike(data):
    arr = np.zeros((4, 20), dtype=np.float32)
    arr[1, 7] = 50.0
    pooled, extent = pool_window(arr, 0, 20, 0, 4, max_cols=5, max_rows=2, method="max")
    assert pooled.shape == (2, 5)
    assert pooled[0, 1] == 50.0          # spike lands in block row 0 (rows 0-1), col 1 (cols 4-7)
    assert pooled.sum() == 50.0
    assert extent == (0, 20, 0, 4)


def test_mean_pooling_averages_blocks(data):
    pooled, _ = pool_window(data, 0, 20, 0, 8, max_cols=10, max_rows=4, method="mean")
    expected = data.reshape(4, 2, 10, 2).mean(axis=(1, 3))
    np.testing.assert_allclose(pooled, expected)


def test_nearest_pooling_takes_first_sample_of_block(data):
    pooled, _ = pool_window(data, 0, 20, 0, 8, max_cols=10, max_rows=4, method="nearest")
    np.testing.assert_array_equal(pooled, data[::2, ::2])


def test_bounds_outside_array_are_clamped(data):
    pooled, extent = pool_window(data, -5, 100, -2, 50, max_cols=100, max_rows=100)
    np.testing.assert_array_equal(pooled, data)
    assert extent == (0, 20, 0, 8)


def test_trailing_partial_block_is_trimmed_from_extent(data):
    # width 20 into max 3 cols -> factor 7 -> 2 full blocks (14 cols), 6 trailing cols dropped
    pooled, extent = pool_window(data, 0, 20, 0, 8, max_cols=3, max_rows=8)
    assert pooled.shape == (8, 2)
    assert extent == (0, 14, 0, 8)


def test_degenerate_window_yields_one_pixel(data):
    pooled, extent = pool_window(data, 5, 5, 3, 3, max_cols=100, max_rows=100)
    assert pooled.shape == (1, 1)
    assert pooled[0, 0] == data[3, 5]
    assert extent == (5, 6, 3, 4)


def test_window_at_far_edge_yields_one_pixel(data):
    pooled, extent = pool_window(data, 20, 25, 8, 9, max_cols=100, max_rows=100)
    assert pooled.shape == (1, 1)
    assert pooled[0, 0] == data[7, 19]
    assert extent == (19, 20, 7, 8)


def test_fractional_bounds_are_expanded_to_whole_samples(data):
    pooled, extent = pool_window(data, 2.3, 6.7, 1.9, 4.1, max_cols=100, max_rows=100)
    np.testing.assert_array_equal(pooled, data[1:5, 2:7])
    assert extent == (2, 7, 1, 5)


def test_unknown_method_rejected(data):
    with pytest.raises(ValueError):
        pool_window(data, 0, 20, 0, 8, max_cols=5, max_rows=4, method="median")


# ---------------------------------------------------------------------- #
#  LodImage (matplotlib AxesImage manager)
# ---------------------------------------------------------------------- #

import matplotlib  # noqa: E402

matplotlib.use("Agg")
from matplotlib.backends.backend_agg import FigureCanvasAgg  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402

from dspz_pipeline.gui.lod_image import LodImage  # noqa: E402


@pytest.fixture
def big():
    # 64 channels x 4000 samples, y axis is 16.5..33.0 MHz like the real data
    rng = np.random.default_rng(1)
    return rng.standard_normal((64, 4000)).astype(np.float32)


@pytest.fixture
def ax():
    fig = Figure(figsize=(4, 2), dpi=100)   # axes will be a few hundred px wide
    FigureCanvasAgg(fig)
    return fig.add_subplot(1, 1, 1)


def test_initial_image_is_pooled_to_axes_size(big, ax):
    lod = LodImage(ax, big, x_extent=(0, 4000), y_extent=(16.5, 33.0), oversample=2)
    w_px = ax.get_window_extent().width
    shown = lod.image.get_array()
    assert shown.shape[1] <= 2 * w_px + 1
    assert shown.shape[1] < 4000
    x0, x1, y0, y1 = lod.image.get_extent()
    assert x0 == 0 and 4000 - x1 < 4000 / shown.shape[1]   # only a sub-pixel trailing block dropped
    assert (y0, y1) == (16.5, 33.0)
    assert ax.get_xlim() == (0, 4000)
    assert ax.get_ylim() == (16.5, 33.0)


def test_zoom_in_shows_full_resolution_and_keeps_limits(big, ax):
    lod = LodImage(ax, big, x_extent=(0, 4000), y_extent=(16.5, 33.0))
    ax.set_xlim(1000, 1050)
    ax.set_ylim(16.5, 16.5 + 8 * (16.5 / 64))   # channels 0..8
    shown = lod.image.get_array()
    np.testing.assert_array_equal(shown, big[0:8, 1000:1050])
    x0, x1, y0, y1 = lod.image.get_extent()
    assert (x0, x1) == (1000, 1050)
    assert y0 == pytest.approx(16.5)
    assert y1 == pytest.approx(16.5 + 8 * (16.5 / 64))
    assert ax.get_xlim() == (1000, 1050)


def test_method_switch_changes_pooling(big, ax):
    lod = LodImage(ax, big, x_extent=(0, 4000), y_extent=(16.5, 33.0), method="max")
    shown_max = lod.image.get_array().copy()
    lod.set_method("mean")
    shown_mean = lod.image.get_array()
    assert shown_max.shape == shown_mean.shape
    assert shown_max.mean() > shown_mean.mean()   # block max of N(0,1) > block mean


def test_set_clim_forwards_to_image(big, ax):
    lod = LodImage(ax, big, x_extent=(0, 4000), y_extent=(16.5, 33.0))
    lod.set_clim(-2.0, 7.0)
    assert lod.image.get_clim() == (-2.0, 7.0)


def test_render_at_temporarily_increases_resolution(big, ax):
    lod = LodImage(ax, big, x_extent=(0, 4000), y_extent=(16.5, 33.0))
    before = lod.image.get_array().shape
    with lod.render_at(4000, 64):
        inside = lod.image.get_array()
        np.testing.assert_array_equal(inside, big)
    after = lod.image.get_array().shape
    assert after == before


def test_refresh_without_change_does_not_repool(big, ax):
    lod = LodImage(ax, big, x_extent=(0, 4000), y_extent=(16.5, 33.0))
    shown = lod.image.get_array()
    lod.refresh()
    assert lod.image.get_array() is shown          # set_data not called again
    ax.set_xlim(1000, 1050)
    assert lod.image.get_array() is not shown      # a real change still re-pools
