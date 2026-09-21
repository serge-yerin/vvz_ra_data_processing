"""Smoke tests for the Transient Search window (needs a display)."""

import time

import numpy as np
import pytest

from dspz_pipeline.io.dmt import write_dmt

tk = pytest.importorskip("tkinter")


@pytest.fixture
def dmt_file(tmp_path):
    n_dm, n_t = 51, 70000
    rng = np.random.default_rng(0)
    arr = rng.standard_normal((n_dm, n_t)).astype(np.float32)
    arr[25, 30000] += 40.0
    path = tmp_path / "tiny.dmt"
    write_dmt(path, arr, n_dm, n_t)
    return str(path)


@pytest.fixture
def app(dmt_file):
    from dspz_pipeline.gui.trans_search import TransSearchApp
    try:
        a = TransSearchApp(dmt_file, 12.872)
    except tk.TclError as exc:
        pytest.skip(f"Tk unavailable: {exc}")
    yield a
    a.root.destroy()


def _pump(root, seconds=0.15):
    end = time.time() + seconds
    while time.time() < end:
        root.update()
        time.sleep(0.01)


def test_sixteen_panels_built_once(app):
    assert len(app.panel_images) == 16
    imgs = list(app.panel_images)
    app._redraw()
    assert list(app.panel_images) == imgs        # same AxesImage objects reused


def test_slider_updates_clim_in_place(app):
    imgs = list(app.panel_images)
    app.scl_min.set(-3.0)
    app.scl_max.set(7.0)
    app._on_scale_change()
    _pump(app.root)
    assert list(app.panel_images) == imgs
    assert all(im.get_clim() == (-3.0, 7.0) for im in app.panel_images)


def test_dm_step_changes_title_only(app):
    imgs = list(app.panel_images)
    data_before = app.panel_images[0].get_array().copy()
    app._adj_dm(+1)
    assert "DM step=+1" in app.fig._suptitle.get_text()
    assert list(app.panel_images) == imgs
    np.testing.assert_array_equal(app.panel_images[0].get_array(), data_before)


def test_smpar_change_updates_image_data(app):
    imgs = list(app.panel_images)
    data_before = app.panel_images[7].get_array().copy()
    app._adj_smpar(+2)       # 4 -> 6 (widths 4 and 5 share the same half-window)
    assert list(app.panel_images) == imgs
    assert not np.array_equal(app.panel_images[7].get_array(), data_before)


def test_parts_overlay_follows_selection(app):
    def visible():
        return [i for i, span in enumerate(app.panel_spans)
                if span is not None and span.get_width() > 0]
    assert visible() == list(range(16))     # 1 part: everything highlighted
    app._adj_nofp(+1)                       # 2 parts, part 0 -> samples 0..32768
    _pump(app.root)
    assert visible() == list(range(8))
