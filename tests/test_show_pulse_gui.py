"""Smoke tests for the Cleaned-data panel of ShowPulseApp (needs a display)."""

import time

import numpy as np
import pytest

from dspz_pipeline.config import HEADER_SIZE_BYTES

tk = pytest.importorskip("tkinter")


@pytest.fixture
def ucd_file(tmp_path):
    n_time, wofsg = 600, 4096
    rng = np.random.default_rng(0)
    arr = rng.standard_normal((n_time, wofsg)).astype(np.float32)
    arr[300, :] += 30.0          # a bright "pulse"
    path = tmp_path / "tiny.ucd"
    with open(path, "wb") as fh:
        fh.write(b"\0" * HEADER_SIZE_BYTES)
        fh.write(arr.tobytes())
    return str(path), n_time


@pytest.fixture
def app(ucd_file):
    from dspz_pipeline.gui.show_pulse import ShowPulseApp
    path, picsize = ucd_file
    try:
        a = ShowPulseApp(path + ".dmt", 12.872, 25, ns=200, picsize=picsize, smpar=4)
    except tk.TclError as exc:          # no display available
        pytest.skip(f"Tk unavailable: {exc}")
    yield a
    a._on_close()


def _pump(root, seconds=0.15):
    end = time.time() + seconds
    while time.time() < end:
        root.update()
        time.sleep(0.01)


def test_default_pooling_is_max(app):
    assert app.data_lod.method == "max"


def test_slider_change_preserves_zoom(app):
    ax = app.data_ax
    ax.set_xlim(100, 150)
    ax.set_ylim(20.0, 21.0)
    app.vmin_var.set(0.5)
    app._on_vmin_change("0.5")
    _pump(app.root)
    assert ax.get_xlim() == (100, 150)
    assert ax.get_ylim() == (20.0, 21.0)
    # clim is applied in raw data units: mean + v * std
    lo, hi = app.data_lod.image.get_clim()
    assert lo == pytest.approx(app._data_mean + 0.5 * app._data_std)


def test_dm_change_does_not_touch_data_panel(app):
    img_before = app.data_lod.image
    ax = app.data_ax
    ax.set_xlim(100, 150)
    app._adj_dm(+1)
    assert app.data_lod.image is img_before
    assert ax.get_xlim() == (100, 150)


def test_pooling_switch_updates_lod(app):
    app.pool_var.set("mean")
    app._on_pool_change()
    assert app.data_lod.method == "mean"


def test_axes_fill_window_after_resize(app):
    win = app.data_win
    win.geometry("1400x900")
    _pump(app.root, 0.4)          # let the debounced resize + redraw run
    app.data_canvas.draw()
    pos = app.data_ax.get_position()
    assert pos.height > 0.85, f"axes only {pos.height:.2f} of figure height"
    assert pos.width > 0.85, f"axes only {pos.width:.2f} of figure width"
