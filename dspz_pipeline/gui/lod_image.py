"""
Level-of-detail rendering helpers for very large 2-D images in matplotlib.

Matplotlib's ``imshow`` normalises and resamples the *whole* array on every
redraw, regardless of the current zoom.  For a 4096 x 44000 float32
spectrogram that is ~10 s per draw.  The helpers here keep the full-resolution
array in memory but hand matplotlib only the currently visible sub-window,
block-pooled down to roughly the size of the axes in pixels, so every redraw
(colour-limit change, zoom, pan) takes tens of milliseconds.
"""

from __future__ import annotations

import math

import numpy as np

POOL_METHODS = ("max", "mean", "nearest")


def pool_window(
    data: np.ndarray,
    col_lo: float,
    col_hi: float,
    row_lo: float,
    row_hi: float,
    max_cols: int,
    max_rows: int,
    method: str = "max",
) -> tuple[np.ndarray, tuple[int, int, int, int]]:
    """Extract ``data[row_lo:row_hi, col_lo:col_hi]`` pooled to at most
    ``(max_rows, max_cols)`` pixels.

    Fractional bounds are expanded outward to whole samples and clamped to the
    array.  A degenerate (empty) window is widened to one sample.  Each output
    pixel covers an integer block of ``fy x fx`` input samples; a trailing
    partial block is dropped, so the returned extent ``(col_lo, col_hi,
    row_lo, row_hi)`` (in index units) describes exactly the samples the
    pooled array covers.

    ``method``: ``"max"`` (block maximum), ``"mean"`` (block average) or
    ``"nearest"`` (first sample of each block, like plain sub-sampling).
    """
    if method not in POOL_METHODS:
        raise ValueError(f"method must be one of {POOL_METHODS}, got {method!r}")

    n_rows, n_cols = data.shape
    c0, c1 = _clamp_span(col_lo, col_hi, n_cols)
    r0, r1 = _clamp_span(row_lo, row_hi, n_rows)

    fx = max(1, math.ceil((c1 - c0) / max(1, max_cols)))
    fy = max(1, math.ceil((r1 - r0) / max(1, max_rows)))
    nx = max(1, (c1 - c0) // fx)
    ny = max(1, (r1 - r0) // fy)
    c1 = c0 + nx * fx
    r1 = r0 + ny * fy

    sub = data[r0:r1, c0:c1]
    if fx == 1 and fy == 1:
        pooled = sub
    elif method == "nearest":
        pooled = sub[::fy, ::fx]
    else:
        blocks = sub.reshape(ny, fy, nx, fx)
        pooled = blocks.max(axis=(1, 3)) if method == "max" else blocks.mean(axis=(1, 3))

    return pooled, (c0, c1, r0, r1)


def _clamp_span(lo: float, hi: float, n: int) -> tuple[int, int]:
    """Expand ``[lo, hi)`` outward to whole indices inside ``[0, n)``,
    guaranteeing at least one index."""
    lo_i = max(0, min(int(math.floor(lo)), n - 1))
    hi_i = max(lo_i + 1, min(int(math.ceil(hi)), n))
    return lo_i, hi_i


class LodImage:
    """A matplotlib ``AxesImage`` that always shows a pooled view of the
    currently visible part of a large 2-D array.

    The image is created once; afterwards only ``set_data``/``set_extent``
    are called on it, so the axes limits (zoom, pan, toolbar history) are
    never disturbed.  Any change of the axes limits triggers a re-pool of
    the visible window.

    ``x_extent``/``y_extent`` give the data-coordinate span of the whole
    array along columns and rows (``origin="lower"`` convention, i.e. row 0
    sits at ``y_extent[0]``).
    """

    def __init__(
        self,
        ax,
        data: np.ndarray,
        x_extent: tuple[float, float],
        y_extent: tuple[float, float],
        method: str = "max",
        oversample: float = 1.0,
        **imshow_kwargs,
    ):
        if method not in POOL_METHODS:
            raise ValueError(f"method must be one of {POOL_METHODS}, got {method!r}")
        self.ax = ax
        self.data = data
        self.x_extent = (float(x_extent[0]), float(x_extent[1]))
        self.y_extent = (float(y_extent[0]), float(y_extent[1]))
        self.method = method
        self.oversample = oversample
        self._busy = False
        self._last_key = None

        n_rows, n_cols = data.shape
        self._dx = (self.x_extent[1] - self.x_extent[0]) / n_cols
        self._dy = (self.y_extent[1] - self.y_extent[0]) / n_rows

        imshow_kwargs.setdefault("aspect", "auto")
        imshow_kwargs.setdefault("interpolation", "nearest")
        imshow_kwargs["origin"] = "lower"
        ax.set_xlim(*self.x_extent)
        ax.set_ylim(*self.y_extent)
        ax.set_autoscale_on(False)

        pooled, idx_extent = self._pool_for(self.x_extent, self.y_extent, self._axes_px())
        self.image = ax.imshow(pooled, extent=self._to_data_extent(idx_extent), **imshow_kwargs)
        self._last_key = self._key(self.x_extent, self.y_extent, self._axes_px())

        ax.callbacks.connect("xlim_changed", self._on_lim_changed)
        ax.callbacks.connect("ylim_changed", self._on_lim_changed)

    # -- public ---------------------------------------------------------- #

    def set_clim(self, vmin: float, vmax: float) -> None:
        self.image.set_clim(vmin, vmax)

    def set_method(self, method: str) -> None:
        if method not in POOL_METHODS:
            raise ValueError(f"method must be one of {POOL_METHODS}, got {method!r}")
        self.method = method
        self.refresh(force=True)

    def refresh(self, target_px: tuple[float, float] | None = None, force: bool = False) -> None:
        """Re-pool the visible window; ``target_px`` overrides the axes size.

        Does nothing if the view, target size and method are unchanged since
        the last pool (unless ``force``), so repeated resize/limit events are
        cheap.
        """
        if self._busy:
            return
        self._busy = True
        try:
            xlim = sorted(self.ax.get_xlim())
            ylim = sorted(self.ax.get_ylim())
            target_px = target_px or self._axes_px()
            key = self._key(xlim, ylim, target_px)
            if key == self._last_key and not force:
                return
            self._last_key = key
            pooled, idx_extent = self._pool_for(xlim, ylim, target_px)
            self.image.set_data(pooled)
            self.image.set_extent(self._to_data_extent(idx_extent))
        finally:
            self._busy = False

    def render_at(self, width_px: float, height_px: float):
        """Context manager: show the visible window pooled for an output of
        ``width_px x height_px`` (e.g. for ``savefig``), then restore."""
        return _RenderAt(self, (width_px, height_px))

    # -- internals ------------------------------------------------------- #

    def _key(self, xlim, ylim, target_px):
        return (tuple(xlim), tuple(ylim), tuple(target_px), self.method)

    def _axes_px(self) -> tuple[float, float]:
        bbox = self.ax.get_window_extent()
        return max(1.0, bbox.width * self.oversample), max(1.0, bbox.height * self.oversample)

    def _pool_for(self, xlim, ylim, target_px):
        col_lo = (xlim[0] - self.x_extent[0]) / self._dx
        col_hi = (xlim[1] - self.x_extent[0]) / self._dx
        row_lo = (ylim[0] - self.y_extent[0]) / self._dy
        row_hi = (ylim[1] - self.y_extent[0]) / self._dy
        return pool_window(
            self.data, col_lo, col_hi, row_lo, row_hi,
            max_cols=int(math.ceil(target_px[0])), max_rows=int(math.ceil(target_px[1])),
            method=self.method,
        )

    def _to_data_extent(self, idx_extent):
        c0, c1, r0, r1 = idx_extent
        return (
            self.x_extent[0] + c0 * self._dx, self.x_extent[0] + c1 * self._dx,
            self.y_extent[0] + r0 * self._dy, self.y_extent[0] + r1 * self._dy,
        )

    def _on_lim_changed(self, _ax):
        self.refresh()


class _RenderAt:
    def __init__(self, lod: LodImage, target_px):
        self.lod = lod
        self.target_px = target_px

    def __enter__(self):
        self.lod.refresh(self.target_px)
        return self.lod

    def __exit__(self, *exc):
        self.lod.refresh()
        return False
