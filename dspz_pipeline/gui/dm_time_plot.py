"""
Visualisation of the DM-time plane produced by IndSearch.

Can be run as a standalone script to view any existing .dmt file::

    python -m dspz_pipeline.gui.dm_time_plot path/to/file.dmt 12.872

Or called programmatically by indsearch_main after the dedispersion run.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from dspz_pipeline.io.dmt import read_dmt


def plot_dm_time(
    dm_time_plane: np.ndarray,
    dm_const: float,
    dm_step: float = 0.004,
    save_path: str | Path | None = None,
) -> None:
    """Display the DM-time plane as a 2D intensity image plus a 1D time series
    at the DM index with the highest peak signal.

    Parameters
    ----------
    dm_time_plane : ndarray, shape (n_dm, n_time)
        Dedispersed DM-time data as written to the .dmt file.
    dm_const : float
        Central DM [pc/cm^3].
    dm_step : float
        DM step size [pc/cm^3] (default: 0.004).
    save_path : str or Path or None
        Full path to save the PNG. If None the plot is shown but not saved.
    """
    n_dm, n_time = dm_time_plane.shape
    dm_min = dm_const - (n_dm // 2) * dm_step
    dm_max = dm_min + (n_dm - 1) * dm_step

    # 16:9 canvas; GridSpec with a narrow colorbar column only in the top row
    fig = plt.figure(figsize=(16, 9))
    gs = fig.add_gridspec(2, 2, width_ratios=[50, 1], height_ratios=[3, 1])

    ax0 = fig.add_subplot(gs[0, 0])
    cax = fig.add_subplot(gs[0, 1])
    ax1 = fig.add_subplot(gs[1, 0], sharex=ax0)

    # --- DM-time image ---
    im = ax0.imshow(
        dm_time_plane,
        aspect="auto",
        origin="lower",
        extent=[0, n_time - 1, dm_min, dm_max],
        cmap="Greys",
        interpolation="nearest",
    )
    ax0.set_ylabel("DM [pc/cm\u00b3]")
    ax0.set_title("DM\u2013Time plane (IndSearch output)")
    ax0.tick_params(labelbottom=False)   # x labels shown by ax1 below

    # --- Colorbar: tight against ax0, title above, small tick labels ---
    cbar = fig.colorbar(im, cax=cax)
    cbar.ax.tick_params(labelsize=7)
    cax.set_title("Intensity", fontsize=8, pad=3)

    # --- Time series at peak DM ---
    peak_dm_idx = int(dm_time_plane.max(axis=1).argmax())
    peak_dm_val = dm_min + peak_dm_idx * dm_step

    ax1.plot(dm_time_plane[peak_dm_idx, :], linewidth=0.6, color="steelblue")
    ax1.set_xlim(0, n_time - 1)
    ax1.set_xlabel(f"Time sample\nTime series at peak DM = {peak_dm_val:.4f} pc/cm\u00b3")
    ax1.set_ylabel("Intensity")

    fig.subplots_adjust(left=0.06, right=0.97, top=0.93, bottom=0.10, hspace=0, wspace=0)

    if save_path is not None:
        plt.savefig(save_path, dpi=300)
        print(f"Plot saved to: {save_path}")

    plt.show()


def main() -> None:
    """CLI entry point for standalone use."""
    parser = argparse.ArgumentParser(
        description="Plot DM-time plane from a .dmt file",
    )
    parser.add_argument("dmt_file", help="Path to the .dmt file")
    parser.add_argument("dm", type=float, help="Central DM [pc/cm^3]")
    parser.add_argument(
        "--dm-step", type=float, default=0.004,
        help="DM step size [pc/cm^3] (default: 0.004)",
    )
    args = parser.parse_args()

    dmt_path = Path(args.dmt_file)
    dm_time_plane, _, _ = read_dmt(dmt_path)

    save_path = dmt_path.with_suffix(".png")
    plot_dm_time(dm_time_plane, args.dm, dm_step=args.dm_step, save_path=save_path)


if __name__ == "__main__":
    main()
