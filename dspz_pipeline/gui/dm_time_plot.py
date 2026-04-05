"""
Visualisation of the DM-time plane produced by IndSearch (Stage 2).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


def plot_dm_time(
    dm_time_plane: np.ndarray,
    dm_const: float,
    dm_step_numb: int,
    dm_step: float = 0.004,
    output_dir: str | Path = "_output",
) -> None:
    """
    Display the DM-time plane as a 2D intensity image plus a 1D time series
    at the DM index with the highest peak signal.

    Parameters
    ----------
    dm_time_plane : ndarray, shape (dm_step_numb+1, n_time)
        Dedispersed, DM-summed data as written to the .dmt file.
    dm_const : float
        Central DM [pc/cm^3].
    dm_step_numb : int
        Number of DM steps.
    dm_step : float
        DM step size [pc/cm^3].
    output_dir : str or Path
        Directory where the plot image is saved.
    """
    n_dm, n_time = dm_time_plane.shape
    dm_min = dm_const - (dm_step_numb // 2) * dm_step
    dm_max = dm_min + dm_step_numb * dm_step

    fig, axes = plt.subplots(
        2, 1,
        figsize=(11, 5.5),
        gridspec_kw={"height_ratios": [3, 1]},
    )

    # --- DM-time image ---
    ax0 = axes[0]
    im = ax0.imshow(
        dm_time_plane,
        aspect="auto",
        origin="lower",
        extent=[0, n_time - 1, dm_min, dm_max],
        cmap="Grays",
        interpolation="nearest",
    )
    ax0.set_xlabel("Time sample")
    ax0.set_ylabel("DM [pc/cm\u00b3]")
    ax0.set_title("DM\u2013Time plane (IndSearch output)")
    plt.colorbar(im, ax=ax0, label="Intensity")

    # --- Time series at peak DM ---
    peak_dm_idx = int(dm_time_plane.max(axis=1).argmax())
    peak_dm_val = dm_min + peak_dm_idx * dm_step

    ax1 = axes[1]
    ax1.plot(dm_time_plane[peak_dm_idx, :], linewidth=0.6, color="steelblue")
    ax1.set_xlabel("Time sample")
    ax1.set_ylabel("Intensity")
    ax1.set_title(f"Time series at peak DM = {peak_dm_val:.4f} pc/cm\u00b3")

    plt.tight_layout()

    # Save to output directory
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    save_path = out_path / "dm_time_plane.png"
    plt.savefig(save_path, dpi=300)
    print(f"Plot saved to: {save_path}")

    plt.show()
