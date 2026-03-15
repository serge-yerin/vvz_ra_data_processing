"""
Visualisation of the DM-time plane produced by IndSearch.
"""

import numpy as np
import matplotlib.pyplot as plt


def plot_dm_time(dm_time_plane, dm_const, dm_step_numb, dm_step=0.004):
    """
    Display the DM-time plane as a 2D intensity image plus a 1D time series
    at the DM index with the highest peak signal.

    Parameters
    ----------
    dm_time_plane : ndarray, shape (dm_step_numb+1, n_time)
        Dedispersed, DM-summed data as written to the .dmt file.
    dm_const      : float  Central DM [pc/cm^3]
    dm_step_numb  : int    Number of DM steps (DM_STEP_NUMB constant)
    dm_step       : float  DM step size [pc/cm^3]
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
    plt.savefig("dm_time_plane.png", dpi=300)
    plt.show()


if __name__ == "__main__":
    
    # output_path = "data/Cleaned_ PSRB0834p06A141010_032001.jds.ucd.dmt"
    # output_path = "data//Cleaned_ PSRB0834p06A141010_032001.jds DM=12.872.ucd.dmt"
    output_path = "../../../Survey_Processing/Cleaned_ PSRB0834p06A141010_032001.jds.ucd.dmt"
    
    
    with open(output_path, "rb") as fout:
        fout.seek(2 * 4)  # Skip the first 2 int32 values (header)
        dm_time_plane = np.fromfile(fout, dtype="<f4").T.reshape((51, 65536), order='F')
    
    plot_dm_time(dm_time_plane, 12.872, dm_step_numb=50, dm_step=0.004)