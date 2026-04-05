"""
Repeating pulse analysis via FFT.

Translated from IDL ``repeatinganalisys.pro``.

Computes the FFT power spectrum of the selected time window from the
dedispersed data to search for periodic (repeating) pulsar signals.
Displays S/N vs time and the FFT power spectrum.

Usage::

    Called from TransSearch GUI when the REP button is active, or::

        python -m dspz_pipeline.gui.repeating_analysis <dmt_file> <dm>
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("TkAgg")

import matplotlib.pyplot as plt
import numpy as np

from dspz_pipeline.config import DEFAULT_DM_CONST
from dspz_pipeline.utils import smooth_edge


def repeating_analysis(
    acc_dm_norm: np.ndarray,
    dm_pos: int,
    minscl: float,
    maxscl: float,
    p_viz: np.ndarray,
) -> None:
    """Analyse repeating pulses via FFT of the selected time window.

    Translates IDL ``RepeatingAnalisys`` from repeatinganalisys.pro.

    Parameters
    ----------
    acc_dm_norm : np.ndarray, shape (dm_stepnumb, picsize)
        Normalised dedispersed data.
    dm_pos : int
        DM step index to analyse.
    minscl, maxscl : float
        Clipping range for the time series.
    p_viz : np.ndarray, shape (picsize,) or (N_PANELS * KADR,)
        Visualisation mask — non-zero values mark the selected time window.
    """
    picsize = acc_dm_norm.shape[1]

    # Extract the time series at the selected DM, clipped to [minscl, maxscl]
    timeline = np.clip(acc_dm_norm[dm_pos, :], minscl, maxscl).copy()

    # Build the windowing function from pViz
    n_total = 65536
    win_long = np.zeros(n_total, dtype=np.float64)

    n_panels = 16
    kadr = 4096
    p_viz_trimmed = p_viz[:n_panels * kadr] if len(p_viz) >= n_panels * kadr else np.pad(p_viz, (0, n_panels * kadr - len(p_viz)))

    win_fun = np.mean(p_viz_trimmed.reshape(n_panels, -1), axis=1)
    for i in range(n_panels):
        win_long[i * kadr:(i + 1) * kadr] = win_fun[i] / 200.0

    # Find the active window
    active = np.where(win_long > 0)[0]
    if active.size == 0:
        print("RepeatingAnalysis: no active time window selected.")
        return
    plot_min = int(active[0])
    plot_max = int(active[-1])

    # --- Window 1: S/N vs Time ---------------------------------------- #
    fig, (ax_sn, ax_fft) = plt.subplots(1, 2, figsize=(10, 4))

    active_timeline = timeline[active]
    ax_sn.plot(active, active_timeline, ".", markersize=2)
    ax_sn.set_xlabel("Time sample")
    ax_sn.set_ylabel("S/N")
    ax_sn.set_title("S/N vs Time")

    # --- Window 2: FFT ------------------------------------------------ #
    windowed = timeline[:n_total] * win_long[:min(n_total, picsize)]
    fft_line = np.abs(np.fft.fft(np.clip(windowed, minscl, maxscl))) ** 2

    # Zero out the lowest bins (DC and very low frequencies)
    fft_line[:41] = fft_line[8000:8041] if len(fft_line) > 8040 else 0

    ax_fft.plot(fft_line[:4001])
    ax_fft.set_xlim(0, 4000)
    ax_fft.set_xlabel("FFT bin")
    ax_fft.set_ylabel("Power")
    ax_fft.set_title("FFT power spectrum")

    fig.tight_layout()
    fig.show()


def main() -> None:
    """CLI entry point for ``python -m dspz_pipeline.gui.repeating_analysis``."""
    parser = argparse.ArgumentParser(description="Repeating pulse FFT analysis")
    parser.add_argument("dmt_file", help="Path to the .dmt data file")
    parser.add_argument("dm", type=float, help="Central DM in pc/cm^3")
    parser.add_argument("--dm-pos", type=int, default=25, help="DM step index (default: 25)")
    args = parser.parse_args()

    from dspz_pipeline.io.dmt import read_dmt
    acc_dm_file, dm_stepnumb, picsize = read_dmt(args.dmt_file)

    # Simple smoothing for standalone use
    acc_dm = np.zeros_like(acc_dm_file)
    smpar = 4
    smpar_bg = 512
    for j in range(dm_stepnumb):
        acc_dm[j, :] = smooth_edge(acc_dm_file[j, :], smpar) - smooth_edge(acc_dm_file[j, :], smpar_bg)
    mn = np.mean(acc_dm)
    sd = np.std(acc_dm)
    if sd == 0:
        sd = 1.0
    acc_dm_norm = (acc_dm - mn) / sd

    # Use the full time range as the window
    p_viz = np.full(picsize, 200.0, dtype=np.float32)

    repeating_analysis(
        acc_dm_norm, args.dm_pos,
        float(np.min(acc_dm_norm)), float(np.max(acc_dm_norm)),
        p_viz,
    )
    plt.show()


if __name__ == "__main__":
    main()
