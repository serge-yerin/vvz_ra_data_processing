"""
Interactive transient / pulsar search GUI.

Translated from IDL ``TransSearch.pro``.

Displays dedispersed data (.dmt file) as 16 horizontal strips of
DM-vs-time spectrograms. The user can:
  - Adjust high-pass smoothing (smpar) and low-pass background (smpar_b)
  - Shift the DM offset
  - Toggle individual pulse analysis (IND) and repeating analysis (REP)
  - Select a time window and number of parts for FFT analysis
  - Click on the spectrogram to inspect pulses

Usage::

    python -m dspz_pipeline.gui.trans_search path/to/file.dmt 12.872
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("TkAgg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

try:
    import tkinter as tk
    from tkinter import ttk
except ImportError:
    raise SystemExit("tkinter is required for the GUI — install the tk package.")

from dspz_pipeline.config import DEFAULT_DM_CONST
from dspz_pipeline.io.dmt import read_dmt
from dspz_pipeline.utils import smooth_edge


def _compute_accDM(acc_dm_file, dm_stepnumb, picsize, smpar, smpar_background):
    """Apply high-pass and low-pass smoothing to produce accDM.

    IDL: accDM[j,*] = smooth(accDMfile[j,*], smpar, /EDG)
                     - smooth(accDMfile[j,*], smpar_background, /EDG)
    """
    acc_dm = np.zeros_like(acc_dm_file)
    for j in range(dm_stepnumb):
        row = acc_dm_file[j, :]
        acc_dm[j, :] = smooth_edge(row, smpar) - smooth_edge(row, smpar_background)
    return acc_dm


class TransSearchApp:
    """Tkinter application for interactive transient/pulsar search.

    Translates the interactive loop and GUI from IDL ``TransSearch.pro``.
    """

    KADR = 4096     # samples per panel strip
    N_PANELS = 16   # number of horizontal strips
    NSFRAME = 50    # half-width of pulse selection window

    def __init__(self, filename: str, dm_const: float):
        self.filename = filename
        self.dm_const = dm_const

        # --- Load data ---------------------------------------------------- #
        self.acc_dm_file, self.dm_stepnumb, self.picsize = read_dmt(filename)
        print(f"Loaded: {self.dm_stepnumb} DM steps x {self.picsize} samples\n")

        # --- State variables (matching IDL globals) ----------------------- #
        self.smpar = 4
        self.smpar_pow = 9
        self.smpar_background = 2 ** self.smpar_pow
        self.dm_pos = 25          # index into DM steps (0..dm_stepnumb-1)
        self.ind_mode = False     # individual pulse analysis toggle
        self.rep_mode = False     # repeating analysis toggle
        self.nofp_pow = 0
        self.nofp = 1             # number of parts (2^nofp_pow)
        self.p_num = 0            # current part number

        self.minscl = None
        self.maxscl = None

        # Computed arrays
        self.acc_dm = None
        self.acc_dm_norm = None
        self.maxnorm = 0.0
        self.minnorm = 0.0

        # --- Build GUI ---------------------------------------------------- #
        self.root = tk.Tk()
        self.root.title("Transient / Pulsar Search")
        self.root.protocol("WM_DELETE_WINDOW", self._on_close)

        self._build_controls()
        self._build_canvas()

        # Initial computation and display
        self._recompute_smoothing()
        self._redraw()

    # ------------------------------------------------------------------ #
    #  GUI construction
    # ------------------------------------------------------------------ #

    def _build_controls(self):
        """Create the control panel with buttons matching the IDL interface."""
        ctrl = ttk.Frame(self.root)
        ctrl.pack(side=tk.TOP, fill=tk.X, padx=4, pady=4)

        # Row 1: smoothing params
        r1 = ttk.Frame(ctrl)
        r1.pack(fill=tk.X, pady=2)

        ttk.Button(r1, text="Close", command=self._on_close).pack(side=tk.LEFT, padx=4)

        ttk.Label(r1, text="  smpar:").pack(side=tk.LEFT)
        ttk.Button(r1, text="<", width=2, command=lambda: self._adj_smpar(-1)).pack(side=tk.LEFT)
        self.lbl_smpar = ttk.Label(r1, text=str(self.smpar), width=3)
        self.lbl_smpar.pack(side=tk.LEFT)
        ttk.Button(r1, text=">", width=2, command=lambda: self._adj_smpar(+1)).pack(side=tk.LEFT)

        ttk.Label(r1, text="  smpar_b:").pack(side=tk.LEFT)
        ttk.Button(r1, text="<", width=2, command=lambda: self._adj_smpar_b(-1)).pack(side=tk.LEFT)
        self.lbl_smpar_b = ttk.Label(r1, text=str(self.smpar_background), width=5)
        self.lbl_smpar_b.pack(side=tk.LEFT)
        ttk.Button(r1, text=">", width=2, command=lambda: self._adj_smpar_b(+1)).pack(side=tk.LEFT)

        ttk.Label(r1, text="  DM step:").pack(side=tk.LEFT)
        ttk.Button(r1, text="<", width=2, command=lambda: self._adj_dm(-1)).pack(side=tk.LEFT)
        self.lbl_dm = ttk.Label(r1, text=self._dm_label(), width=6)
        self.lbl_dm.pack(side=tk.LEFT)
        ttk.Button(r1, text=">", width=2, command=lambda: self._adj_dm(+1)).pack(side=tk.LEFT)

        self.btn_ind = ttk.Button(r1, text="Individual", width=14, command=self._toggle_ind)
        self.btn_ind.pack(side=tk.LEFT, padx=4)
        self.btn_rep = ttk.Button(r1, text="Repetitive", width=14, command=self._toggle_rep)
        self.btn_rep.pack(side=tk.LEFT, padx=4)

        ttk.Label(r1, text="  parts:").pack(side=tk.LEFT)
        ttk.Button(r1, text="<", width=2, command=lambda: self._adj_nofp(-1)).pack(side=tk.LEFT)
        self.lbl_nofp = ttk.Label(r1, text=str(self.nofp), width=3)
        self.lbl_nofp.pack(side=tk.LEFT)
        ttk.Button(r1, text=">", width=2, command=lambda: self._adj_nofp(+1)).pack(side=tk.LEFT)

        ttk.Label(r1, text="  N of parts:").pack(side=tk.LEFT)
        ttk.Button(r1, text="<", width=2, command=lambda: self._adj_pnum(-1)).pack(side=tk.LEFT)
        self.lbl_pnum = ttk.Label(r1, text=str(self.p_num), width=3)
        self.lbl_pnum.pack(side=tk.LEFT)
        ttk.Button(r1, text=">", width=2, command=lambda: self._adj_pnum(+1)).pack(side=tk.LEFT)

        # Row 2: scale sliders
        r2 = ttk.Frame(ctrl)
        r2.pack(fill=tk.X, pady=2)
        ttk.Label(r2, text="Min scale:").pack(side=tk.LEFT)
        self.scl_min = tk.Scale(r2, from_=-50, to=50, resolution=0.1,
                                orient=tk.HORIZONTAL, length=200,
                                command=self._on_scale_change)
        self.scl_min.pack(side=tk.LEFT, padx=4)
        ttk.Label(r2, text="Max scale:").pack(side=tk.LEFT)
        self.scl_max = tk.Scale(r2, from_=-50, to=50, resolution=0.1,
                                orient=tk.HORIZONTAL, length=200,
                                command=self._on_scale_change)
        self.scl_max.pack(side=tk.LEFT, padx=4)

    def _build_canvas(self):
        """Create the matplotlib canvas for the 16-strip spectrogram."""
        self.fig = Figure(figsize=(11, 9), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.root)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        self.canvas.mpl_connect("button_press_event", self._on_click)

    # ------------------------------------------------------------------ #
    #  Computation
    # ------------------------------------------------------------------ #

    def _recompute_smoothing(self):
        """Recompute accDM and normalised array after parameter change."""
        self.acc_dm = _compute_accDM(
            self.acc_dm_file, self.dm_stepnumb, self.picsize,
            self.smpar, self.smpar_background,
        )
        # Normalise: (accDM - mean) / stddev
        mn = np.mean(self.acc_dm)
        sd = np.std(self.acc_dm)
        if sd == 0:
            sd = 1.0
        self.acc_dm_norm = (self.acc_dm - mn) / sd
        self.maxnorm = float(np.max(self.acc_dm_norm))
        self.minnorm = float(np.min(self.acc_dm_norm))

        if self.minscl is None:
            self.minscl = self.minnorm
            self.maxscl = self.maxnorm

        # Update scale sliders range and value
        self.scl_min.config(from_=self.minnorm, to=self.maxnorm)
        self.scl_max.config(from_=self.minnorm, to=self.maxnorm)
        self.scl_min.set(self.minscl)
        self.scl_max.set(self.maxscl)

    # ------------------------------------------------------------------ #
    #  Drawing
    # ------------------------------------------------------------------ #

    def _redraw(self):
        """Redraw the 16-strip spectrogram display."""
        self.fig.clear()

        kadr = self.KADR
        n_panels = self.N_PANELS
        minscl = self.minscl
        maxscl = self.maxscl

        # Compute part visualization bar (pViz equivalent)
        end_p = self.p_num * n_panels + (n_panels * n_panels) // self.nofp if self.nofp > 0 else n_panels * n_panels
        beg_p = self.p_num * n_panels
        if end_p > n_panels * n_panels:
            del_p = end_p - n_panels * n_panels
            end_p = n_panels * n_panels
            beg_p = self.p_num * n_panels - del_p
        else:
            del_p = 0

        # Store pViz for RepeatingAnalysis
        self.p_viz = np.zeros(n_panels * kadr, dtype=np.float32)
        viz_start = beg_p * kadr // n_panels
        viz_end = end_p * kadr // n_panels
        viz_start = max(0, min(viz_start, self.picsize))
        viz_end = max(0, min(viz_end, self.picsize))
        if viz_end > viz_start:
            self.p_viz[viz_start:viz_end] = 200.0

        for np_idx in range(n_panels):
            ax = self.fig.add_subplot(n_panels, 1, n_panels - np_idx)

            t_start = kadr * np_idx
            t_end = (np_idx + 1) * kadr
            if t_end > self.picsize:
                t_end = self.picsize

            chunk = self.acc_dm_norm[:, t_start:t_end]

            # Clip to scale range (IDL: >minscl<maxscl)
            display = np.clip(chunk, minscl, maxscl)

            ax.imshow(
                display,
                aspect="auto",
                origin="lower",
                cmap="gray_r",
                vmin=minscl,
                vmax=maxscl,
                extent=[t_start, t_end, 0, self.dm_stepnumb],
                interpolation="nearest",
            )

            # Mark the active part with a coloured overlay on the left
            if viz_start <= t_end and viz_end >= t_start:
                ax.axvspan(
                    max(t_start, viz_start), min(t_end, viz_end),
                    alpha=0.15, color="yellow", zorder=2,
                )

            ax.set_ylabel(f"P{np_idx}", fontsize=7, rotation=0, labelpad=15)
            ax.tick_params(labelsize=6)
            if np_idx > 0:
                ax.set_xticklabels([])

        self.fig.suptitle(
            f"smpar={self.smpar}  smpar_b={self.smpar_background}  "
            f"DM step={self.dm_pos - 25:+d}  "
            f"DM={self.dm_const + (self.dm_pos - 25) * 0.04:.3f}  "
            f"IND={'ON' if self.ind_mode else 'off'}  "
            f"REP={'ON' if self.rep_mode else 'off'}",
            fontsize=9,
        )
        self.fig.tight_layout(rect=[0, 0, 1, 0.97])
        self.canvas.draw()

    # ------------------------------------------------------------------ #
    #  Control callbacks
    # ------------------------------------------------------------------ #

    def _dm_label(self) -> str:
        sign = "+" if (self.dm_pos - 25) >= 0 else ""
        return f"{sign}{self.dm_pos - 25}"

    def _adj_smpar(self, delta):
        self.smpar = max(1, min(8, self.smpar + delta))
        self.lbl_smpar.config(text=str(self.smpar))
        self._recompute_smoothing()
        self._redraw()

    def _adj_smpar_b(self, delta):
        self.smpar_pow = max(5, min(12, self.smpar_pow + delta))
        self.smpar_background = 2 ** self.smpar_pow
        self.lbl_smpar_b.config(text=str(self.smpar_background))
        self._recompute_smoothing()
        self._redraw()

    def _adj_dm(self, delta):
        self.dm_pos = max(0, min(self.dm_stepnumb - 1, self.dm_pos + delta))
        self.lbl_dm.config(text=self._dm_label())
        self._redraw()

    def _toggle_ind(self):
        self.ind_mode = not self.ind_mode
        self.btn_ind.config(text=f"Individual {'ON' if self.ind_mode else 'OFF'}")
        print(f"IND = {self.ind_mode}")

    def _toggle_rep(self):
        self.rep_mode = not self.rep_mode
        self.btn_rep.config(text=f"Repetitive {'ON' if self.rep_mode else 'OFF'}")
        print(f"REP = {self.rep_mode}")

    def _adj_nofp(self, delta):
        self.nofp_pow = max(0, min(4, self.nofp_pow + delta))
        self.nofp = 2 ** self.nofp_pow
        self.lbl_nofp.config(text=str(self.nofp))
        self._redraw()

    def _adj_pnum(self, delta):
        self.p_num = max(0, min(15, self.p_num + delta))
        self.lbl_pnum.config(text=str(self.p_num))
        self._redraw()

    def _on_scale_change(self, _val=None):
        self.minscl = float(self.scl_min.get())
        self.maxscl = float(self.scl_max.get())
        self._redraw()

    # ------------------------------------------------------------------ #
    #  Click handling — select a point in the spectrogram
    # ------------------------------------------------------------------ #

    def _on_click(self, event):
        """Handle click on the spectrogram to select a time sample."""
        if event.inaxes is None:
            return

        x_data = event.xdata
        if x_data is None:
            return

        ns = int(round(x_data))
        ns = max(self.NSFRAME, min(ns, self.picsize - self.NSFRAME - 1))

        print(f"Selected NS={ns}")

        # Show pulse selection window (surface + image)
        self._show_pulse_selection(ns)

        # Dispatch to sub-analyses if toggled
        if self.rep_mode:
            from dspz_pipeline.gui.repeating_analysis import repeating_analysis
            repeating_analysis(
                self.acc_dm_norm, self.dm_pos, self.minscl, self.maxscl,
                self.p_viz,
            )

        if self.ind_mode:
            from dspz_pipeline.gui.show_pulse import show_pulse_gui
            show_pulse_gui(
                self.filename, self.dm_const, self.dm_pos, ns,
                self.picsize, self.smpar, self.acc_dm_norm, self.dm_stepnumb,
            )

    def _show_pulse_selection(self, ns: int):
        """Show the pulse-selection window (3-D surface + 2-D image)."""
        nsf = self.NSFRAME
        chunk = self.acc_dm_norm[:, ns - nsf: ns + nsf + 1]
        chunk_clipped = np.clip(chunk, self.minscl, self.maxscl)

        fig_sel, (ax_img, ax_surf) = plt.subplots(
            1, 2, figsize=(8, 4),
            subplot_kw={"projection": None},
        )
        # 2-D image
        ax_img.imshow(
            chunk_clipped,
            aspect="auto", origin="lower", cmap="gray_r",
            extent=[ns - nsf, ns + nsf, 0, self.dm_stepnumb],
        )
        ax_img.set_xlabel("Time sample")
        ax_img.set_ylabel("DM step")
        ax_img.set_title("Pulse selection")

        # 3-D surface
        ax_surf.remove()
        ax_surf = fig_sel.add_subplot(1, 2, 2, projection="3d")
        t_axis = np.arange(2 * nsf + 1)
        dm_axis = np.arange(self.dm_stepnumb)
        T, D = np.meshgrid(t_axis, dm_axis)
        ax_surf.plot_surface(T, D, chunk_clipped, cmap="viridis", linewidth=0, antialiased=False)
        ax_surf.set_xlabel("Time")
        ax_surf.set_ylabel("DM step")
        ax_surf.set_zlabel("S/N")
        ax_surf.set_title("3D view")

        fig_sel.tight_layout()
        fig_sel.show()

    # ------------------------------------------------------------------ #
    #  Main loop
    # ------------------------------------------------------------------ #

    def _on_close(self):
        self.root.destroy()

    def run(self):
        """Start the tkinter main loop."""
        self.root.mainloop()


def trans_search_gui(filename: str, dm_const: float) -> None:
    """Entry point to launch the TransSearch GUI.

    Parameters
    ----------
    filename : str
        Path to the ``.dmt`` file.
    dm_const : float
        Central dispersion measure in pc/cm^3.
    """
    app = TransSearchApp(filename, dm_const)
    app.run()


def main() -> None:
    """CLI entry point for ``python -m dspz_pipeline.gui.trans_search``."""
    parser = argparse.ArgumentParser(description="Interactive transient/pulsar search")
    parser.add_argument("dmt_file", help="Path to the .dmt data file")
    parser.add_argument("dm", type=float, help="Central DM in pc/cm^3")
    args = parser.parse_args()
    trans_search_gui(args.dmt_file, args.dm)


if __name__ == "__main__":
    main()
