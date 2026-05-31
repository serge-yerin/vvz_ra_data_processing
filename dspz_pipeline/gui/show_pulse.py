"""
Individual pulse viewer with interactive DM / time-shift / bandwidth adjustment.

Translated from IDL ``showpulse.pro``.

Displays:
  - The cleaned spectrogram from the .ucd file around the selected time
  - The pulse profile (sum across frequency) with sub-band breakdowns
  - A spectrum of the pulse (S/N vs frequency)
  - Interactive buttons for adjusting DM, time shift, and frequency resolution

Usage::

    python -m dspz_pipeline.gui.show_pulse <ucd_file> <dm> --ns 32768
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("TkAgg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

try:
    import tkinter as tk
    from tkinter import ttk
except ImportError:
    raise SystemExit("tkinter is required for the GUI.")

from dspz_pipeline.config import HEADER_SIZE_BYTES
from dspz_pipeline.io.dmt import compute_dm_delays
from dspz_pipeline.cleaning.robust_stats import erov
from dspz_pipeline.utils import smooth_edge


class ShowPulseApp:
    """Interactive individual pulse viewer.

    Translates IDL ``ShowPulse`` procedure from showpulse.pro.
    """

    def __init__(
        self,
        filename: str,
        dm_const: float,
        dm_pos: int,
        ns: int,
        picsize: int,
        smpar: int,
        acc_dm_norm: np.ndarray | None = None,
        dm_stepnumb: int = 51,
        trans_app=None,
    ):
        # The filename from TransSearch points to the .dmt file;
        # ShowPulse reads the .ucd file (strip .dmt extension)
        self.dmt_filename = filename
        ucd_name = filename
        if ucd_name.endswith(".dmt"):
            ucd_name = ucd_name[:-4]
        self.ucd_filename = ucd_name
        self.dm_const = dm_const
        self.dm_pos = dm_pos
        self.ns = ns
        self.picsize = picsize
        self.smpar = smpar
        self.acc_dm_norm = acc_dm_norm
        self.dm_stepnumb = dm_stepnumb
        self.trans_app = trans_app

        # Parameters from IDL showpulse.pro
        self.time_res = 64 * 8192.0 / 66_000_000.0
        self.wofsg = 4096
        self.nsframe = 50
        self.dm_step = 0.002
        self.sh_t = 0          # time shift in samples
        self.smfreq = 8        # frequency smoothing exponent

        # Compute current DM
        self.dm = dm_const + 0.04 * (dm_pos - 25)

        # Save initial values for reset
        self.dm_init = self.dm
        self.sh_t_init = self.sh_t
        self.smfreq_init = self.smfreq

        # Load data from .ucd file
        self._load_ucd_data()

        # Build GUI
        self.root = tk.Toplevel() if acc_dm_norm is not None else tk.Tk()
        self.root.title("Individual Pulse Viewer")
        self.root.protocol("WM_DELETE_WINDOW", self._on_close)

        self._build_controls()
        self._build_canvas()
        self._update_display()

    def _load_ucd_data(self):
        """Load a window of data from the .ucd file around self.ns."""
        nsf = self.nsframe
        max_nofs = 44000

        # Clamp ns to valid range
        if self.ns < 2 * nsf:
            self.ns = 2 * nsf
        if self.ns >= self.picsize - 2 * nsf - 1:
            self.ns = self.picsize - 2 * nsf - 1

        # Read a chunk from the .ucd file using file offset
        offset_bytes = HEADER_SIZE_BYTES + 4 * self.wofsg * (self.ns - 2 * nsf)
        n_time = min(max_nofs, self.picsize - (self.ns - 2 * nsf))

        self.dat_ucd = np.memmap(
            self.ucd_filename,
            dtype=np.float32,
            mode="r",
            offset=offset_bytes,
            shape=(n_time, self.wofsg),
        ).T.copy()  # shape: (wofsg, n_time)

    # ------------------------------------------------------------------ #
    #  GUI
    # ------------------------------------------------------------------ #

    def _build_controls(self):
        ctrl = ttk.Frame(self.root)
        ctrl.pack(side=tk.TOP, fill=tk.X, padx=4, pady=4)

        ttk.Button(ctrl, text="ESC", command=self._on_close).pack(side=tk.LEFT, padx=2)
        ttk.Button(ctrl, text="Reset", command=self._reset).pack(side=tk.LEFT, padx=2)
        ttk.Button(ctrl, text="Save PNGs", command=self._save_pngs).pack(side=tk.LEFT, padx=2)

        ttk.Label(ctrl, text="  DM:").pack(side=tk.LEFT)
        ttk.Button(ctrl, text="-", width=2, command=lambda: self._adj_dm(-1)).pack(side=tk.LEFT)
        self.lbl_dm = ttk.Label(ctrl, text=f"{self.dm:.4f}", width=10)
        self.lbl_dm.pack(side=tk.LEFT)
        ttk.Button(ctrl, text="+", width=2, command=lambda: self._adj_dm(+1)).pack(side=tk.LEFT)

        ttk.Label(ctrl, text="  Shift:").pack(side=tk.LEFT)
        ttk.Button(ctrl, text="<<", width=3, command=lambda: self._adj_shift(-5)).pack(side=tk.LEFT)
        ttk.Button(ctrl, text="<", width=2, command=lambda: self._adj_shift(-1)).pack(side=tk.LEFT)
        self.lbl_shift = ttk.Label(ctrl, text=str(self.sh_t), width=4)
        self.lbl_shift.pack(side=tk.LEFT)
        ttk.Button(ctrl, text=">", width=2, command=lambda: self._adj_shift(+1)).pack(side=tk.LEFT)
        ttk.Button(ctrl, text=">>", width=3, command=lambda: self._adj_shift(+5)).pack(side=tk.LEFT)

        ttk.Label(ctrl, text="    Subband adjust:").pack(side=tk.LEFT)
        ttk.Button(ctrl, text="Finer", command=lambda: self._adj_smfreq(+1)).pack(side=tk.LEFT, padx=2)
        ttk.Button(ctrl, text="Wider", command=lambda: self._adj_smfreq(-1)).pack(side=tk.LEFT, padx=2)
        self.lbl_band = ttk.Label(ctrl, text="", width=12)
        self.lbl_band.pack(side=tk.LEFT)

        self.highlight_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(
            ctrl, text="Highlight",
            variable=self.highlight_var,
            command=self._update_display,
        ).pack(side=tk.LEFT, padx=6)

        ttk.Label(ctrl, text="  STD pts:").pack(side=tk.LEFT)
        self.std_npts_var = tk.IntVar(value=20)
        ttk.Entry(ctrl, textvariable=self.std_npts_var, width=5).pack(side=tk.LEFT)
        ttk.Button(ctrl, text="Calc STD", command=self._update_display).pack(side=tk.LEFT, padx=2)

    def _build_canvas(self):
        # Main window: 2x2 grid of pulse analysis plots
        self.fig = Figure(figsize=(11, 7), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.root)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        # Secondary window: cleaned data overview
        self.data_win = tk.Toplevel(self.root)
        self.data_win.title("Cleaned data (from .ucd)")
        self.data_win.protocol("WM_DELETE_WINDOW", self._on_close)

        toolbar_frame = ttk.Frame(self.data_win)
        toolbar_frame.pack(side=tk.TOP, fill=tk.X)

        main_frame = ttk.Frame(self.data_win)
        main_frame.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        sliders_frame = ttk.Frame(main_frame)
        sliders_frame.pack(side=tk.LEFT, fill=tk.Y, padx=2, pady=2)

        col_high = ttk.Frame(sliders_frame)
        col_high.pack(side=tk.LEFT, fill=tk.Y, padx=1)
        ttk.Label(col_high, text="vmax").pack(side=tk.TOP)
        self.vmax_var = tk.DoubleVar(value=10.0)
        self.scale_high = tk.Scale(
            col_high, from_=15.0, to=-5.0, resolution=0.1,
            orient=tk.VERTICAL, variable=self.vmax_var,
            command=self._on_vmax_change, showvalue=True,
            length=240, width=10, sliderlength=18,
        )
        self.scale_high.pack(side=tk.TOP, fill=tk.Y, expand=True)

        col_low = ttk.Frame(sliders_frame)
        col_low.pack(side=tk.LEFT, fill=tk.Y, padx=1)
        ttk.Label(col_low, text="vmin").pack(side=tk.TOP)
        self.vmin_var = tk.DoubleVar(value=-1.0)
        self.scale_low = tk.Scale(
            col_low, from_=15.0, to=-5.0, resolution=0.1,
            orient=tk.VERTICAL, variable=self.vmin_var,
            command=self._on_vmin_change, showvalue=True,
            length=240, width=10, sliderlength=18,
        )
        self.scale_low.pack(side=tk.TOP, fill=tk.Y, expand=True)

        self.data_fig = Figure(figsize=(11, 3), dpi=100)
        self.data_canvas = FigureCanvasTkAgg(self.data_fig, master=main_frame)
        self.data_canvas.get_tk_widget().pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

        NavigationToolbar2Tk(self.data_canvas, toolbar_frame).update()

    # ------------------------------------------------------------------ #
    #  Display update
    # ------------------------------------------------------------------ #

    def _update_display(self):
        """Recompute dedispersion and redraw all panels."""
        self.fig.clear()

        nsf = self.nsframe
        wofsg = self.wofsg
        nofch = 2 ** self.smfreq

        # Compute dispersion delays for current DM
        dt = compute_dm_delays(self.dm, 33.0, 16.5, wofsg, 16.5 / wofsg)

        # Build the pulse array: shift each freq channel by the delay
        pulse = np.zeros((wofsg, 2 * nsf + 1), dtype=np.float64)
        for j in range(wofsg):
            stsp = int(round(dt[j] / self.time_res)) + self.sh_t + nsf
            stsp = max(0, min(stsp, self.dat_ucd.shape[1] - 2 * nsf - 1))
            end = stsp + 2 * nsf + 1
            if end <= self.dat_ucd.shape[1]:
                pulse[j, :] = self.dat_ucd[j, stsp:end]

        self._redraw_data_plot()

        # Create ax_img first so ax_spec can share its Y axis
        ax_img = self.fig.add_subplot(2, 2, 2)

        # --- Ax 2: Spectrum of pulse (rotated CCW: frequency on Y, STDs on X) #
        ax_spec = self.fig.add_subplot(2, 2, 1, sharey=ax_img)
        op_spec = np.sum(pulse[:, nsf - 10:nsf + 11], axis=1)
        mt_spec, st_spec = erov(op_spec)
        if st_spec == 0:
            st_spec = 1.0
        rebinned = np.mean(
            ((op_spec - mt_spec) / st_spec).reshape(nofch, wofsg // nofch),
            axis=1,
        ) * np.sqrt(nofch / float(wofsg))
        try:
            n_pts = int(self.std_npts_var.get())
        except (tk.TclError, ValueError):
            n_pts = 20
        n_pts = max(2, min(n_pts, len(rebinned)))
        norm_mean = float(np.mean(rebinned[:n_pts]))
        norm_std = float(np.std(rebinned[:n_pts]))
        if norm_std == 0:
            norm_std = 1.0
        rebinned = (rebinned - norm_mean) / norm_std
        freq_axis = 16.5 + np.arange(nofch) / (nofch - 1.0) * 16.5
        ax_spec.step(rebinned, freq_axis, where="mid")
        ax_spec.invert_xaxis()
        ax_spec.set_xlabel("STDs")
        ax_spec.set_ylabel("Frequency (MHz)")

        if getattr(self, "highlight_var", None) is not None and self.highlight_var.get():
            half_step = (16.5 / (nofch - 1.0)) / 2.0 if nofch > 1 else 0.0
            y_bot = freq_axis[0] - half_step
            y_top = freq_axis[n_pts - 1] + half_step
            ax_spec.axhspan(y_bot, y_top, alpha=0.15, color="yellow", zorder=2)
        band_khz = 33000.0 * wofsg / float(nofch) / 8192.0
        self.lbl_band.config(text=f"{band_khz:.1f} kHz")
        ax_spec.set_title(f"Spectrum of pulse (Bandwidth: {band_khz:.1f} kHz)")

        # --- Ax 3: Pulse image -------------------------------------------- #
        pulse_sm = pulse.copy()
        for j in range(wofsg):
            pulse_sm[j, :] = smooth_edge(pulse[j, :], self.smpar)

        # --- Ax 4: Pulse profile (smoothed, total) ----------------------- #
        ax_prof = self.fig.add_subplot(2, 2, 3)
        profile = np.sum(pulse_sm, axis=0)
        bg_mean = np.mean(np.sum(pulse_sm[:, :31], axis=0))
        bg_std = np.std(np.sum(pulse_sm[:, :31], axis=0))
        if bg_std == 0:
            bg_std = 1.0
        profile_sn = (profile - bg_mean) / bg_std
        ax_prof.plot(profile_sn)
        ax_prof.set_xlabel("Time sample")
        ax_prof.set_ylabel("S/N")
        ax_prof.set_title(f"Pulse profile  DM={self.dm:.4f}")

        ax_img.imshow(
            pulse_sm,
            aspect="auto", origin="lower", cmap="gray_r",
            interpolation="nearest",
            extent=[0, pulse_sm.shape[1], 16.5, 33.0],
        )
        ax_img.set_xlabel("Time sample")
        ax_img.set_ylabel("Frequency (MHz)")
        ax_img.set_title(f"Dedispersed pulse (Spectrum # {self.ns} in file)")

        if getattr(self, "highlight_var", None) is not None and self.highlight_var.get():
            x_mid = pulse_sm.shape[1] / 2.0
            x_half = 0.10 * pulse_sm.shape[1]
            ax_img.axvspan(
                x_mid - x_half, x_mid + x_half,
                alpha=0.15, color="yellow", zorder=2,
            )

        # --- Ax 5: Sub-band profiles ------------------------------------- #
        ax_sub = self.fig.add_subplot(2, 2, 4)
        sub_ranges = [(0, 1024), (1024, 2048), (2048, 3072), (3072, 4096)]
        freq_labels = ["16.50-20.63", "20.63-24.75", "24.75-28.88", "28.88-33.00"]
        sub_colors = ["tab:blue", "tab:green", "tab:orange", "tab:red"]
        for (lo, hi), lbl, color in zip(sub_ranges, freq_labels, sub_colors):
            hi = min(hi, wofsg)
            sub = np.sum(pulse_sm[lo:hi, :], axis=0)
            sub_bg_m = np.mean(np.sum(pulse_sm[lo:hi, :31], axis=0))
            sub_bg_s = np.std(np.sum(pulse_sm[lo:hi, :31], axis=0))
            if sub_bg_s == 0:
                sub_bg_s = 1.0
            ax_sub.plot((sub - sub_bg_m) / sub_bg_s, label=f"{lbl} MHz", color=color)
        ax_sub.legend(fontsize=7)
        ax_sub.set_xlabel("Time sample")
        ax_sub.set_ylabel("S/N")
        ax_sub.set_title(f"Sub-band profiles (Shift: {self.sh_t})")

        self.fig.tight_layout()
        self.canvas.draw()

    # ------------------------------------------------------------------ #
    #  Control callbacks
    # ------------------------------------------------------------------ #

    def _redraw_data_plot(self):
        """Redraw the cleaned-data overview using current vmin/vmax sliders."""
        if not hasattr(self, "_data_norm"):
            n_show = min(44000, self.dat_ucd.shape[1])
            data_slice = self.dat_ucd[:, :n_show]
            data_mean = float(np.mean(data_slice))
            data_std = float(np.std(data_slice))
            if data_std == 0:
                data_std = 1.0
            self._data_norm = (data_slice - data_mean) / data_std
            self._data_n_show = n_show

        self.data_fig.clear()
        ax_data = self.data_fig.add_subplot(1, 1, 1)
        ax_data.imshow(
            self._data_norm,
            aspect="auto", origin="lower", cmap="gray_r",
            interpolation="nearest",
            extent=[0, self._data_n_show, 16.5, 33.0],
            vmin=self.vmin_var.get(), vmax=self.vmax_var.get(),
        )
        ax_data.set_title("Cleaned data (from .ucd)")
        ax_data.set_xlabel("Time sample")
        ax_data.set_ylabel("Frequency (MHz)")
        self.data_fig.tight_layout()
        self.data_canvas.draw()

    def _on_vmin_change(self, value):
        v = float(value)
        vmax = self.vmax_var.get()
        if v >= vmax:
            self.vmin_var.set(round(vmax - 0.1, 1))
            return
        self._redraw_data_plot()

    def _on_vmax_change(self, value):
        v = float(value)
        vmin = self.vmin_var.get()
        if v <= vmin:
            self.vmax_var.set(round(vmin + 0.1, 1))
            return
        self._redraw_data_plot()

    def _reset(self):
        self.dm = self.dm_init
        self.sh_t = self.sh_t_init
        self.smfreq = self.smfreq_init
        self.lbl_dm.config(text=f"{self.dm:.4f}")
        self.lbl_shift.config(text=str(self.sh_t))
        self._update_display()

    def _adj_dm(self, direction):
        self.dm += direction * self.dm_step
        self.lbl_dm.config(text=f"{self.dm:.4f}")
        self._update_display()

    def _adj_shift(self, direction):
        if -self.nsframe < self.sh_t - direction < self.nsframe:
            self.sh_t -= direction
        self.lbl_shift.config(text=str(self.sh_t))
        self._update_display()

    def _adj_smfreq(self, direction):
        self.smfreq = max(2, min(10, self.smfreq + direction))
        self._update_display()

    @staticmethod
    def _save_figure_fixed(fig, path, target_width=16.0, dpi=200):
        """Save a figure at a fixed width in inches (preserving aspect ratio),
        independent of the current window size."""
        orig_w, orig_h = fig.get_size_inches()
        aspect = orig_h / orig_w
        fig.set_size_inches(target_width, target_width * aspect)
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        fig.set_size_inches(orig_w, orig_h)
        fig.canvas.draw_idle()

    def _save_pngs(self):
        ucd_path = Path(self.ucd_filename)
        folder = ucd_path.parent / f"{ucd_path.name}_individual_pulse_ns={self.ns}"
        folder.mkdir(exist_ok=True)

        self._save_figure_fixed(self.fig,      folder / "individual_pulse_viewer.png")
        print(f"Saved: {folder / 'individual_pulse_viewer.png'}")

        self._save_figure_fixed(self.data_fig, folder / "cleaned_data.png", target_width=32.0)
        print(f"Saved: {folder / 'cleaned_data.png'}")

        if self.trans_app is not None:
            self._save_figure_fixed(self.trans_app.fig, folder / "transient_search.png", target_width=32.0)
            print(f"Saved: {folder / 'transient_search.png'}")

            fig_sel = getattr(self.trans_app, "fig_sel", None)
            if fig_sel is not None:
                self._save_figure_fixed(fig_sel, folder / "pulse_selection.png")
                print(f"Saved: {folder / 'pulse_selection.png'}")

        print(f"All PNGs saved to: {folder}")

    def _on_close(self):
        self.data_win.destroy()
        self.root.destroy()

    def run(self):
        """Start the tkinter main loop (only when run standalone)."""
        if isinstance(self.root, tk.Tk):
            self.root.mainloop()


def show_pulse_gui(
    filename: str,
    dm_const: float,
    dm_pos: int,
    ns: int,
    picsize: int,
    smpar: int,
    acc_dm_norm: np.ndarray | None = None,
    dm_stepnumb: int = 51,
    trans_app=None,
) -> None:
    """Launch the individual pulse viewer."""
    app = ShowPulseApp(
        filename, dm_const, dm_pos, ns, picsize, smpar,
        acc_dm_norm, dm_stepnumb, trans_app,
    )
    if acc_dm_norm is None:
        app.run()


def main() -> None:
    """CLI entry point for ``python -m dspz_pipeline.gui.show_pulse``."""
    parser = argparse.ArgumentParser(description="Individual pulse viewer")
    parser.add_argument("ucd_file", help="Path to the .ucd data file")
    parser.add_argument("dm", type=float, help="Central DM in pc/cm^3")
    parser.add_argument("--dm-pos", type=int, default=25, help="DM step index (default: 25)")
    parser.add_argument("--ns", type=int, default=1000, help="Time sample to inspect (default: 1000)")
    parser.add_argument("--picsize", type=int, default=None, help="Total time samples (auto-detected if omitted)")
    parser.add_argument("--smpar", type=int, default=4, help="Smoothing parameter (default: 4)")
    args = parser.parse_args()

    # Auto-detect picsize from file
    if args.picsize is None:
        import os
        fsize = os.path.getsize(args.ucd_file)
        args.picsize = (fsize - HEADER_SIZE_BYTES) // (4 * 4096)

    # For standalone run, create a fake .dmt filename so the .ucd path derivation works
    fake_dmt = args.ucd_file + ".dmt"

    app = ShowPulseApp(
        fake_dmt, args.dm, args.dm_pos, args.ns, args.picsize, args.smpar,
    )
    app.run()


if __name__ == "__main__":
    main()
