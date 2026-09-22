"""
Main pipeline script: read .jds files, clean, dedisperse, launch analysis.

Translated from IDL ``process_survey.pro``.

Pipeline steps:
  1. Select one or more .jds data files.
  2. Read the DSPZ binary header and extract instrument parameters.
  3. For each file, for each frame (1024 spectra):
     a. Decode the custom mantissa/exponent format to float64.
     b. Apply RFI cleaning (flatten + SumThreshold + patrol).
     c. Write cleaned data to a ``.ucd`` output file.
  4. Run incoherent dedispersion (:func:`~dspz_pipeline.analysis.dedispersion.ind_search`)
     to produce a ``.dmt`` file.
  5. Optionally launch the interactive TransSearch GUI.

Usage::

    python -m dspz_pipeline.process_survey --files _data/*.jds --dm 12.872

"""

from __future__ import annotations

import argparse
import os
import sys
from collections import deque
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from tqdm import tqdm

import numpy as np

from dspz_pipeline.config import (
    DEFAULT_DM_CONST,
    DEFAULT_MODE,
    DEFAULT_NOFS,
    DEFAULT_PERIOD,
    DEFAULT_PULSAR_LABEL,
    SHIFT_AB,
)
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg

from dspz_pipeline.analysis.dedispersion import ind_search
from dspz_pipeline.io.jds_reader import JdsFile, decode_raw_frame, write_ucd_header
from dspz_pipeline.cleaning.rfi_cleaning import adr_cleaning


def default_workers() -> int:
    """Leave one core for the reader/writer process; cap at 8."""
    return max(1, min(8, (os.cpu_count() or 2) - 1))


def clean_frame_task(
    raw_bytes: bytes,
    wofsg: int,
    nofs: int,
    avrs: int,
    mode: int,
    png_args: tuple | None,
) -> bytes:
    """Decode + clean one frame; optionally save its PNG.

    Runs in a worker process (frames are independent) and returns the
    cleaned frame as float32 bytes in Fortran order, ready to append to
    the .ucd file.
    """
    imdat = decode_raw_frame(raw_bytes, wofsg, nofs, avrs, mode)
    imdat, mask = adr_cleaning(imdat, pat=True, nar=True, wid=True)
    if png_args is not None:
        ucd_path, frame_num, total_frames = png_args
        save_frame_png(imdat, mask, Path(ucd_path), frame_num, total_frames)
    return imdat.astype(np.float32).tobytes(order="F")


def save_frame_png(imdat, mask, ucd_path: Path, frame_num: int, total_frames: int) -> None:
    """Save imdat and mask side-by-side as a PNG next to the .ucd file.

    Parameters
    ----------
    imdat : ndarray
        Cleaned data array for this frame.
    mask : ndarray
        RFI mask array for this frame.
    ucd_path : Path
        Path to the output .ucd file (used to derive folder and filename).
    frame_num : int
        1-based frame number (i + 1).
    total_frames : int
        Total number of frames across all JDS files (for zero-padding and suptitle).
    """
    out_dir = ucd_path.parent / ucd_path.stem
    out_dir.mkdir(exist_ok=True)

    n_digits = len(str(total_frames))
    fname = f"{ucd_path.stem}_{frame_num:0{n_digits}d}.png"

    fig = Figure(figsize=(12, 5))
    FigureCanvasAgg(fig)

    fig.suptitle(f"{ucd_path.stem}  frame {frame_num} of {total_frames}", fontsize=10)

    ax_data = fig.add_subplot(1, 2, 1)
    ax_data.imshow(imdat, aspect="auto", cmap="Greys", origin="lower")
    ax_data.set_title("Data")

    ax_mask = fig.add_subplot(1, 2, 2)
    ax_mask.imshow(mask, aspect="auto", cmap="Greys_r", origin="lower")
    ax_mask.set_title("Mask")

    fig.tight_layout()
    fig.savefig(out_dir / fname, dpi=150)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    p = argparse.ArgumentParser(
        description="DSPZ pipeline: .jds \u2192 clean \u2192 .ucd \u2192 dedisperse \u2192 .dmt",
    )
    p.add_argument(
        "--indir", default=".",
        help="Directory containing the .jds input files (default: current directory).",
    )
    p.add_argument(
        "--files", nargs="+", required=True,
        help="One or more .jds file names (relative to --indir) to process (in order).",
    )
    p.add_argument(
        "--outdir", default="_output",
        help="Directory for output files (default: _output/).",
    )
    p.add_argument(
        "--dm", type=float, default=DEFAULT_DM_CONST,
        help=f"Central dispersion measure in pc/cm^3 (default: {DEFAULT_DM_CONST}).",
    )
    p.add_argument(
        "--period", type=float, default=DEFAULT_PERIOD,
        help=f"Pulsar period in seconds (default: {DEFAULT_PERIOD}).",
    )
    p.add_argument(
        "--label", default=DEFAULT_PULSAR_LABEL,
        help=f"Pulsar / source label for filenames (default: {DEFAULT_PULSAR_LABEL}).",
    )
    p.add_argument(
        "--mode", type=int, default=DEFAULT_MODE, choices=[0, 1, 2],
        help="Channel mode: 0=ch0-ch1, 1=ch0, 2=ch1 (default: 1).",
    )
    p.add_argument(
        "--nofs", type=int, default=DEFAULT_NOFS,
        help=f"Spectra per frame (default: {DEFAULT_NOFS}).",
    )
    p.add_argument(
        "--no-gui", action="store_true",
        help="Skip launching the interactive TransSearch GUI.",
    )
    p.add_argument(
        "--save_cleaning_mask", action="store_true",
        help="Save PNG images of the cleaned data and RFI mask for each frame.",
    )
    p.add_argument(
        "--workers", type=int, default=default_workers(),
        help="Worker processes for RFI cleaning; 1 = sequential in this process "
             f"(default: {default_workers()}).",
    )
    return p.parse_args(argv)


def run_pipeline(args: argparse.Namespace) -> None:
    """Execute the full processing pipeline."""
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    indir = Path(args.indir)
    jds_files = [indir / f for f in args.files]
    shortnames = [f.name for f in jds_files]
    n_in_list = len(jds_files)

    print(f"\nProcessing {n_in_list} file(s): {[s for s in shortnames]}")
    print(f"\nDM = {args.dm} pc/cm^3,  label = {args.label}")

    # ---- Read header from first file to set up output ---- #
    with JdsFile(jds_files[0], nofs=args.nofs) as jds0:
        first_header = jds0.header

    # Total frame count across all JDS files (progress bar, PNG naming/padding)
    total_frames = 0
    for jds_path in jds_files:
        with JdsFile(jds_path, nofs=args.nofs) as jds:
            total_frames += jds.nframe
    if args.save_cleaning_mask:
        print(f"\nImages of data and RFI mask will be saved for each frame. !!! This may take additional time !!!\n")

    # Output name: Cleaned_<label>_<first input file>.ucd
    cleaned_filename = outdir / f"Cleaned_{args.label}_{shortnames[0]}.ucd"

    # Write .ucd header (copy of .jds header with nofs patched)
    write_ucd_header(cleaned_filename, first_header, args.nofs)

    workers = max(1, args.workers)
    print(f"\nRFI cleaning with {workers} worker process(es)")

    # Frames are cleaned in worker processes; the main process reads raw
    # frames, keeps a bounded number in flight, and appends the results to
    # the .ucd strictly in frame order.
    pbar = tqdm(total=total_frames, desc="  Cleaning", unit="frame", dynamic_ncols=True)
    # Plain print() would land in the middle of the progress bar (the bar
    # redraws on stderr, print writes to stdout, and the two interleave);
    # tqdm.write clears the bar, prints the line, then redraws the bar.
    say = pbar.write
    pool = ProcessPoolExecutor(max_workers=workers) if workers > 1 else None
    max_inflight = 2 * workers
    pending: deque = deque()
    global_frame = 0

    def submit(task_args):
        """Return a zero-arg callable that yields the cleaned frame bytes."""
        if pool is None:
            result = clean_frame_task(*task_args)
            return lambda: result
        return pool.submit(clean_frame_task, *task_args).result

    with open(cleaned_filename, "ab") as ucd_fh:
        try:
            for n_file in range(n_in_list):
                jds_path = jds_files[n_file]
                say(f"\n--- File {n_file + 1} of {n_in_list}: {jds_path.name} ---")

                with JdsFile(jds_path, nofs=args.nofs) as jds:
                    hdr = jds.header
                    say(f"  Name:  {hdr.sname}")
                    say(f"  Local: {hdr.stime}")
                    say(f"  UTC:   {hdr.sgmtt}")
                    say(f"  Mode:  {'waveform' if hdr.data_mode == 0 else 'spectra' if hdr.data_mode == 1 else 'correlation'}")
                    say(f"  Fmin: {hdr.fmin_mhz:.1f} MHz,  Fmax: {hdr.fmax_mhz:.1f} MHz,  "
                        f"wofsg: {hdr.wofsg},  avrs: {hdr.avrs},  "
                        f"Time resolution: {hdr.time_res_s * 1000:.3f} ms")
                    say(f"  Total number of frames in file: {jds.nframe}")

                    for _i, raw_bytes in jds.raw_frames():
                        global_frame += 1
                        png_args = (
                            (str(cleaned_filename), global_frame, total_frames)
                            if args.save_cleaning_mask else None
                        )
                        pending.append(submit(
                            (raw_bytes, hdr.wofsg, args.nofs, hdr.avrs, args.mode, png_args)
                        ))
                        while len(pending) >= max_inflight:
                            ucd_fh.write(pending.popleft()())
                            pbar.update(1)

            while pending:
                ucd_fh.write(pending.popleft()())
                pbar.update(1)
        finally:
            pbar.close()
            if pool is not None:
                pool.shutdown(wait=True, cancel_futures=True)

    print(f"\nCleaned data written to: {cleaned_filename}")

    # ---- Dedispersion ---- #
    print("\n\n--- Running incoherent dedispersion ---\n")
    dmt_path = ind_search(
        cleaned_filename,
        args.dm,
        fmax_mhz=first_header.fmax_mhz,
        fmin_mhz=first_header.fmin_mhz,
    )
    print(f"Dedispersed data written to: {dmt_path}")

    # ---- Interactive analysis (optional) ---- #
    if not args.no_gui:
        print("\n\n--- Launching TransSearch GUI ---\n")
        from dspz_pipeline.gui.trans_search import trans_search_gui
        trans_search_gui(str(dmt_path), args.dm)


def main() -> None:
    """Entry point for ``python -m dspz_pipeline.process_survey``."""
    args = parse_args()
    run_pipeline(args)


if __name__ == "__main__":
    main()
