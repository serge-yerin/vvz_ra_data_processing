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
import sys
import time
from pathlib import Path

import numpy as np

from dspz_pipeline.config import (
    DEFAULT_DM_CONST,
    DEFAULT_MODE,
    DEFAULT_NOFS,
    DEFAULT_PERIOD,
    DEFAULT_PULSAR_LABEL,
    SHIFT_AB,
)
from dspz_pipeline.analysis.dedispersion import ind_search
from dspz_pipeline.io.jds_reader import JdsFile, write_ucd_header
from dspz_pipeline.cleaning.rfi_cleaning import adr_cleaning


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    p = argparse.ArgumentParser(
        description="DSPZ pipeline: .jds \u2192 clean \u2192 .ucd \u2192 dedisperse \u2192 .dmt",
    )
    p.add_argument(
        "--files", nargs="+", required=True,
        help="One or more .jds data files to process (in order).",
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
    return p.parse_args(argv)


def run_pipeline(args: argparse.Namespace) -> None:
    """Execute the full processing pipeline."""
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    jds_files = [Path(f) for f in args.files]
    shortnames = [f.name for f in jds_files]
    n_in_list = len(jds_files)

    print(f"\nProcessing {n_in_list} file(s): {[s for s in shortnames]}")
    print(f"\nDM = {args.dm} pc/cm^3,  label = {args.label}")

    # ---- Read header from first file to set up output ---- #
    with JdsFile(jds_files[0], nofs=args.nofs) as jds0:
        first_header = jds0.header

    # Construct output filename (matches IDL convention)
    cleaned_filename = outdir / f"Cleaned_ {args.label}{shortnames[0]}.ucd"

    # Write .ucd header (copy of .jds header with nofs patched)
    write_ucd_header(cleaned_filename, first_header, args.nofs)

    # Open .ucd for appending data after the header
    ucd_fh = open(cleaned_filename, "ab")

    try:
        for n_file in range(n_in_list):
            jds_path = jds_files[n_file]
            print(f"\n--- File {n_file + 1} of {n_in_list}: {jds_path.name} ---")

            with JdsFile(jds_path, nofs=args.nofs) as jds:
                hdr = jds.header
                print(f"  Name:  {hdr.sname}")
                print(f"  Local: {hdr.stime}")
                print(f"  UTC:   {hdr.sgmtt}")
                print(f"  Mode:  {'waveform' if hdr.data_mode == 0 else 'spectra' if hdr.data_mode == 1 else 'correlation'}")
                print(f"  Fmin: {hdr.fmin_mhz:.1f} MHz, Fmax: {hdr.fmax_mhz:.1f} MHz, "
                      f"wofsg: {hdr.wofsg}, avrs: {hdr.avrs}, "
                      f"TimeRes: {hdr.time_res_s * 1000:.3f} ms")
                print(f"  Frames: {jds.nframe}")

                tt = time.time()

                for i, imdat in jds.frames(mode=args.mode):
                    # Clean the frame
                    imdat, mask = adr_cleaning(imdat, pat=True, nar=True, wid=True)

                    # Write cleaned data as float32 to .ucd
                    ucd_fh.write(imdat.astype(np.float32).tobytes(order="F"))

                    elapsed = time.time() - tt
                    print(f"  Frame {i + 1}/{jds.nframe}, Time = {elapsed:.2f}s")
                    tt = time.time()

    finally:
        ucd_fh.close()

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
