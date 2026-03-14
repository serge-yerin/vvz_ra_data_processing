#!/usr/bin/env python3
"""
IndSearch — Incoherent pulsar dedispersion search for UTR-2 .ucd files.

Translated from IDL (IndSearch.pro) to Python.

Usage
-----
    python main.py <file.ucd> <DM>  [--no-plot]

Arguments
---------
    file.ucd   Path to the input observation file
    DM         Central dispersion measure [pc/cm^3]

Options
-------
    --no-plot  Write the .dmt output file but skip the matplotlib figure
"""

import argparse
import sys
import numpy as np

from indsearch import read_header, count_frames, process


def main():
    parser = argparse.ArgumentParser(
        description="Pulsar DM search — IndSearch (UTR-2 .ucd files)",
    )
    parser.add_argument("filename", help="Input data file (.ucd)")
    parser.add_argument("dm",       type=float, help="Central DM [pc/cm^3]")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip the output plot")
    args = parser.parse_args()

    # --- Read header and print observation info ---
    with open(args.filename, "rb") as f:
        params = read_header(f)

    params["n_kadr"] = count_frames(args.filename, params["wofsg"], params["nofs"])

    # Replicate IDL print statements
    print(
        f"File name: {params['sname']}"
        f"   Local time: {params['stime']}"
        f"   UTC time: {params['sgmtt']}"
        f"   Record: {params['ssysn']}"
        f"   Place: {params['splace']}"
        f"   Descriptor: {params['sdesc']}"
    )
    # IDL: print, HStr.Sdspp[6:32]  — IDL slices are inclusive, so elements 6..32
    print(params["sdspp"][6:33])

    print(
        f"Fmin={params['fmin']:.4f} MHz"
        f"  Fmax={params['fmax']:.4f} MHz"
        f"  Channels={params['wofsg']}"
        f"  Spectra/frame={params['nofs']}"
        f"  Averages={params['avrs']}"
        f"  TimeRes={params['time_res']:.6e} s"
        f"  Frames={params['n_kadr']}"
    )

    # --- Run the DM search ---
    with open(args.filename, "rb") as f:
        acc_dm, max_shift, pic_size, dm_step_numb = process(f, params, args.dm)

    # --- Write binary output (.dmt), matching IDL writeu byte layout ---
    output_path = args.filename + ".dmt"
    out_slice = acc_dm[:, max_shift:max_shift + pic_size]   # shape (51, 65536)

    with open(output_path, "wb") as fout:
        # IDL: writeu,5, Dmstepnumb+1L, picsize
        # Each is a single 32-bit signed integer (IDL LONG), little-endian
        np.array([dm_step_numb + 1, pic_size], dtype="<i4").tofile(fout)

        # IDL: writeu,5, accDM[*, maxshift:picsize+maxshift-1]
        # IDL writeu writes 2D arrays in Fortran (column-major) order:
        # element [dm, t] at byte offset dm + t*(dm_step_numb+1) within the block.
        # Transposing and then writing C-order reproduces that byte layout.
        out_slice.T.astype("<f4").tofile(fout)

    print(f"\nOutput written to: {output_path}")

    # --- Optional plot ---
    if not args.no_plot:
        from indsearch.plot import plot_dm_time
        plot_dm_time(out_slice, args.dm, dm_step_numb)


if __name__ == "__main__":
    main()
