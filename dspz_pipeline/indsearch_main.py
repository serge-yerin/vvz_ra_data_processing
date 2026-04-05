#!/usr/bin/env python3
"""
IndSearch — Incoherent pulsar dedispersion search for UTR-2 .ucd files.

Translated from IDL (IndSearch.pro) to Python.

Usage
-----
    python main.py <file.ucd> <DM>  [--no-plot]
    python main.py "data/Cleaned_ PSRB0834p06A141010_032001.jds.ucd" "12.8579"

Arguments
---------
    file.ucd   Path to the input observation file
    DM         Central dispersion measure [pc/cm^3]

Options
-------
    --no-plot  Write the .dmt output file but skip the matplotlib figure
"""

import argparse

import numpy as np

from dspz_pipeline.io.ucd import read_ucd_header_as_dict, count_frames
from dspz_pipeline.analysis.indsearch import process
from dspz_pipeline.io.dmt import write_dmt


def main():
    parser = argparse.ArgumentParser(
        description="Pulsar DM search \u2014 IndSearch (UTR-2 .ucd files)",
    )
    parser.add_argument("filename", help="Input data file (.ucd)")
    parser.add_argument("dm", type=float, help="Central DM [pc/cm^3]")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip the output plot")
    args = parser.parse_args()

    # --- Read header and print observation info ---
    params = read_ucd_header_as_dict(args.filename)
    params["n_kadr"] = count_frames(args.filename, params["wofsg"], params["nofs"])

    # Replicate IDL print statements
    print(
        f"\n   File name:  {params['sname']} \n"
        f"   Local time: {params['stime']}  \n"
        f"   UTC time:   {params['sgmtt']} \n"
        f"   Record:     {params['ssysn']} \n"
        f"   Place:      {params['splace']} \n"
        f"   Descriptor: {params['sdesc']} \n"
    )
    # IDL: print, HStr.Sdspp[6:32]  — IDL slices are inclusive, so elements 6..32
    print(params["sdspp"][6:33])

    print(
        f"  Fmin =           {params['fmin']:.4f} MHz \n"
        f"  Fmax =           {params['fmax']:.4f} MHz \n"
        f"  Channels =       {params['wofsg']} \n"
        f"  Spectra/frame =  {params['nofs']} \n"
        f"  Averages =       {params['avrs']} \n"
        f"  TimeRes =        {params['time_res']:.6e} s \n"
        f"  Frames =         {params['n_kadr']} \n"
    )

    # --- Run the DM search ---
    with open(args.filename, "rb") as f:
        acc_dm, max_shift, pic_size, dm_step_numb = process(f, params, args.dm)

    # --- Write binary output (.dmt) ---
    output_path = args.filename + ".dmt"
    out_slice = acc_dm[:, max_shift:max_shift + pic_size]   # shape (51, 65536)

    write_dmt(output_path, out_slice, dm_step_numb + 1, pic_size)
    print(f"\nOutput written to: {output_path}")

    # --- Optional plot ---
    if not args.no_plot:
        from dspz_pipeline.gui.dm_time_plot import plot_dm_time
        plot_dm_time(out_slice, args.dm, dm_step_numb)


if __name__ == "__main__":
    main()
