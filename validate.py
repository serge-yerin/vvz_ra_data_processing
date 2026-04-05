"""
Validation script: compare Python pipeline output against IDL reference files.

Usage:
    python validate.py
"""

import os
import sys
import struct
import numpy as np


def compare_binary_files(our_path, ref_path, label, dtype=np.float32, header_skip=0):
    """Compare two binary files element by element."""
    print(f"\n{'='*60}")
    print(f"  {label}")
    print(f"{'='*60}")
    print(f"  Our file:  {our_path} ({os.path.getsize(our_path):,} bytes)")
    print(f"  Ref file:  {ref_path} ({os.path.getsize(ref_path):,} bytes)")

    our_size = os.path.getsize(our_path)
    ref_size = os.path.getsize(ref_path)

    if our_size != ref_size:
        print(f"  *** SIZE MISMATCH: {our_size} vs {ref_size} ***")
        return False

    # Compare header
    if header_skip > 0:
        with open(our_path, "rb") as f:
            our_hdr = f.read(header_skip)
        with open(ref_path, "rb") as f:
            ref_hdr = f.read(header_skip)
        if our_hdr == ref_hdr:
            print(f"  Header ({header_skip} bytes): IDENTICAL")
        else:
            n_diff = sum(a != b for a, b in zip(our_hdr, ref_hdr))
            print(f"  Header ({header_skip} bytes): {n_diff} bytes differ")

    # Compare data
    n_elements = (our_size - header_skip) // np.dtype(dtype).itemsize
    our_data = np.fromfile(our_path, dtype=dtype, offset=header_skip, count=n_elements)
    ref_data = np.fromfile(ref_path, dtype=dtype, offset=header_skip, count=n_elements)

    identical = np.sum(our_data == ref_data)
    total = our_data.size
    pct = 100 * identical / total

    diff = np.abs(our_data - ref_data)
    max_diff = diff.max()
    mean_diff = diff.mean()

    print(f"  Data elements: {total:,}")
    print(f"  Identical:     {identical:,} / {total:,} ({pct:.4f}%)")
    print(f"  Max abs diff:  {max_diff:.6e}")
    print(f"  Mean abs diff: {mean_diff:.6e}")
    print(f"  Differ > 1e-6: {np.sum(diff > 1e-6):,}")
    print(f"  Differ > 1e-4: {np.sum(diff > 1e-4):,}")
    print(f"  Differ > 1e-2: {np.sum(diff > 1e-2):,}")

    # Zero pattern match
    our_z = our_data == 0
    ref_z = ref_data == 0
    zero_match = np.sum(our_z == ref_z)
    print(f"  Zero pattern:  {zero_match:,} / {total:,} match")

    if pct == 100.0:
        print(f"  RESULT: PERFECT MATCH")
        return True
    elif pct > 99.9:
        print(f"  RESULT: EXCELLENT MATCH (>99.9%)")
        return True
    elif pct > 99.0:
        print(f"  RESULT: GOOD MATCH (>99%)")
        return True
    else:
        print(f"  RESULT: DIFFERENCES DETECTED")
        return False


def compare_dmt_files(our_path, ref_path):
    """Compare two .dmt files, handling different picsize (IDL hardcodes 65536)."""
    from dspz_pipeline.io.dmt import read_dmt

    print(f"\n{'='*60}")
    print(f"  DMT file comparison")
    print(f"{'='*60}")
    print(f"  Our file:  {our_path} ({os.path.getsize(our_path):,} bytes)")
    print(f"  Ref file:  {ref_path} ({os.path.getsize(ref_path):,} bytes)")

    our_data, our_dm, our_pic = read_dmt(our_path)
    ref_data, ref_dm, ref_pic = read_dmt(ref_path)

    print(f"  Our:  DMstepnumb={our_dm}, picsize={our_pic}")
    print(f"  Ref:  DMstepnumb={ref_dm}, picsize={ref_pic}")

    if our_dm != ref_dm:
        print(f"  *** DM step count mismatch ***")
        return False

    # Compare overlapping region
    min_pic = min(our_pic, ref_pic)
    our_trunc = our_data[:, :min_pic]
    ref_trunc = ref_data[:, :min_pic]

    if our_pic != ref_pic:
        print(f"  NOTE: picsize differs (IDL hardcodes 65536). "
              f"Comparing first {min_pic} samples.")

    diff = np.abs(our_trunc - ref_trunc)
    identical = np.sum(our_trunc == ref_trunc)
    total = our_trunc.size
    pct = 100 * identical / total

    print(f"  Data elements: {total:,}")
    print(f"  Identical:     {identical:,} / {total:,} ({pct:.4f}%)")
    print(f"  Max abs diff:  {diff.max():.6e}")
    print(f"  Mean abs diff: {diff.mean():.6e}")
    print(f"  Differ > 1e-2: {np.sum(diff > 1e-2):,}")
    print(f"  Differ > 1e-1: {np.sum(diff > 1e-1):,}")
    print(f"  Differ > 1:    {np.sum(diff > 1):,}")

    # Per DM step summary (sample 5 steps)
    print(f"  Per-DM-step sample:")
    for j in [0, ref_dm // 4, ref_dm // 2, 3 * ref_dm // 4, ref_dm - 1]:
        d = np.abs(our_trunc[j] - ref_trunc[j])
        step_pct = 100 * np.sum(our_trunc[j] == ref_trunc[j]) / min_pic
        print(f"    j={j:2d}: max_diff={d.max():.4e}, mean_diff={d.mean():.4e}, "
              f"identical={step_pct:.2f}%")

    if diff.max() < 2.0:
        print(f"  RESULT: GOOD MATCH (max diff < 2.0, from input .ucd rounding)")
        return True
    else:
        print(f"  RESULT: DIFFERENCES DETECTED")
        return False


def main():
    our_ucd = "_output/Cleaned_ PSRB0834p06A141010_032001.jds.ucd"
    ref_ucd = "sample_inter/Cleaned_ PSRB0834p06A141010_032001.jds.ucd"
    our_dmt = "_output/Cleaned_ PSRB0834p06A141010_032001.jds.ucd.dmt"
    ref_dmt = "sample_inter/Cleaned_ PSRB0834p06A141010_032001.jds.ucd.dmt"

    results = []

    if os.path.exists(our_ucd) and os.path.exists(ref_ucd):
        results.append(compare_binary_files(our_ucd, ref_ucd, "UCD file comparison", header_skip=1024))
    else:
        print(f"UCD files not found. Run the pipeline first.")
        if not os.path.exists(our_ucd):
            print(f"  Missing: {our_ucd}")

    if os.path.exists(our_dmt) and os.path.exists(ref_dmt):
        results.append(compare_dmt_files(our_dmt, ref_dmt))
    else:
        print(f"\nDMT files not found. Run dedispersion first.")
        if not os.path.exists(our_dmt):
            print(f"  Missing: {our_dmt}")

    if results:
        print(f"\n{'='*60}")
        if all(results):
            print("  ALL CHECKS PASSED")
        else:
            print("  SOME CHECKS FAILED")
        print(f"{'='*60}")


if __name__ == "__main__":
    main()
