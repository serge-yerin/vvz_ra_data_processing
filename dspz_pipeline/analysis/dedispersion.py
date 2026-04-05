"""
Incoherent dedispersion for Stage 1 of the pipeline.

Reads a .ucd file, applies dedispersion over a range of DM trial values,
and writes the result to a .dmt file.

Translated from the dedispersion portion of IDL ``process_survey.pro``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from dspz_pipeline.config import (
    DM_HALF_STEPS,
    DM_STEP_SIZE,
    DM_TOTAL_STEPS,
    HEADER_SIZE_BYTES,
)
from dspz_pipeline.io.dmt import compute_dm_delays, write_dmt
from dspz_pipeline.io.jds_reader import read_ucd_header


def ind_search(
    ucd_path: str | Path,
    dm_const: float,
    fmax_mhz: float = 33.0,
    fmin_mhz: float = 16.5,
) -> Path:
    """Incoherent dedispersion search over a range of DM trial values.

    Translated from IDL ``IndSearch.pro``.

    The routine reads the .ucd file (1024-byte header + float32 spectrogram),
    applies dedispersion for each of 51 trial DM values
    (``dm_const - 25*0.004`` to ``dm_const + 25*0.004``), producing a
    dedispersed time series for each DM by shifting and summing across
    frequency channels.

    Parameters
    ----------
    ucd_path : str or Path
        Path to the cleaned ``.ucd`` file.
    dm_const : float
        Central dispersion measure in pc cm^-3.
    fmax_mhz : float
        Upper frequency edge in MHz.
    fmin_mhz : float
        Lower frequency edge in MHz.

    Returns
    -------
    dmt_path : Path
        Path to the output ``.dmt`` file.
    """
    ucd_path = Path(ucd_path)
    dmt_path = Path(str(ucd_path) + ".dmt")

    # Read .ucd header to get wofsg
    hdr = read_ucd_header(ucd_path)
    wofsg = hdr.wofsg
    time_res = hdr.time_res_s
    df_mhz = (fmax_mhz - fmin_mhz) / wofsg

    # Determine total number of spectra in the file
    file_size = ucd_path.stat().st_size
    picsize = (file_size - HEADER_SIZE_BYTES) // (4 * wofsg)

    print(f"IndSearch: wofsg={wofsg}, picsize={picsize}, "
          f"DM_const={dm_const}, DM range="
          f"[{dm_const - DM_HALF_STEPS * DM_STEP_SIZE:.3f}, "
          f"{dm_const + DM_HALF_STEPS * DM_STEP_SIZE:.3f}] \n")

    # Memory-map the .ucd data (skip 1024-byte header)
    ucd_data = np.memmap(
        ucd_path, dtype=np.float32, mode="r",
        offset=HEADER_SIZE_BYTES,
        shape=(picsize, wofsg),
    )

    # Allocate output: (DMstepnumb, picsize)
    dm_stepnumb = DM_TOTAL_STEPS
    acc_dm = np.zeros((dm_stepnumb, picsize), dtype=np.float32)

    for j in range(dm_stepnumb):
        dm = dm_const + (j - DM_HALF_STEPS) * DM_STEP_SIZE

        # Compute delays for this DM
        dt = compute_dm_delays(dm, fmax_mhz, fmin_mhz, wofsg, df_mhz)
        shifts = np.round(dt / time_res).astype(np.int64)
        max_shift = int(np.max(shifts))

        # Dedisperse: shift each frequency channel and sum
        dedispersed = np.zeros(picsize, dtype=np.float64)
        for ch in range(wofsg):
            s = int(shifts[ch])
            if s >= picsize:
                continue
            dedispersed[:picsize - s] += ucd_data[s:, ch].astype(np.float64)

        acc_dm[j, :] = dedispersed.astype(np.float32)

        if (j + 1) % 10 == 0 or j == 0 or j == dm_stepnumb - 1:
            print(f"  DM step {j + 1}/{dm_stepnumb}: DM={dm:.3f} pc/cm^3, "
                  f"max_shift={max_shift} samples")

    # Write .dmt file
    write_dmt(dmt_path, acc_dm, dm_stepnumb, picsize)

    print(f"\nIndSearch: wrote {dmt_path} "
          f"({dm_stepnumb} DM steps x {picsize} samples)")

    return dmt_path
