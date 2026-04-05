"""
Subband-based incoherent dedispersion search loop (Stage 2).

Replicates the body of IDL procedure IndSearch.pro.

This is the stage 2 algorithm that operates on .ucd files with a
subband decomposition approach (8 subbands of 256 channels each).
"""

import time as _time

import numpy as np

from dspz_pipeline.config import (
    DM_HALF_STEPS,
    DM_STEP_SIZE,
    DM_TOTAL_STEPS,
)
from dspz_pipeline.io.dmt import compute_dm_delays
from dspz_pipeline.io.ucd import read_frame

# Constants that are hard-coded in the IDL source
RANGE = 256          # frequency channels per subband  (IDL: range = 256)
MAX_SHIFT = 65536    # time-axis padding length         (IDL: maxshift = 65536L)
PIC_SIZE = 65536     # output time samples             (IDL: picsize = 65536L)


def process(f, params, dm_const):
    """
    Perform subband-based incoherent dedispersion over a grid of DM values.

    Parameters
    ----------
    f        : binary file object (open for reading)
    params   : dict returned by read_ucd_header_as_dict(), augmented with 'n_kadr'
    dm_const : float  central DM value [pc/cm^3]

    Returns
    -------
    acc_dm       : ndarray, shape (DM_TOTAL_STEPS, MAX_SHIFT + n_kadr*nofs)
    max_shift    : int   (= MAX_SHIFT)
    pic_size     : int   (= PIC_SIZE)
    dm_step_numb : int   (= DM_HALF_STEPS * 2)
    """
    nofs = params["nofs"]
    fmin = params["fmin"]
    fmax = params["fmax"]
    wofsg = params["wofsg"]
    time_res = params["time_res"]
    n_kadr = params["n_kadr"]
    df = (fmax - fmin) / wofsg

    n_subbands = wofsg // RANGE
    total_time = MAX_SHIFT + n_kadr * nofs
    dm_step_numb = DM_HALF_STEPS * 2  # 50

    # Pre-allocate arrays — shapes mirror IDL declarations
    temp = np.zeros((RANGE, total_time), dtype=np.float32)
    acc_dm = np.zeros((dm_step_numb + 1, total_time), dtype=np.float32)
    acc_sb = np.zeros((n_subbands, total_time), dtype=np.float32)

    # Outer loop: DM steps
    for l in range(dm_step_numb + 1):

        dm = dm_const - (dm_step_numb // 2) * DM_STEP_SIZE + l * DM_STEP_SIZE
        dt = compute_dm_delays(dm, fmax, fmin, wofsg, df)

        t0 = _time.time()

        # Inner loop: subbands
        for k in range(n_subbands):
            rbeg = k * RANGE

            temp[:] = 0.0

            # Load all frames into the working buffer
            for j in range(n_kadr):
                data = read_frame(f, j, wofsg, nofs)
                col_start = MAX_SHIFT + j * nofs
                col_end = MAX_SHIFT + (j + 1) * nofs
                temp[:, col_start:col_end] = data[rbeg:rbeg + RANGE, :]

            # Dedispersion: circular-shift each channel by its DM delay
            shifts = -np.round(dt[rbeg:rbeg + RANGE] / time_res).astype(int)
            for i in range(RANGE):
                if shifts[i] != 0:
                    temp[i, :] = np.roll(temp[i, :], int(shifts[i]))

            # Integrate over frequency channels within the subband
            acc_sb[k, :] = temp.sum(axis=0)

            elapsed = _time.time() - t0
            print(
                f"DM step {l:3d}  from {dm_step_numb}"
                f"  subband {k:3d}"
                f"  processing time {elapsed:.2f}s"
            )
            t0 = _time.time()

        # Integrate over subbands
        acc_dm[l, :] = acc_sb.sum(axis=0)

    return acc_dm, MAX_SHIFT, PIC_SIZE, dm_step_numb
