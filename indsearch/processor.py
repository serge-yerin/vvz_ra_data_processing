"""
Main incoherent dedispersion search loop.
Replicates the body of IDL procedure IndSearch.
"""

import time as _time
import numpy as np

from .dispersion import compute_dm_delays
from .reader     import read_frame

# Constants that are hard-coded in the IDL source
RANGE        = 256      # frequency channels per subband  (IDL: range = 256)
MAX_SHIFT    = 65536    # time-axis padding length         (IDL: maxshift = 65536L)
DM_STEP      = 0.004    # DM step [pc/cm^3]               (IDL: DMstep = 0.004)
DM_STEP_NUMB = 50       # number of DM steps              (IDL: Dmstepnumb = 50)
PIC_SIZE     = 65536    # output time samples             (IDL: picsize = 65536L)


def process(f, params, dm_const):
    """
    Perform incoherent dedispersion over a grid of DM values.

    Parameters
    ----------
    f        : binary file object (open for reading)
    params   : dict returned by read_header(), augmented with 'n_kadr'
    dm_const : float  central DM value [pc/cm^3]

    Returns
    -------
    acc_dm       : ndarray, shape (DM_STEP_NUMB+1, MAX_SHIFT + n_kadr*nofs)
    max_shift    : int   (= MAX_SHIFT)
    pic_size     : int   (= PIC_SIZE)
    dm_step_numb : int   (= DM_STEP_NUMB)
    """
    nofs     = params["nofs"]
    fmin     = params["fmin"]
    fmax     = params["fmax"]
    wofsg    = params["wofsg"]
    time_res = params["time_res"]
    n_kadr   = params["n_kadr"]
    df       = (fmax - fmin) / wofsg

    n_subbands = wofsg // RANGE
    total_time = MAX_SHIFT + n_kadr * nofs

    # Pre-allocate arrays — shapes mirror IDL declarations
    # IDL: temp   = fltarr(range, maxshift + i_kadr*nofs)
    # IDL: accDM  = fltarr(Dmstepnumb+1, maxshift + i_kadr*nofs)
    # IDL: accSB  = fltarr(wofsg/range, maxshift + i_kadr*nofs)
    temp   = np.zeros((RANGE, total_time), dtype=np.float32)
    acc_dm = np.zeros((DM_STEP_NUMB + 1, total_time), dtype=np.float32)
    acc_sb = np.zeros((n_subbands, total_time), dtype=np.float32)

    # Outer loop: DM steps
    # IDL: for l = 0, Dmstepnumb do begin
    for l in range(DM_STEP_NUMB + 1):

        # IDL: DM = DM_Const - Dmstepnumb/2 * DMstep + l * DMstep
        # Dmstepnumb/2 is integer division in IDL (50/2 = 25)
        dm = dm_const - (DM_STEP_NUMB // 2) * DM_STEP + l * DM_STEP
        dt = compute_dm_delays(dm, fmax, fmin, wofsg, df)

        t0 = _time.time()

        # Inner loop: subbands
        # IDL: fqbeg=0  fqend=wofsg/range
        # IDL: for k = fqbeg, fqend-1 do begin
        for k in range(n_subbands):
            rbeg = k * RANGE

            # IDL: temp[*,*] = 0.
            temp[:] = 0.0

            # Load all frames into the working buffer
            # IDL: for j = 0, i_kadr-1 do begin
            #        data = datdspz[j]
            #        temp[0:range-1, maxshift+j*nofs : maxshift+(j+1)*nofs-1]
            #              = data[rbeg:rbeg+range-1, *]
            #      endfor
            for j in range(n_kadr):
                data = read_frame(f, j, wofsg, nofs)
                # data shape: (wofsg, nofs)  — data[i_freq, i_time]
                # temp shape: (RANGE, total_time)
                # IDL inclusive end rbeg+range-1 → Python exclusive end rbeg+RANGE
                # IDL inclusive end maxshift+(j+1)*nofs-1 → Python maxshift+(j+1)*nofs
                col_start = MAX_SHIFT + j * nofs
                col_end   = MAX_SHIFT + (j + 1) * nofs
                temp[:, col_start:col_end] = data[rbeg:rbeg + RANGE, :]

            # Dedispersion: circular-shift each channel by its DM delay
            # IDL: for i = 0, range-1 do begin
            #        temp[i, *] = shift(temp[i, *], -round(dt[i+rbeg] / TimeRes))
            #      endfor
            # IDL shift(arr, n) ≡ np.roll(arr, n) — same sign convention
            shifts = -np.round(dt[rbeg:rbeg + RANGE] / time_res).astype(int)
            for i in range(RANGE):
                if shifts[i] != 0:
                    temp[i, :] = np.roll(temp[i, :], int(shifts[i]))

            # Integrate over frequency channels within the subband
            # IDL: accSB[k, *] = total(temp, 1)
            # total(arr, 1) sums along IDL dim 1 = numpy axis 0
            acc_sb[k, :] = temp.sum(axis=0)

            elapsed = _time.time() - t0
            print(
                f"DM step {l:3d}  from {DM_STEP_NUMB}"
                f"  subband {k:3d}"
                f"  processing time {elapsed:.2f}s"
            )
            t0 = _time.time()

        # Integrate over subbands
        # IDL: accDM[l, *] = total(accSB, 1)
        acc_dm[l, :] = acc_sb.sum(axis=0)

    return acc_dm, MAX_SHIFT, PIC_SIZE, DM_STEP_NUMB
