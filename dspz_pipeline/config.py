"""
Shared constants and default parameters for the DSPZ pipeline.

These values are extracted from the original IDL code where they appeared
as hardcoded literals or COMMON-block variables.
"""

# --------------------------------------------------------------------------- #
#  Instrument constants (DSPZ receiver)
# --------------------------------------------------------------------------- #

# ADC clock frequency in Hz — used to convert averaging count to time resolution
CLOCK_FREQ_HZ = 66_000_000

# Number of ADC samples per FFT window
FFT_LENGTH = 8192

# Nyquist bandwidth of the receiver in MHz
NYQUIST_BANDWIDTH_MHZ = 33.0

# Number of FFT channels that span the Nyquist bandwidth
TOTAL_FFT_CHANNELS = 8192

# --------------------------------------------------------------------------- #
#  DSPZ binary header layout
# --------------------------------------------------------------------------- #

# Total header size in bytes (fixed for all DSPZ file formats)
HEADER_SIZE_BYTES = 1024

# Offset into the Sdspp array where instrument parameters begin.
# The first 16 UINT32 slots are reserved / used for other purposes.
HEAD_OFFSET = 16

# Byte widths of the text fields in the header (in order)
HEADER_TEXT_FIELDS = {
    "sname": 32,
    "stime": 32,
    "sgmtt": 32,
    "ssysn": 32,
    "ssyst": 32,
    "splace": 96,
    "sdesc": 256,
}

# The Sdspp parameter array: 128 × UINT32 (512 bytes), follows the text fields
SDSPP_COUNT = 128

# --------------------------------------------------------------------------- #
#  DSPZ data encoding bit-masks
# --------------------------------------------------------------------------- #

# Custom floating-point encoding used in .jds spectral data:
#   bits 5–31 = mantissa,  bits 0–4 = exponent
MANTISSA_MASK = 0xFFFF_FFE0   # 4294967232
EXPONENT_MASK = 0x0000_001F   # 31

# Scaling constant applied after decoding mantissa/exponent
# Derivation: 4 * 2 * 1024 / 2^32  =  8192 / 4294967296
DECODE_SCALE = 8192.0 / 4_294_967_296.0

# --------------------------------------------------------------------------- #
#  Default pipeline parameters (from process_survey.pro)
# --------------------------------------------------------------------------- #

# Number of spectra per processing frame
DEFAULT_NOFS = 1024

# Processing mode: 0 = multiplication (ch0*ch1), 1 = channel 1, 2 = channel 2
DEFAULT_MODE = 1

# Filename shift for date parsing (ShiftAB = 1 for DSPZ AB format)
SHIFT_AB = 1

# --------------------------------------------------------------------------- #
#  Dedispersion parameters
# --------------------------------------------------------------------------- #

# Number of DM trial steps on each side of the central DM
DM_HALF_STEPS = 25

# DM step size in pc cm^-3 (IDL IndSearch.pro: DMstep=0.004)
DM_STEP_SIZE = 0.004

# Total number of DM trial values (2 * 25 + 1 = 51)
DM_TOTAL_STEPS = 2 * DM_HALF_STEPS + 1

# Dispersion constant: delay(s) = DM / 2.41e-4 * (1/f_lo^2 - 1/f_hi^2)
# where frequencies are in MHz.  The IDL code uses 2.4103 with f in units
# that give 1e4/f^2, which is equivalent.
DISPERSION_CONSTANT = 2.41e-4  # s MHz^2 pc^-1 cm^3

# --------------------------------------------------------------------------- #
#  TransSearch GUI defaults
# --------------------------------------------------------------------------- #

# Signal-to-noise threshold for event detection
DEFAULT_SN_LIMIT = 5.5

# Default high-pass smoothing parameter (number of samples)
DEFAULT_SMPAR = 4

# Default low-pass smoothing exponent (background = 2^smpar_pow)
DEFAULT_SMPAR_POW = 9

# Frame size for TransSearch display (spectra per panel)
TRANS_KADR = 4096

# Number of spectra on each side of a selected point for pulse window
DEFAULT_NSFRAME = 50

# --------------------------------------------------------------------------- #
#  Default pulsar parameters (PSR B0834+06 — the sample data target)
# --------------------------------------------------------------------------- #

DEFAULT_DM_CONST = 12.872        # pc cm^-3
DEFAULT_PERIOD = 1.29224132384   # seconds
DEFAULT_PULSAR_LABEL = "PSRB0834p06"
