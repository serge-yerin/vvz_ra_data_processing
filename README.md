# DSPZ Pipeline

A Python pipeline for processing low-frequency radio astronomy data from DSPZ receivers (0-40 MHz range) to search for pulsar and transient signals. It reads raw `.jds` binary data files, removes radio frequency interference (RFI), performs incoherent dedispersion to search for pulsar and transient signals, and provides interactive visualization tools for inspecting the results.

This is a direct translation of the original IDL codebase written by Dr. Vyacheslav Zakharenko for the DSPZ (Digital Spectro-Polarimeter "Z") instrument used at the UTR-2 radio telescope and URAN VLBI system.

## Data Processing Pipeline

The pipeline has two stages that run sequentially:

```
Stage 1                                    Stage 2
.jds files --> RFI cleaning --> .ucd --> dedispersion search --> .dmt --> plot
(raw data)                     (clean)  (IndSearch)             (DM-time plane)
```

**Stage 1** (`process_survey`) reads raw `.jds` binary data, cleans RFI, and
writes cleaned spectrogram data to `.ucd` files. It also includes its own
dedispersion step and interactive Tkinter GUIs for pulse inspection.

**Stage 2** (`indsearch_main`) takes the `.ucd` files from Stage 1 and performs a
subband-based incoherent dedispersion search over a grid of DM trial values,
producing a `.dmt` file and an optional DM-time plane plot.

## Project Structure

```
dspz_pipeline/                          # Python package
├── __init__.py
├── config.py                           # All shared constants and default parameters
├── utils.py                            # Shared utilities (smoothing, text decoding)
│
├── io/                                 # File I/O modules
│   ├── jds_reader.py                   # .jds file reading, decoding, JdsHeader dataclass
│   ├── ucd.py                          # .ucd header reading, frame reading (shared)
│   └── dmt.py                          # .dmt read/write, DM delay computation (shared)
│
├── cleaning/                           # RFI cleaning pipeline
│   ├── robust_stats.py                 # Iterative sigma-clipping (erov, background)
│   ├── flatten.py                      # Per-channel bandpass normalisation
│   ├── sumthr.py                       # SumThreshold RFI flagging algorithm
│   └── rfi_cleaning.py                 # Full cleaning driver (adr_cleaning)
│
├── analysis/                           # Dedispersion algorithms
│   ├── dedispersion.py                 # Stage 1: per-channel shift-and-sum
│   └── indsearch.py                    # Stage 2: subband-based search loop
│
├── gui/                                # Interactive visualization tools
│   ├── trans_search.py                 # DM-vs-time spectrogram viewer (Tkinter)
│   ├── show_pulse.py                   # Individual pulse inspector (Tkinter)
│   ├── repeating_analysis.py           # FFT periodicity analysis
│   └── dm_time_plot.py                 # DM-time plane matplotlib plot
│
├── process_survey.py                   # Stage 1 entry point: .jds --> .ucd --> .dmt
└── indsearch_main.py                   # Stage 2 entry point: .ucd --> .dmt

_data/                                  # Sample .jds input files
_output/                                # Generated output files (.ucd, .dmt, plots)
idl_code/                               # Original IDL source code (for reference)
validate.py                             # Utility to compare output against IDL reference
.vscode/launch.json                     # VSCode debug configurations
```

## Installation

### With uv (recommended)

```bash
uv venv
uv pip install -r requirements.txt
```

### With pip

```bash
python3.12 -m venv .venv
# On Linux/macOS:
source .venv/bin/activate
# On Windows:
.venv\Scripts\activate

pip install -r requirements.txt
```

### As an editable package (for development)

```bash
pip install -e .
```

This installs the command-line entry points (`dspz-process`, `dspz-indsearch`,
etc.) so you can run them from anywhere.

## How to Run

### Stage 1: Full pipeline (.jds to .ucd to .dmt)

```bash
python -m dspz_pipeline.process_survey --files _data/A141010_032001.jds _data/A141010_032843.jds --dm 12.872 --label "PSRB0834p06" --outdir _output
```

**What it does:**
- Reads each `.jds` file frame by frame
- Decodes the DSPZ binary format
- Cleans each frame (removes RFI)
- Writes the cleaned data to `_output/Cleaned_ PSRB0834p06A141010_032001.jds.ucd`
- Runs incoherent dedispersion over 51 DM trial values
- Writes dedispersed data to `_output/Cleaned_ PSRB0834p06A141010_032001.jds.ucd.dmt`
- Launches the interactive TransSearch GUI

**Expected runtime:** ~7 minutes per file (128 frames total for 2 files),
plus ~10 minutes for dedispersion. Progress is printed to the console.

Add `--no-gui` to skip the GUI and just produce the output files.

### Stage 2: IndSearch dedispersion (.ucd to .dmt)

```bash
python -m dspz_pipeline.indsearch_main "_output/Cleaned_ PSRB0834p06A141010_032001.jds.ucd"  12.872
```

**What it does:**
- Reads the `.ucd` file produced by Stage 1
- Performs subband-based incoherent dedispersion over 51 DM trial values
  (8 subbands of 256 channels each)
- Writes the DM-time plane to a `.dmt` file
- Displays a DM-time plane plot saved to `_output/dm_time_plane.png`
  (skip with `--no-plot`)

### Stage 1 interactive tools (on existing output files)

#### TransSearch GUI (DM-vs-time spectrogram viewer)

```bash
python -m dspz_pipeline.gui.trans_search "_output/Cleaned_ PSRB0834p06A141010_032001.jds.ucd.dmt" 12.872
```

**Controls:**
- **smpar < >** -- adjust high-pass smoothing (1-8)
- **smpar_b < >** -- adjust low-pass background removal (32-4096)
- **DM step < >** -- shift DM offset from central value
- **IND** -- toggle individual pulse analysis (click spectrogram to inspect)
- **REP** -- toggle repeating pulse FFT analysis
- **parts / N of parts** -- select time window for FFT analysis
- **Min/Max scale sliders** -- adjust display contrast
- **Click on spectrogram** -- select a time point for pulse inspection

#### Individual pulse viewer

```bash
python -m dspz_pipeline.gui.show_pulse "_output/Cleaned_ PSRB0834p06A141010_032001.jds.ucd" 12.872 --ns 32768
```

**Controls:**
- **DM +/-** -- fine-tune dispersion measure (step = 0.002 pc/cm^3)
- **Shift < >** -- shift the dedispersion in time
- **Narr / Wide** -- adjust frequency resolution for spectrum display

#### Repeating pulse FFT analysis

```bash
python -m dspz_pipeline.gui.repeating_analysis "_output/Cleaned_ PSRB0834p06A141010_032001.jds.ucd.dmt" 12.872
```

### Running the full pipeline end-to-end

To run both stages sequentially from raw data to final analysis:

```bash
# Stage 1: clean raw data and produce .ucd
python -m dspz_pipeline.process_survey --files _data/A141010_032001.jds _data/A141010_032843.jds --dm 12.872 --label "PSRB0834p06" --outdir _output --no-gui

# Stage 2: run IndSearch dedispersion on the .ucd output
python -m dspz_pipeline.indsearch_main "_output/Cleaned_ PSRB0834p06A141010_032001.jds.ucd" 12.872
```

## Command-Line Arguments Reference

### process_survey (Stage 1)

| Argument    | Default         | Description                                    |
|-------------|-----------------|------------------------------------------------|
| `--files`   | (required)      | One or more `.jds` files to process            |
| `--outdir`  | `_output`       | Output directory                               |
| `--dm`      | `12.872`        | Central dispersion measure (pc/cm^3)           |
| `--period`  | `1.29224132384` | Pulsar period (seconds)                        |
| `--label`   | `PSRB0834p06`   | Source label for output filenames              |
| `--mode`    | `1`             | Channel: 0 = ch0-ch1, 1 = ch0, 2 = ch1       |
| `--nofs`    | `1024`          | Spectra per frame                              |
| `--no-gui`  | (flag)          | Skip launching the interactive GUI             |

### indsearch_main (Stage 2)

| Argument    | Default    | Description                                    |
|-------------|------------|------------------------------------------------|
| `filename`  | (required) | Path to `.ucd` input file                      |
| `dm`        | (required) | Central dispersion measure (pc/cm^3)           |
| `--no-plot` | (flag)     | Skip the matplotlib output figure              |

## File Formats

### .jds (DSPZ raw data)

The `.jds` format is used by the DSPZ digital receiver. Each file is typically
~2 GB and represents a continuous segment of observation.

**Header (1024 bytes):**
- Bytes 0-31: File name (ASCII, null-padded)
- Bytes 32-63: Local time string
- Bytes 64-95: UTC time string
- Bytes 96-127: System name
- Bytes 128-159: System type
- Bytes 160-255: Observation place
- Bytes 256-511: Description
- Bytes 512-1023: Parameter array (`Sdspp`, 128 x UINT32)
  - `Sdspp[24]`: Data mode (0=waveform, 1=spectra, 2=correlation)
  - `Sdspp[28]`: Low frequency channel number
  - `Sdspp[29]`: High frequency channel number
  - `Sdspp[30]`: Number of frequency channels (`wofsg`)
  - `Sdspp[31]`: Averaging count (`avrs`)

**Data (after header):**
- Sequence of frames, each containing `UINT32[2 x wofsg x nofs]`
- Two channels (typically left/right circular polarisation)
- Custom floating-point encoding: bits 5-31 = mantissa, bits 0-4 = exponent
- Decode formula: `float(mantissa) / 2^exponent * 8192 / 4294967296 / avrs`

**Typical parameters:**
- 4096 frequency channels spanning 16.5-33.0 MHz
- 64x averaging -> ~7.94 ms time resolution
- 1024 spectra per frame x 64 frames per file = 65536 spectra per file

### .ucd (cleaned data)

Intermediate file produced by Stage 1 after RFI cleaning.

**Header (1024 bytes):** Identical structure to `.jds`, but `Sdspp[26]` is
patched with the actual frame size used during processing.

**Data:** float32 spectrograms in Fortran (column-major) order, one frame
of shape `(wofsg, nofs)` per write.

### .dmt (dedispersed DM-time plane)

Output file from dedispersion (both Stage 1 and Stage 2).

**Structure:**
- Bytes 0-3: int32 -- number of DM steps (typically 51)
- Bytes 4-7: int32 -- number of time samples (`picsize`)
- Remaining: float32 array of shape `(n_dm, picsize)` in Fortran (column-major) order

**Content:** 51 dedispersed time series, one per DM trial value. DM values
range from `dm_const - 25 * 0.004` to `dm_const + 25 * 0.004` pc/cm^3.

## Running and Debugging in VSCode

### Recommended Extensions

- **Python** (`ms-python.python`) -- Python language support
- **Pylance** (`ms-python.vscode-pylance`) -- type checking and IntelliSense
- **Debugpy** -- built-in Python debugger (included with the Python extension)

### Debug Configurations

The `.vscode/launch.json` file includes pre-configured debug configurations:

1. **Stage 1: Process Survey (full pipeline)** -- runs the complete pipeline with GUI
2. **Stage 1: Process Survey (no GUI)** -- runs the pipeline without the interactive GUI
3. **Stage 1: TransSearch GUI** -- opens the interactive analysis on existing `.dmt` output
4. **Stage 1: ShowPulse (standalone)** -- opens the pulse viewer on existing `.ucd` output
5. **Stage 1: Repeating Analysis (standalone)** -- runs FFT analysis on existing `.dmt` output
6. **Stage 2: IndSearch** -- runs subband dedispersion on existing `.ucd` output
7. **Stage 2: IndSearch (no plot)** -- same as above, without the matplotlib figure

To use: open the Run and Debug panel (Ctrl+Shift+D), select a configuration
from the dropdown, and press F5.

## Troubleshooting

### "ModuleNotFoundError: No module named 'dspz_pipeline'"

Make sure you are running from the project root directory, or install the
package with `pip install -e .`. If using `python -m`, the project root
must be in your Python path:

```bash
# On Linux/macOS:
export PYTHONPATH=.
# On Windows:
set PYTHONPATH=.
```

### "tkinter not found" or "No module named '_tkinter'"

tkinter is required for the interactive GUIs (Stage 1 only). It is included
with most Python installations but may need to be installed separately:

```bash
# Ubuntu/Debian:
sudo apt install python3-tk
# Fedora:
sudo dnf install python3-tkinter
# macOS (Homebrew):
brew install python-tk@3.12
# Windows: included with the official Python installer
```

### Processing is slow

**Stage 1:** The RFI cleaning step processes ~6 seconds per frame. For two
2 GB files (128 frames), expect ~13 minutes for cleaning plus ~10 minutes
for dedispersion.

**Stage 2:** IndSearch processes 51 DM steps x 8 subbands, reading all frames
for each combination. Runtime depends on file size and disk speed.

Both stages print progress to the console.

### Output files don't match the reference exactly

Small numerical differences (< 1e-6 relative) are expected due to
floating-point arithmetic differences between IDL and Python/numpy.
The scientific results should be identical.

### "MemoryError" when processing large files

The pipeline processes data frame by frame and should not require more
than ~200 MB of RAM. If you encounter memory issues, ensure you are
using a 64-bit Python installation.
