# Pipeline for manual search of pursar and transient emission

A Python code for processing low-frequency radio astronomy data from DSPZ receivers (0-40 MHz range) to search for pulsar and transient signals. It reads `.jds` binary data files recorded in spectra mode, removes low-level radio frequency interference (RFI), performs incoherent dedispersion to search for pulsar and transient signals, and provides interactive visualization tools for inspecting the results.

This is a direct translation of the original IDL codebase written by Dr. Vyacheslav Zakharenko for the DSPZ (Digital Spectro-Polarimeter "Z") instrument used at the UTR-2 radio telescope and URAN VLBI system. Pay attention, many variables are hardcoded for particular parameters observations (like 16.5 - 33.0 MHz range)

You can find more on the receivers in:


[1]	V. Zakharenko et al., “Digital Receivers for Low-Frequency Radio Telescopes UTR-2, URAN, GURT,” J. Astron. Instrum., vol. 05, no. 04, p. 1641010, Dec. 2016, doi: 10.1142/S2251171716410105.


[2]	V. B. Ryabov et al., “A low-noise, high-dynamic-range, digital receiver for radio astronomy applications: an efficient solution for observing radio-bursts from Jupiter, the Sun, pulsars, and other astrophysical plasmas below 30 MHz,” Astron. Astrophys., vol. 510, p. A16, Feb. 2010, doi: 10.1051/0004-6361/200913335.


## Data Processing Pipeline

The code consists of Full main pipeline (`process_survey`) and a part that uses an intermeiate result of the full pipline (.ucd file):

```
(raw data)                            (clean & normalized)                      (DM-time plane)
Full main pipeline                            
`.jds` --> Normalizing-->  RFI cleaning --> `.ucd` --> Multiple DM Dedispersion --> `.dmt` --> Plot & GUI analysis
Individual search branch                          |
                                                   --> Individual Search
```

**Full pipeline** (`process_survey`) reads raw `.jds` binary data (typically 2 consequtive files of 2 GB each recorded in spectra mode), cleans RFI, and
writes cleaned spectrogram data to `.ucd` files. It also includes its own
dedispersion step and interactive Tkinter GUIs for pulse inspection.

**Individual Search** (`indsearch_main`) takes the `.ucd` files from Full pipeline and performs a
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
│   ├── dedispersion.py                 # Full pipeline: per-channel shift-and-sum
│   └── indsearch.py                    # Individual Search: subband-based search loop
│
├── gui/                                # Interactive visualization tools
│   ├── trans_search.py                 # DM-vs-time spectrogram viewer (Tkinter)
│   ├── show_pulse.py                   # Individual pulse inspector (Tkinter)
│   ├── lod_image.py                    # Level-of-detail image rendering for large spectrograms
│   ├── repeating_analysis.py           # FFT periodicity analysis
│   └── dm_time_plot.py                 # DM-time plane matplotlib plot
│
├── process_survey.py                   # Full pipeline entry point: .jds --> .ucd --> .dmt
└── indsearch_main.py                   # Individual Search entry point: .ucd --> .dmt

_data/                                  # Sample .jds input files
_output/                                # Generated output files (.ucd, .dmt, plots)
idl_code/                               # Original IDL source code (for reference)
tests/                                  # pytest suite (pooling helpers + GUI smoke tests)
validate.py                             # Utility to compare output against IDL reference
.vscode/launch.json                     # VSCode debug configurations
```

## Installation

#### With uv (recommended)

```bash
uv venv
uv pip install -r requirements.txt
```

#### With pip

```bash
python3.12 -m venv .venv
# On Linux/macOS:
source .venv/bin/activate
# On Windows:
.venv\Scripts\activate

pip install -r requirements.txt
```

#### As an editable package (for development)

```bash
pip install -e .
```

This installs the command-line entry points (`dspz-process`, `dspz-indsearch`, etc.) so you can run them from anywhere.

## How to Run

### Full main pipeline 'Process survey' (.jds to .ucd to .dmt)

```bash
python -m dspz_pipeline.process_survey --indir _data_0834 --files A141010_032001.jds A141010_032843.jds --dm 12.872 --label "B0834p06" --outdir _output_0834 --save_cleaning_mask


python -m dspz_pipeline.process_survey --indir _data_1133 --files C231121_032738.jds C231121_033619.jds --dm 4.8471 --label "B1133p16" --outdir _output_1133 --save_cleaning_mask
```

**What it does:**
- Reads each `.jds` file frame by frame
- Decodes the DSPZ binary format
- Cleans (removes RFI) and normalizes response in each frame 
- Saves images of masks and cleaned data if `--save_cleaning_mask` flag applied 
- Writes the cleaned data to `_output_XXXX/Cleaned_<label>_<first .jds file>.ucd`
- Runs incoherent dedispersion over 51 DM trial values
- Writes dedispersed data to `_output_XXXX/Cleaned_<label>_<first .jds file>.ucd.dmt`
- Launches the interactive Transient Search GUI to analyze dedispersed data array

**Expected runtime:** about 3-4 minutes for two 2 GB files (128 frames) on an
8-core machine with the default `--workers`, including dedispersion
(see *Parallel cleaning* below). Progress is printed to the console.

Add `--no-gui` to skip the GUI and just produce the output files.

Add `--save_cleaning_mask` to save a PNG image of the cleaned data and RFI mask for every frame. Images are stored in a subfolder next to the `.ucd` file, named after the `.ucd` file (without extension). Each PNG is named `<ucd_stem>_NNN.png` (zero-padded by total frame count) and shows the Data and Mask arrays side-by-side with a Greys colormap. By default this is off to avoid the extra I/O overhead during processing.

#### Parallel cleaning: the `--workers` option

RFI cleaning is by far the most expensive step (~6 s of CPU per 1024-spectra
frame, i.e. ~13 min of CPU for two 2 GB files). Every frame is cleaned
independently of the others, so the pipeline cleans several frames at the
same time in separate **worker processes**:

- The main process reads the `.jds` files frame by frame and hands each raw
  frame to a worker. A worker decodes the frame, cleans it and (with
  `--save_cleaning_mask`) renders its PNG.
- The main process collects the results and appends them to the `.ucd` file
  **strictly in frame order**, keeping at most `2 x workers` frames in flight.
- Each frame goes through exactly the same code as in a single-process run,
  so **the `.ucd` and `.dmt` files are byte-identical for any value of
  `--workers`**. The option only changes how long the run takes.
- The following dedispersion step runs in the main process only; it is
  fast (~30 s) and not affected by `--workers`.

**Default:** number of CPU cores minus one (one core stays free for the
reading/writing main process and the OS), capped at 8. The value used is
printed at the start of the run (`RFI cleaning with N worker process(es)`).

**How to choose the value for your machine:**

| Situation | Recommended `--workers` |
|---|---|
| Typical desktop / laptop, 8+ cores, 16+ GB RAM | leave the default |
| 4 cores | `3` (the default gives this automatically) |
| 2 cores or a busy shared machine | `1` or `2` |
| Little RAM (8 GB or less) | `2`-`4`; watch memory use |
| Debugging, profiling, or a `MemoryError` / very slow run | `1` (sequential, no worker processes at all) |

**Memory:** each worker needs roughly **250 MB** while cleaning a frame (raw
frame + decoded float64 copy + intermediate arrays), or about **450 MB** with
`--save_cleaning_mask` (PNG rendering). In addition, dedispersion loads the
whole `.ucd` file into RAM (**~2 GB** for two 2 GB input files),
independently of `--workers`. So a run with 8 workers and PNGs needs roughly
8 x 0.45 + 2 = ~5.6 GB of free RAM at peak; if the machine starts swapping
(disk light on, run much slower than expected) or Python raises
`MemoryError`, lower `--workers`.

**CPU:** more workers than physical cores does not help; the cleaning is pure
CPU work. On a laptop, many workers may also trigger thermal throttling; if
the run is slower than expected, try one or two fewer.

**Speed-up is not perfectly linear:** the reader/writer and the operating
system need some CPU too, and PNG rendering adds ~1.7 s per frame of work to
the workers. Expect roughly 6 s / `workers` per frame overall.

`--workers 1` runs the cleaning inside the main process without creating any
worker processes. Use it when in doubt: it is the simplest configuration,
needs the least memory, and produces the same files (just slower).




### Transient Search GUI (DM-vs-time spectrogram viewer)

The last step of full pipeline that can be runned separately on calculated .dmt file  to start visual interactive analysis.
To run use command:

```bash
python -m dspz_pipeline.gui.trans_search "_output_0834/Cleaned_B0834p06_A141010_032001.jds.ucd.dmt" 12.872

python -m dspz_pipeline.gui.trans_search "_output_1133/Cleaned_B1133p16_C231121_032738.jds.ucd.dmt" 4.8471
```

#### Graphic User Interface description

In the end of the main pipeline or by running **dspz_pipeline.gui.dm_time_plot** on the existing .dmt file the main window of the user interface opens. It shows the 16 line plot of the data with compensated dispersion delay with DM value in the preset ranage (usually __ values). The script is intended to analyze data from 2 consecutively recoded 2GB .jds files (typically from UTR-2 pulsar-transient sky survey).
The main window has 2 modes of operation:
- Single pulse mode (button __ is ON)
- Repetitive pulse mode (button __ is ON)




**Controls** (grouped by what they affect):
- **Close** / **Save PNG** -- window-level actions
- **Smoothing** group -- shapes the displayed data and everything derived from it:
  - **smpar < >** -- adjust high-pass smoothing (1-8)
  - **smpar_b < >** -- adjust low-pass background removal (32-4096)
- **On plot click** group -- what a click on the spectrogram analyses. The pulse-selection
  window (2-D image + 3-D surface) always opens; the two checkboxes are independent and
  can both be on:
  - **DM step < >** -- DM offset from the central value used by both analyses below
  - **[ ] Individual pulse viewer** -- also open the individual pulse viewer and the Cleaned data window
  - **[ ] Repetitive (FFT) analysis** -- also run the repeating-pulse FFT analysis on the selected time window
  - **parts / N of parts** -- time window for the FFT analysis (shown as yellow highlighting of the spectrogram strips)
- **Min/Max scale sliders** -- adjust display contrast
- **Click on spectrogram** -- select a time point for pulse inspection


### Single pulse mode
When the button __ is ON user can 

### Repetitive pulse mode


### Individual pulse viewer 

GUI window opened in **Individual** analysis mode after click on spectrogram.
Shows ...

**Controls:**
- **DM +/-** -- fine-tune dispersion measure (step = 0.002 pc/cm^3)
- **Shift < >** -- shift the dedispersion in time
- **Narr / Wide** -- adjust frequency resolution for spectrum display


### Repeating pulse FFT analysis  
Plot window opened in **Repetitive** analysis mode after click on spectrogram.
Shows the ...






## Running the full pipeline end-to-end

To run both stages sequentially from raw data to final analysis:

```bash
# Full pipeline: clean raw data and produce .ucd
python -m dspz_pipeline.process_survey --indir _data --files A141010_032001.jds A141010_032843.jds --dm 12.872 --label "PSRB0834p06" --outdir _output --no-gui

# Individual Search run IndSearch dedispersion on the .ucd output
python -m dspz_pipeline.indsearch_main "_output/Cleaned_PSRB0834p06_A141010_032001.jds.ucd" 12.872
```




### Individual Search dedispersion (.ucd to .dmt)

```bash
python -m dspz_pipeline.indsearch_main "_output_0834/Cleaned_B0834p06_A141010_032001.jds.ucd"  12.872

python -m dspz_pipeline.indsearch_main "_output_1133/Cleaned_B1133p16_C231121_032738.jds.ucd"  4.8471
```

**What it does:**
It makes a dedisperion of the data in a range of DM points and shows the full data on a plot where dynamic spectrum and integrated over frequency time profile of individual pulses is shown, step by step:

- Reads the `.ucd` file produced by Full pipeline
- Performs subband-based incoherent dedispersion over 51 DM trial values  (8 subbands of 256 channels each)
- Writes the DM-time plane to a `.dmt` file
- Saves a DM-time plane PNG next to the `.dmt` file with the same stem (skip with `--no-plot`)

### DM-time plane viewer (standalone, on existing .dmt files)

If you already have a `.dmt` file and want to (re-)plot it without rerunning dedispersion:

```bash
python -m dspz_pipeline.gui.dm_time_plot "_output_0834/Cleaned_B0834p06_A141010_032001.jds.ucd.dmt" 12.872

python -m dspz_pipeline.gui.dm_time_plot "_output_1133/Cleaned_B1133p16_C231121_032738.jds.ucd.dmt" 4.8471
```

**What it does:**
- Reads the `.dmt` file from disk
- Displays the DM-time plane image and peak-DM time series
- Saves the PNG next to the `.dmt` file with the same name

### Full pipeline interactive tools (on existing output files)

#### TransientSearch GUI (DM-vs-time spectrogram viewer)

If you already calculated .dmt by the Full pipeline and want to start visual interactive analysis, you can run command:

```bash
python -m dspz_pipeline.gui.trans_search "_output_0834/Cleaned_PSRB0834p06_A141010_032001.jds.ucd.dmt" 12.872

python -m dspz_pipeline.gui.trans_search "_output_1133/Cleaned_B1133p16_C231121_032738.jds.ucd.dmt" 4.8471
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

**Cleaned data window** (opened together with the pulse viewer) shows the raw cleaned spectrogram
from the `.ucd` file around the selected pulse (up to 44000 samples x 4096 channels).

Layout: the matplotlib toolbar (home / back / forward / pan / zoom / save) is in the top row;
the left column holds the **Overview pooling** selector and the **vmax** / **vmin** contrast
sliders; the plot fills the rest of the window and re-flows when the window is resized or
maximised.

Rendering: only the currently visible part of the array is drawn, pooled down to screen
resolution (see `gui/lod_image.py`), so the sliders and toolbar zoom/pan respond in a fraction of a
second instead of ~10 s, and the current zoom/pan is kept when a slider is moved. When several
samples fall onto one screen pixel they are combined according to **Overview pooling**:
`max` (default; narrow pulses and RFI stay visible at any zoom), `mean` (smoother overview) or
`nearest` (plain sub-sampling, i.e. what plain `imshow` does). Zoom in far enough and the
individual samples are shown as they are. The slider values are thresholds in units of the data
STD around its mean and are applied *after* pooling. **Save PNGs** re-renders the visible window
at the PNG's resolution, not the screen's.

#### Individual pulse viewer (we do not run it separately)

```bash
python -m dspz_pipeline.gui.show_pulse "_output/Cleaned_PSRB0834p06_A141010_032001.jds.ucd" 12.872 --ns 32768
```

**Controls:**
- **DM +/-** -- fine-tune dispersion measure (step = 0.002 pc/cm^3)
- **Shift < >** -- shift the dedispersion in time
- **Narr / Wide** -- adjust frequency resolution for spectrum display

#### Repeating pulse FFT analysis  (we do not run it separately)

```bash
python -m dspz_pipeline.gui.repeating_analysis "_output/Cleaned_PSRB0834p06_A141010_032001.jds.ucd.dmt" 12.872
```

### Running the full pipeline end-to-end

To run both stages sequentially from raw data to final analysis:

```bash
# Full pipeline: clean raw data and produce .ucd
python -m dspz_pipeline.process_survey --indir _data --files A141010_032001.jds A141010_032843.jds --dm 12.872 --label "PSRB0834p06" --outdir _output --no-gui

# Individual Search run IndSearch dedispersion on the .ucd output
python -m dspz_pipeline.indsearch_main "_output/Cleaned_PSRB0834p06_A141010_032001.jds.ucd" 12.872
```

## Command-Line Arguments Reference

### process_survey (Full pipeline)

| Argument    | Default         | Description                                    |
|-------------|-----------------|------------------------------------------------|
| `--indir`   | `.`             | Directory containing the `.jds` input files    |
| `--files`   | (required)      | One or more `.jds` file names (relative to `--indir`) |
| `--outdir`  | `_output`       | Output directory                               |
| `--dm`      | `12.872`        | Central dispersion measure (pc/cm^3)           |
| `--period`  | `1.29224132384` | Pulsar period (seconds)                        |
| `--label`   | `PSRB0834p06`   | Source label for output filenames              |
| `--mode`    | `1`             | Channel: 0 = ch0-ch1, 1 = ch0, 2 = ch1         |
| `--nofs`    | `1024`          | Spectra per frame                              |
| `--no-gui`  | (flag)          | Skip launching the interactive GUI             |
| `--save_cleaning_mask` | (flag) | Save PNG images of cleaned data and RFI mask for each frame (see below) |
| `--workers` | cores - 1 (max 8) | Worker processes for RFI cleaning; `1` = sequential |

### indsearch_main (Individual Search)

| Argument    | Default    | Description                                    |
|-------------|------------|------------------------------------------------|
| `filename`  | (required) | Path to `.ucd` input file                      |
| `dm`        | (required) | Central dispersion measure (pc/cm^3)           |
| `--no-plot` | (flag)     | Skip the matplotlib output figure              |

### dm_time_plot (standalone DM-time viewer)

| Argument      | Default  | Description                                    |
|---------------|----------|------------------------------------------------|
| `dmt_file`    | (required) | Path to the `.dmt` file to plot              |
| `dm`          | (required) | Central dispersion measure (pc/cm^3)         |
| `--dm-step`   | `0.004`  | DM step size [pc/cm^3]                         |



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

Intermediate file produced by Full pipeline after RFI cleaning.

**Header (1024 bytes):** Identical structure to `.jds`, but `Sdspp[26]` is
patched with the actual frame size used during processing.

**Data:** float32 spectrograms in Fortran (column-major) order, one frame
of shape `(wofsg, nofs)` per write.

### .dmt (dedispersed DM-time plane)

Output file from dedispersion (both Full pipeline and Individual Search).

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

1. **Full pipeline: Process Survey (full pipeline)** -- runs the complete pipeline with GUI
2. **Full pipeline: Process Survey (no GUI)** -- runs the pipeline without the interactive GUI
3. **Full pipeline: TransSearch GUI** -- opens the interactive analysis on existing `.dmt` output
4. **Full pipeline: ShowPulse (standalone)** -- opens the pulse viewer on existing `.ucd` output
5. **Full pipeline: Repeating Analysis (standalone)** -- runs FFT analysis on existing `.dmt` output
6. **Individual Search: IndSearch** -- runs subband dedispersion on existing `.ucd` output
7. **Individual Search: IndSearch (no plot)** -- same as above, without the matplotlib figure

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

tkinter is required for the interactive GUIs (Full pipeline only). It is included
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
On macOS install .venv from the main python path:

```bash
# Find brew python path
brew --prefix python

# Use the full path to brew's python
/opt/homebrew/opt/python@3.12/bin/python3.12 -m venv .venv
# adjust version number to what you have installed

# Activate venv and install dependencies
source .venv/bin/activate
pip install -r requirements.txt

# Reselect interpreter in VSCode
# Ctrl+Shift+P → Python: Select Interpreter → pick the .venv one

# In your terminal:
brew install python-tk@3.13

# Verify it works
python -c "import tkinter; tkinter._test()"
```



### Running the tests

```bash
python -m pytest tests
```

The GUI smoke tests open real Tk windows for a few seconds; they are skipped automatically
if no display is available.

### Processing is slow

**Full pipeline:** RFI cleaning takes ~6 s of CPU per frame, but frames are
cleaned in parallel (`--workers`), so two 2 GB files (128 frames) take ~3 min
on 8 workers plus ~30 s for dedispersion. If the run is much slower than that,
check the number of workers printed at the start (`RFI cleaning with N worker
process(es)`), make sure the machine is not swapping, and see *Parallel
cleaning: the `--workers` option* above for how to pick a value that fits
your CPU and RAM.

**Run fails with `MemoryError` or the machine becomes unresponsive:** too many
workers for the available RAM. Re-run with a smaller `--workers` (`1` always
works), see the memory notes in the `--workers` section.

**Individual Search:** IndSearch processes 51 DM steps (by default) x 8 subbands, reading all frames for each combination. Runtime depends on file size and disk speed.

Both stages print progress to the console.

### Output files don't match the reference exactly

Small numerical differences (< 1e-6 relative) are expected due to floating-point arithmetic differences between IDL and Python/numpy.
The scientific results should be identical.

### "MemoryError" when processing large files

The pipeline processes data frame by frame and should not require more than ~200 MB of RAM. If you encounter memory issues, ensure you are using a 64-bit Python installation.
