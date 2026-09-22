"""End-to-end pipeline on tiny synthetic .jds files: the parallel path must
produce exactly the same .ucd / .dmt bytes as the sequential one."""

import hashlib

import numpy as np
import pytest

from dspz_pipeline.config import HEAD_OFFSET, NYQUIST_BANDWIDTH_MHZ, TOTAL_FFT_CHANNELS
from dspz_pipeline.io.jds_reader import JdsHeader, write_ucd_header
from dspz_pipeline.process_survey import parse_args, run_pipeline

WOFSG, NOFS, AVRS = 64, 1024, 64


def _write_jds(path, n_frames, seed):
    hdr = JdsHeader(sname="synthetic")
    hdr.sdspp[8 + HEAD_OFFSET] = 1                                          # spectra mode
    hdr.sdspp[12 + HEAD_OFFSET] = int(16.5 * TOTAL_FFT_CHANNELS / NYQUIST_BANDWIDTH_MHZ)
    hdr.sdspp[13 + HEAD_OFFSET] = int(33.0 * TOTAL_FFT_CHANNELS / NYQUIST_BANDWIDTH_MHZ)
    hdr.sdspp[14 + HEAD_OFFSET] = WOFSG
    hdr.sdspp[15 + HEAD_OFFSET] = AVRS
    write_ucd_header(path, hdr, NOFS)          # .jds and .ucd share the header format
    rng = np.random.default_rng(seed)
    # mantissa in the high bits, small exponent in the low 5 bits
    mant = rng.integers(2**20, 2**26, size=(n_frames, NOFS, WOFSG, 2), dtype=np.uint32) << 5
    expn = rng.integers(0, 4, size=mant.shape, dtype=np.uint32)
    with open(path, "ab") as fh:
        fh.write((mant | expn).tobytes())


def _sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _run(tmp_path, tag, workers, save_png=False):
    indir = tmp_path / "in"
    indir.mkdir(exist_ok=True)
    for k, name in enumerate(["a.jds", "b.jds"]):
        if not (indir / name).exists():
            _write_jds(indir / name, n_frames=2, seed=k)
    outdir = tmp_path / f"out_{tag}"
    argv = ["--indir", str(indir), "--files", "a.jds", "b.jds", "--dm", "12.872",
            "--label", "T", "--outdir", str(outdir), "--no-gui", "--workers", str(workers)]
    if save_png:
        argv.append("--save_cleaning_mask")
    run_pipeline(parse_args(argv))
    return outdir / "Cleaned_T_a.jds.ucd", outdir / "Cleaned_T_a.jds.ucd.dmt"


def test_parallel_output_identical_to_sequential(tmp_path):
    ucd1, dmt1 = _run(tmp_path, "seq", workers=1)
    ucd2, dmt2 = _run(tmp_path, "par", workers=2, save_png=True)
    assert ucd1.stat().st_size == 1024 + 4 * NOFS * WOFSG * 4      # 4 frames
    assert _sha(ucd1) == _sha(ucd2)
    assert _sha(dmt1) == _sha(dmt2)
    pngs = sorted((ucd2.parent / ucd2.stem).glob("*.png"))
    assert [p.name for p in pngs] == [f"{ucd2.stem}_{i}.png" for i in range(1, 5)]


def test_workers_default_is_positive():
    args = parse_args(["--files", "x.jds"])
    assert args.workers >= 1


def test_file_headers_printed_once_through_the_progress_bar(tmp_path, capsys):
    _run(tmp_path, "hdr", workers=1)
    out = capsys.readouterr().out
    assert out.count("--- File 1 of 2: a.jds ---") == 1
    assert out.count("--- File 2 of 2: b.jds ---") == 1
    assert out.count("Total number of frames in file:") == 2
