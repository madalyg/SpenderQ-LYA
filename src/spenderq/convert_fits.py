"""
convert_fits.py

Converts every dereddened text spectrum (.dr16, .dr9, .txt, ...) in a parent folder into
an observed-frame FITS file that SpenderQ can read, so init_directories.py can sort them
into per-quasar folders.

Each input file has three columns: wavelength, flux, error. The output FITS has a COADD
table with loglam, flux and ivar = 1 / error^2 (ivar = 0 where the error is missing or
not positive). Flux is left in its dereddened units; analyze_ehvo.py normalizes it.

Files are matched to the CSV by PLATE-MJD-FIBER in their names (e.g. 6138-57328-0746-dered.dr16,
spec-6122-56246-0696norm.dr16). Rest-frame files are shifted to the observed frame with the
redshift from the CSV; files that are already observed-frame are written as they are.

Usage:
    python convert_fits.py [CSV_PATH] [INPUT_DIR] [OUTPUT_DIR]

    CSV_PATH    quasar CSV (default: dat/Variability/NV/NV_variability_list_fits.csv)
    INPUT_DIR   folder of dereddened text spectra (default: <CSV folder>/input-files-unnormalized)
    OUTPUT_DIR  where the FITS files go (default: <CSV folder>/input-files-fits)

Example:
    python convert_fits.py dat/Dense_Lya/J2318_list_fits.csv
"""

import csv
import re
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits

INPUT_SUBDIR = "input-files-unnormalized"
OUTPUT_SUBDIR = "input-files-fits"
OBSERVATION_COLUMNS = ("Observation 1", "Observation 2")
PLATE_MJD_FIBER = re.compile(r"(\d{4,5})-(\d{5})-(\d{3,4})")
# SDSS spectra start near 3600 Å observed; anything bluer than this has to be rest-frame.
MIN_OBSERVED_WAVE = 3000.0


def spec_id(filename):
    """Return (plate, mjd, fiber) ints from any filename containing PLATE-MJD-FIBER."""
    match = PLATE_MJD_FIBER.search(Path(filename).name)
    if match is None:
        return None
    return tuple(int(part) for part in match.groups())


def fits_name(plate, mjd, fiber):
    return f"spec-{plate:04d}-{mjd:05d}-{fiber:04d}.fits"


def read_redshifts(csv_path):
    """Map (plate, mjd, fiber) -> redshift for every observation listed in the CSV."""
    redshifts = {}
    with open(csv_path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            z = float(row["Redshift"].strip())
            for col in OBSERVATION_COLUMNS:
                key = spec_id((row.get(col) or "").strip())
                if key is None:
                    continue
                if key in redshifts and redshifts[key] != z:
                    print(f"  [warn] {fits_name(*key)} listed with z={redshifts[key]} and z={z}; using {z}")
                redshifts[key] = z
    return redshifts


def text_to_fits(path, out_dir, redshifts):
    """Convert one text spectrum to an observed-frame FITS file. Returns the output path or None."""
    key = spec_id(path)
    if key is None:
        print(f"  [skip] {path.name}: no PLATE-MJD-FIBER in the filename")
        return None
    z_qso = redshifts.get(key)

    data = np.loadtxt(path)
    wave, flux, err = data[:, 0], data[:, 1], data[:, 2]

    if np.nanmin(wave) < MIN_OBSERVED_WAVE:
        if z_qso is None:
            print(f"  [skip] {path.name}: rest-frame wavelengths but {fits_name(*key)} is not in the CSV")
            return None
        wave = wave * (1.0 + z_qso)
        frame = f"rest frame, z={z_qso}"
    else:
        frame = "observed frame"

    ivar = np.zeros_like(err)
    good = np.isfinite(err) & (err > 0) & np.isfinite(flux)
    ivar[good] = 1.0 / err[good] ** 2

    primary = fits.PrimaryHDU()
    if z_qso is not None:
        primary.header["Z"] = z_qso
    primary.header["DERED"] = (True, "Galactic extinction corrected")
    primary.header["SRCFILE"] = path.name
    table = fits.BinTableHDU.from_columns([
        fits.Column(name="loglam", format="E", array=np.log10(wave)),
        fits.Column(name="flux", format="E", array=flux),
        fits.Column(name="ivar", format="E", array=ivar),
    ])
    table.name = "COADD"

    outpath = out_dir / fits_name(*key)
    fits.HDUList([primary, table]).writeto(outpath, overwrite=True)
    print(f"  {path.name} ({frame}) -> {outpath}")
    return outpath


def convert_folder(csv_path, input_dir, output_dir):
    csv_path, input_dir, output_dir = Path(csv_path), Path(input_dir), Path(output_dir)
    if not input_dir.is_dir():
        raise FileNotFoundError(f"Input folder not found: {input_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)

    redshifts = read_redshifts(csv_path)
    inputs = sorted(p for p in input_dir.iterdir() if p.is_file() and not p.name.startswith("."))
    print(f"Converting {len(inputs)} file(s) from {input_dir} -> {output_dir}\n")

    written = [out for out in (text_to_fits(p, output_dir, redshifts) for p in inputs) if out]
    missing = sorted(set(redshifts) - {spec_id(p) for p in written})
    for key in missing:
        print(f"  [missing] {fits_name(*key)} is in {csv_path.name} but has no file in {input_dir}")
    print(f"\nWrote {len(written)} FITS file(s) to {output_dir.resolve()}")
    return written


if __name__ == "__main__":
    default_csv = Path("dat/Variability/NV/NV_variability_list_fits.csv")
    csv_arg = Path(sys.argv[1]) if len(sys.argv) > 1 else default_csv
    input_arg = Path(sys.argv[2]) if len(sys.argv) > 2 else csv_arg.parent / INPUT_SUBDIR
    output_arg = Path(sys.argv[3]) if len(sys.argv) > 3 else csv_arg.parent / OUTPUT_SUBDIR

    convert_folder(csv_arg, input_arg, output_arg)
