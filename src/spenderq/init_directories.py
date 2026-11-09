5"""
init_directories.py

Reads a variability CSV (format: Quasar Name, Observation 1, Observation 2, Redshift),
creates one subdirectory per quasar, and fills it with that quasar's FITS spectra.

By default the spectra are the dereddened FITS files written by convert_fits.py
(<CSV folder>/input-files-fits/). They are copied into the quasar folders, replacing any
file of the same name (e.g. a raw SAS download); identical files are left alone so
finished analyses are not rerun. With --download, the raw SDSS DR16 spectra are
downloaded from SAS instead.

Usage:
    python init_directories.py [CSV_PATH] [BASE_DIR] [--fits-dir DIR | --download]

    CSV_PATH    path to the CSV file (default: dat/Variability/NV/NV_variability_list_fits.csv)
    BASE_DIR    root directory to create quasar subdirectories in
                (default: same directory as the CSV file)
    --fits-dir  folder of converted FITS files (default: <CSV folder>/input-files-fits)
    --download  download raw spectra from SDSS SAS instead of copying local files

Example:
    python init_directories.py dat/Variability/NV/NV_variability_list_fits.csv
    python init_directories.py dat/Dense_Lya/J2318_list_fits.csv --download
"""

import argparse
import csv
import filecmp
import shutil
import urllib.request
import urllib.error
from pathlib import Path

FITS_SUBDIR = "input-files-fits"

# SDSS DR16 SAS base URLs tried in order for each spectrum.
# Older SDSS-I/II plates (roughly < 3510) live under sdss/spectro/redux/26;
# BOSS/eBOSS plates are under eboss/spectro/redux/v5_13_0.
# Tries both so you do not need to know which survey a plate belongs to.

SDSS_BASE_URLS = [
    "https://data.sdss.org/sas/dr16/eboss/spectro/redux/v5_13_0/spectra/full",
    "https://data.sdss.org/sas/dr16/sdss/spectro/redux/26/spectra/full",
]


def parse_fits_filename(filename):
    """Return (plate, mjd, fiber) ints from 'spec-PLATE-MJD-FIBER.fits'."""
    stem = Path(filename).stem          # e.g. spec-0948-52428-0370
    parts = stem.split("-")             # ['spec', '0948', '52428', '0370']
    if len(parts) != 4 or parts[0] != "spec":
        raise ValueError(f"Unexpected FITS filename format: {filename!r}")
    plate = int(parts[1])
    mjd   = int(parts[2])
    fiber = int(parts[3])
    return plate, mjd, fiber


def sdss_download_urls(filename):
    """Yield candidate SAS URLs for a given spec filename."""
    plate, mjd, fiber = parse_fits_filename(filename)
    fname = f"spec-{plate:04d}-{mjd:05d}-{fiber:04d}.fits"
    plate_dir = f"{plate:04d}"
    for base in SDSS_BASE_URLS:
        yield f"{base}/{plate_dir}/{fname}"


def download_fits(filename, dest_dir):
    """Download *filename* from SDSS SAS into dest_dir. Returns the local path."""
    dest = dest_dir / Path(filename).name
    if dest.exists():
        print(f"  [skip] already exists: {dest}")
        return dest

    last_error = None
    for url in sdss_download_urls(filename):
        try:
            print(f"  downloading {url}")
            urllib.request.urlretrieve(url, dest)
            print(f"  -> saved to {dest}")
            return dest
        except urllib.error.HTTPError as exc:
            last_error = exc
            if exc.code == 404:
                print(f"  [404] not found at {url}, trying next …")
                continue
            raise
        except Exception as exc:
            last_error = exc
            print(f"  [error] {exc}, trying next …")
            continue

    raise RuntimeError(
        f"Could not download {filename!r} from any known SAS location. "
        f"Last error: {last_error}"
    )


def copy_local_fits(filename, src_dir, dest_dir):
    """Copy *filename* from src_dir into dest_dir, replacing a different file of the same name."""
    src = src_dir / Path(filename).name
    dest = dest_dir / Path(filename).name
    if not src.exists():
        raise RuntimeError(f"{src.name} not found in {src_dir} (run convert_fits.py first?)")
    if dest.exists() and filecmp.cmp(src, dest, shallow=False):
        print(f"  [skip] already up to date: {dest}")
        return dest
    action = "replaced" if dest.exists() else "copied"
    # copyfile gives dest a fresh mtime, so analyze_ehvo.py reruns any quasar whose spectrum changed.
    shutil.copyfile(src, dest)
    print(f"  {action} {src} -> {dest}")
    return dest


def quasar_dir_name(quasar_name):
    """Return the directory name: everything up to (not including) the '+' sign."""
    idx = quasar_name.find("+")
    if idx == -1:
        return quasar_name          # no '+' (e.g. negative-dec quasar — keep full name)
    return quasar_name[:idx]


def init_directories(csv_path, base_dir, fits_dir=None):
    """Create the quasar folders and fill them from fits_dir, or from SAS when fits_dir is None."""
    csv_path = Path(csv_path)
    base_dir = Path(base_dir)

    if not csv_path.exists():
        raise FileNotFoundError(f"CSV not found: {csv_path}")
    if fits_dir is not None and not Path(fits_dir).is_dir():
        raise FileNotFoundError(f"FITS folder not found: {fits_dir} (run convert_fits.py first?)")

    base_dir.mkdir(parents=True, exist_ok=True)

    with open(csv_path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        rows = list(reader)

    source = f"local FITS in {fits_dir}" if fits_dir is not None else "SDSS SAS downloads"
    print(f"Found {len(rows)} quasar(s) in {csv_path}; spectra from {source}\n")

    for row in rows:
        name   = row["Quasar Name"].strip()
        obs1   = (row.get("Observation 1") or "").strip()
        obs2   = (row.get("Observation 2") or "").strip()
        z      = row["Redshift"].strip()

        dir_name = quasar_dir_name(name)
        qso_dir  = base_dir / dir_name
        qso_dir.mkdir(parents=True, exist_ok=True)

        print(f"{'='*60}")
        print(f"Quasar : {name}  (z = {z})")
        print(f"Dir    : {qso_dir}")

        for obs in (obs1, obs2):
            if not obs:
                continue
            try:
                if fits_dir is not None:
                    copy_local_fits(obs, Path(fits_dir), qso_dir)
                else:
                    download_fits(obs, qso_dir)
            except RuntimeError as exc:
                print(f"  [FAILED] {exc}")

    print(f"\nQuasar directories are now under: {base_dir.resolve()}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create per-quasar folders and fill them with FITS spectra.")
    parser.add_argument("csv_path", nargs="?", type=Path,
                        default=Path("dat/Variability/NV/NV_variability_list_fits.csv"))
    parser.add_argument("base_dir", nargs="?", type=Path,
                        help="root for the quasar folders (default: the CSV's folder)")
    source = parser.add_mutually_exclusive_group()
    source.add_argument("--fits-dir", type=Path,
                        help=f"converted FITS folder (default: <CSV folder>/{FITS_SUBDIR})")
    source.add_argument("--download", action="store_true",
                        help="download raw spectra from SDSS SAS instead of copying local files")
    args = parser.parse_args()

    base_dir = args.base_dir or args.csv_path.parent
    fits_dir = None if args.download else (args.fits_dir or args.csv_path.parent / FITS_SUBDIR)
    init_directories(args.csv_path, base_dir, fits_dir)
