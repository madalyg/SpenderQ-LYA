"""Pipeline inputs: the quasar CSV, per-quasar folder names, and SDSS-like FITS spectra."""

import csv
from pathlib import Path

import numpy as np
from astropy.io import fits

from spenderq.desi_qso import DESI

OBSERVATION_COLUMNS = ("Observation 1", "Observation 2")


def read_quasar_csv(csv_path):
    """Return one dict per CSV row with keys name, observations, redshift.

    observations lists the non-empty Observation 1 / Observation 2 filenames, so a row
    with only one epoch gives a one-element list.
    """
    rows = []
    with open(csv_path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            observations = [(row.get(col) or "").strip() for col in OBSERVATION_COLUMNS]
            rows.append({
                "name": row["Quasar Name"].strip(),
                "observations": [obs for obs in observations if obs],
                "redshift": float(row["Redshift"].strip()),
            })
    return rows


def quasar_dir(quasar_name):
    """Return the directory stem: everything before the '+' in the quasar name."""
    idx = quasar_name.find("+")
    return quasar_name[:idx] if idx != -1 else quasar_name


def get_mjd(path):
    """Get the MJD from a filename like spec-PLATE-MJD-FIBER.fits."""
    parts = Path(path).stem.split("-")
    # Plates above 9999 are also five digits, so take the MJD by position for SDSS names.
    if parts[0] == "spec" and len(parts) >= 4:
        return parts[2]
    for part in parts:
        if len(part) == 5 and part.isdigit():
            return part
    return Path(path).stem


def get_column(data, *names):
    columns = {name.lower(): name for name in data.names}
    for name in names:
        if name.lower() in columns:
            return np.asarray(data[columns[name.lower()]], dtype=float)
    return None


def load_sdss_like_spectrum(path, z_qso):
    """Load an SDSS spectrum and return observed wavelength, flux, ivar."""
    with fits.open(path) as hdul:
        hdu = None
        for candidate in ("COADD", 1):
            try:
                data = hdul[candidate].data
            except (KeyError, IndexError):
                continue
            if getattr(data, "names", None):
                hdu = data
                break

        if hdu is None:
            raise ValueError(f"No table HDU with spectral columns found in {path}")

        loglam = get_column(hdu, "loglam")
        wave = 10.0**loglam if loglam is not None else get_column(
            hdu, "wave", "wavelength", "lambda"
        )
        flux = get_column(hdu, "flux")
        ivar = get_column(hdu, "ivar", "inverse_variance")
        err = get_column(hdu, "err", "error", "sigma")

        if wave is None or flux is None:
            raise ValueError(f"Missing wavelength or flux column in {path}")
        if ivar is None:
            if err is None:
                ivar = np.ones_like(flux)
            else:
                ivar = np.zeros_like(err)
                good_err = np.isfinite(err) & (err > 0)
                ivar[good_err] = 1.0 / err[good_err] ** 2

    # If wavelengths look rest-frame (too small for observed DESI grid), convert to observed A
    if np.nanpercentile(wave, 95) < np.nanmin(np.asarray(DESI._wave_obs)):
        wave = wave * (1.0 + z_qso)

    return wave, flux, ivar
