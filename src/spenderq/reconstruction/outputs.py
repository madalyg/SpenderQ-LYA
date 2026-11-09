"""Text and CSV outputs written by both the individual-epoch and multi-epoch analyses."""

import csv

import numpy as np

# Column layout of the config CSV that DR16Q/ABSORPTION/absorption.py reads.
ABSORPTION_CSV_HEADER = [
    "NORM SPECTRA FILE NAME", "REDSHIFT", "CALCULATED SNR", "NEEDS RECALCULATION", "Masked Regions",
]


def output_wavelength(wave_obs, z_qso, rest_frame=False):
    """Observed wavelengths, or rest-frame wavelengths (wave_obs / (1 + z)) when rest_frame is True."""
    return wave_obs / (1.0 + z_qso) if rest_frame else wave_obs


def save_recon_norm_txt(epoch, outpath, rest_frame=False):
    """Save flux / SpenderQ continuum of one epoch as a text file.

    Columns: wavelength (observed, or rest frame if rest_frame), flux_over_recon, err.
    """
    valid = np.isfinite(epoch.wave) & np.isfinite(epoch.recon_norm)
    wave = output_wavelength(epoch.wave, epoch.z_qso, rest_frame)
    np.savetxt(
        outpath,
        np.column_stack([wave[valid], epoch.recon_norm[valid], epoch.recon_norm_err[valid]]),
        fmt="%.6f",
        header=f"{'wave_rest' if rest_frame else 'wave_obs'}  flux_over_recon  err",
    )
    return outpath


def save_absorption_csv(out_csv, entries):
    """Write the absorption.py input CSV for entries of (txt_path, z_qso)."""
    with open(out_csv, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(ABSORPTION_CSV_HEADER)
        for txt_path, z in entries:
            writer.writerow([txt_path.name, z, 0, "N", "[]"])
    print(f"\nSaved summary CSV: {out_csv}")
    return out_csv
