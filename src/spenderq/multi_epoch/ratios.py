"""Ratios between two reconstructed epochs of the same quasar, and their text files."""

import numpy as np

from spenderq.multi_epoch.success_criteria import safe_ratio
from spenderq.reconstruction.spectrum_io import get_mjd
from spenderq.spectral_lines import lya_region_observed


def continuum_ratios(a, b):
    """SpenderQ continuum_A / continuum_B and the reverse, on epoch A's grid.

    Dividing mutually independent reconstructions checks their stability between epochs.
    """
    cont_b = np.interp(a.wave, b.wave, b.continuum, left=np.nan, right=np.nan)
    return safe_ratio(a.continuum, cont_b), safe_ratio(cont_b, a.continuum)


def recon_norm_ratios(a, b):
    """(flux/recon)_A / (flux/recon)_B and the reverse; the Lyα forest common to both cancels."""
    return safe_ratio(a.recon_norm, b.recon_norm), safe_ratio(b.recon_norm, a.recon_norm)


def recon_norm_differences(a, b):
    """(flux/recon)_A − (flux/recon)_B and the reverse, sometimes easier to read for broad troughs."""
    return a.recon_norm - b.recon_norm, b.recon_norm - a.recon_norm


def save_continuum_ratio_txt(a, b, ratio_ab, ratio_ba, out_dir):
    """Save continuum ratios over the Lyβ–Lyα forest as <stem_A>_vs_<stem_B>_ratio_norm.txt
    and the reverse (read by summarize_continuum_ratios.py)."""
    lo, hi = lya_region_observed(a.z_qso)
    in_range = (a.wave >= lo) & (a.wave <= hi)
    saved = []
    for num, den, ratio in [(a, b, ratio_ab), (b, a, ratio_ba)]:
        valid = in_range & np.isfinite(ratio)
        outpath = out_dir / f"{num.path.stem}_vs_{den.path.stem}_ratio_norm.txt"
        np.savetxt(
            outpath,
            np.column_stack([a.wave[valid], ratio[valid], np.full(valid.sum(), np.nan)]),
            fmt="%.6f",
        )
        saved.append(outpath)
    return saved


def recon_norm_ratio_paths(out_dir, path_a, path_b):
    """<MJD_A>_over_<MJD_B>_recon.txt and the reverse inside out_dir."""
    mjd_a, mjd_b = get_mjd(path_a), get_mjd(path_b)
    return [
        out_dir / f"{mjd_a}_over_{mjd_b}_recon.txt",
        out_dir / f"{mjd_b}_over_{mjd_a}_recon.txt",
    ]


def save_recon_norm_ratio_txt(a, b, ratio_ab, ratio_ba, out_dir):
    """Save (flux/recon) ratios in the observed frame for absorption.py.

    Columns: observed wavelength, ratio, err (NaN).
    """
    paths = recon_norm_ratio_paths(out_dir, a.path, b.path)
    for outpath, ratio in zip(paths, [ratio_ab, ratio_ba]):
        valid = np.isfinite(a.wave) & np.isfinite(ratio)
        np.savetxt(
            outpath,
            np.column_stack([a.wave[valid], ratio[valid], np.full(valid.sum(), np.nan)]),
            fmt="%.6f",
        )
    return paths
