"""Two-panel plots comparing a pair of reconstructed epochs."""

import matplotlib.pyplot as plt
import numpy as np

from spenderq.plotting.figures import finish_figure, mark_lines
from spenderq.reconstruction.spectrum_io import get_mjd
from spenderq.spectral_lines import LYA, LYB, lya_region_observed

SDSS_WINDOW_OBS = (3800.0, 6000.0)


def plot_continuum_ratio_panels(a, b, ratio_ab, ratio_ba, out_dir, show=False):
    """Plot continuum_A / continuum_B and the reverse over the Lyβ–Lyα forest."""
    z_qso = a.z_qso
    lo, hi = lya_region_observed(z_qso)
    in_range = (a.wave >= lo) & (a.wave <= hi)
    mjd_a, mjd_b = get_mjd(a.path), get_mjd(b.path)

    fig, axes = plt.subplots(2, 1, figsize=(12, 6), constrained_layout=True)
    for ax, ratio, ylabel in zip(
        axes,
        [ratio_ab, ratio_ba],
        [f"{mjd_a} recon / {mjd_b} recon", f"{mjd_b} recon / {mjd_a} recon"],
    ):
        ax.plot(a.wave[in_range], ratio[in_range], color="black", lw=0.9)
        ax.axhline(1.0, color="0.5", ls="--", lw=0.8)
        ax.axvline(LYA * (1.0 + z_qso), color="tab:purple", ls="--", lw=0.9)
        ax.set_ylabel(ylabel, fontsize=12)
        ax.set_xlabel("observed wavelength [Å]")
        ax.set_xlim(lo, hi)
        r = ratio[in_range & np.isfinite(ratio)]
        ylo = (np.nanmin(r) - 0.05) if r.size else 0.85
        ax.set_ylim(ylo, 1.5)

    outpath = out_dir / f"{a.path.stem}_vs_{b.path.stem}_continuum_ratio_panels.png"
    fig.suptitle(
        f"SpenderQ continuum ratios: {a.path.stem} vs {b.path.stem} "
        f"(Lyα forest: {LYB:.0f}–{LYA:.0f} Å rest)",
        fontsize=14,
    )
    return finish_figure(fig, outpath, show)


def _plot_epoch_pair(a, b, top, bottom, ylabel_op, baseline, ylim, kind, out_dir, show):
    """Shared layout for the (flux/recon) ratio and difference plots: A vs B on top, B vs A below."""
    z_qso = a.z_qso
    mjd_a, mjd_b = get_mjd(a.path), get_mjd(b.path)

    fig, axes = plt.subplots(2, 1, figsize=(12, 6), sharex=True, constrained_layout=True)
    for ax, values, num, den in [
        (axes[0], top, mjd_a, mjd_b),
        (axes[1], bottom, mjd_b, mjd_a),
    ]:
        ax.plot(a.wave, values, color="black", lw=0.8)
        ax.axhline(baseline, color="0.5", ls="--", lw=0.8)
        ax.set_ylabel(f"MJD {num} {ylabel_op} MJD {den}")
        ax.set_ylim(*ylim)
    mark_lines(axes[0], z_qso, labeled=("N V", "C IV"))
    mark_lines(axes[1], z_qso, labeled=())
    axes[1].set_xlabel("observed wavelength [Å]")
    axes[1].set_xlim(*SDSS_WINDOW_OBS)

    suffix = "cont_norm_ratio" if kind == "ratio" else "cont_norm_ratio_difference"
    outpath = out_dir / f"{a.path.stem}_vs_{b.path.stem}_{suffix}_observed.png"
    fig.suptitle(f"(flux/cont) {kind} in observed frame: {a.path.stem} vs {b.path.stem}", fontsize=14)
    return finish_figure(fig, outpath, show)


def plot_recon_norm_ratio(a, b, ratio_ab, ratio_ba, out_dir, show=False):
    """Plot (flux/recon)_A / (flux/recon)_B and the reverse in the observed frame."""
    return _plot_epoch_pair(a, b, ratio_ab, ratio_ba, "/", 1.0, (-0.2, 4.0), "ratio", out_dir, show)


def plot_recon_norm_difference(a, b, diff_ab, diff_ba, out_dir, show=False):
    """Plot (flux/recon)_A − (flux/recon)_B and the reverse in the observed frame."""
    return _plot_epoch_pair(a, b, diff_ab, diff_ba, "−", 0.0, (-1.5, 1.5), "difference", out_dir, show)
