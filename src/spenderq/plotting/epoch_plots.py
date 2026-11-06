"""Per-epoch plots of a SpenderQ reconstruction, made in both analysis modes."""

import matplotlib.pyplot as plt
import numpy as np

from spenderq.plotting.figures import finish_figure, mark_lines
from spenderq.spectral_lines import LYA, LYB, RECON_PLOT_XLIM_REST


def masked_forest_pixels(epoch):
    """Lyα-forest pixels SpenderQ weighted down when fitting the continuum."""
    wave_rest = epoch.wave / (1.0 + epoch.z_qso)
    lya_forest = (wave_rest > LYB) & (wave_rest < LYA)
    return lya_forest & (epoch.ivar > 0) & (epoch.weight_clipped == 0)


def plot_reconstruction(epoch, out_dir, observed=False, show=False):
    """Plot the spectrum and its SpenderQ continuum in the rest or observed frame."""
    scale = (1.0 + epoch.z_qso) if observed else 1.0
    x = epoch.wave if observed else epoch.wave / (1.0 + epoch.z_qso)
    masked = masked_forest_pixels(epoch)

    fig, ax = plt.subplots(figsize=(12, 4), constrained_layout=True)
    ax.plot(x, epoch.flux, color="0.25", lw=0.7, label="spectrum")
    ax.plot(x, epoch.continuum, color="tab:blue", lw=1.2, label="SpenderQ continuum")
    ax.scatter(
        x[masked],
        epoch.flux[masked],
        s=5,
        color="tab:orange",
        alpha=0.5,
        label="masked Ly-alpha forest",
    )

    wave_lo, wave_hi = (lim * scale for lim in RECON_PLOT_XLIM_REST)
    in_win = (x > wave_lo) & (x < wave_hi)
    ylim_data = np.concatenate([epoch.flux[in_win], epoch.continuum[in_win]])
    ylim_data = ylim_data[np.isfinite(ylim_data)]
    ylo, yhi = np.nanpercentile(ylim_data, [0.5, 99.5])
    pad = 0.1 * (yhi - ylo)
    ax.set_xlim(wave_lo, wave_hi)
    ax.set_ylim(ylo - pad, yhi + pad)
    ax.set_ylabel("normalized flux")
    ax.set_xlabel("observed wavelength [Å]" if observed else "rest-frame wavelength [Å]")
    ax.legend(loc="upper right", fontsize=14)
    mark_lines(ax, epoch.z_qso if observed else 0.0)

    name = epoch.path.name
    if observed:
        outpath = out_dir / f"{epoch.path.stem}_spenderq_ehvo_observed.png"
        fig.suptitle(f"{name}: SpenderQ continuum (observed frame)", fontsize=14)
    else:
        outpath = out_dir / f"{epoch.path.stem}_spenderq_ehvo.png"
        fig.suptitle(f"{name}: SpenderQ continuum", fontsize=14)
    return finish_figure(fig, outpath, show)


def plot_flux_over_continuum(epoch, out_dir, show=False):
    """Plot flux / SpenderQ continuum in the observed frame."""
    z_qso = epoch.z_qso
    fig, ax = plt.subplots(figsize=(12, 4), constrained_layout=True)
    ax.plot(epoch.wave, epoch.recon_norm, color="black", lw=0.8)
    ax.axhline(1.0, color="0.5", ls="--", lw=0.8, label="spectrum = continuum")
    ax.axhline(0.0, color="0.8", ls=":", lw=0.6)

    ax.set_xlabel("observed wavelength [Å]")
    ax.set_ylabel("flux / SpenderQ continuum")
    ax.set_ylim(0.0, 2.0)
    ax.set_xlim(RECON_PLOT_XLIM_REST[0] * (1.0 + z_qso), RECON_PLOT_XLIM_REST[1] * (1.0 + z_qso))
    mark_lines(ax, z_qso, label_y=ax.get_ylim()[1] * 0.97)
    ax.legend(fontsize=14)

    outpath = out_dir / f"{epoch.path.stem}_flux_over_continuum_observed.png"
    fig.suptitle(f"{epoch.path.name}: flux / SpenderQ continuum (observed frame)", fontsize=14)
    return finish_figure(fig, outpath, show)
